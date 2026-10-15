"""Inverse problem for the interface parameter eps with the (revised) Spectral-PINN.

modes
  adam    : joint Adam training of (network, log eps) -- the procedure of the submitted manuscript
  varpro  : eps^2 is, at every epoch, the least-squares minimiser of the PDE residual for the current network
            (the PDE is linear in eps^2: eps^2 = <u_t + f(u), u_xx>_w / <u_xx, u_xx>_w with the quadrature weights w of the PDE loss), relaxed by 0.9/0.1
  phaseA  : N epochs of IC/BC/data fitting with the PDE switched off before the joint phase
Data: certified spectral reference (refsol.reference), observed at n_obs_x equispaced points and the given times,
optionally corrupted by Gaussian noise (absolute sigma)."""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import sys
import time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import refsol
import r1_core
from common import f_ac, chebyshev_vandermonde
from spectral_pinn import SpectralPINN


def run(job):
    eps_true = float(job["eps"]); ic = job.get("ic", "compat"); K = job.get("K", 32); seed = job.get("seed", 42)
    epochs = job.get("epochs", 4000); lr = job.get("lr", 2e-3)
    lam_pde = job.get("lam_pde", 50.0); lam_ic = job.get("lam_ic", 10.0); lam_bc = job.get("lam_bc", 1.0); lam_data = job.get("lam_data", 2.0)
    lam_dataA = job.get("lam_dataA", 1000.0)
    warm = job.get("warmup", 1000); phaseA = job.get("phaseA", 0); mode = job.get("mode", "adam")
    eps_mult = job.get("eps_lr_mult", 10.0); init_f = job.get("eps_init_f", 2.0); n_obs_x = job.get("n_obs_x", 40)
    noise = job.get("noise", 0.0); noise_seed = job.get("noise_seed", 1); T = job.get("T", 1.0)
    obs_times = job.get("obs_times", [0.2, 0.4, 0.6, 0.8, 1.0])
    grid = job.get("grid", "cgl"); hard = job.get("hard_bc", True)

    torch.manual_seed(seed)
    u0 = r1_core.IC[ic]
    x_obs = np.linspace(-1, 1, n_obs_x)
    ref = refsol.reference(eps_true, ic, max(obs_times), [round(t, 6) for t in obs_times], x_obs)
    u_obs = np.stack([np.asarray(ref[round(t, 6)]) for t in obs_times])
    if noise > 0:
        u_obs = u_obs + np.random.default_rng(noise_seed).normal(0, noise, u_obs.shape)
    u_obs_t = torch.tensor(u_obs); Tm_obs = torch.tensor(chebyshev_vandermonde(x_obs, K))

    model = SpectralPINN(eps_true, K=K, seed=seed, u0=u0, grid=grid, hard_bc=hard)
    eps_val = torch.tensor(eps_true * init_f)
    log_eps = torch.tensor(np.log(eps_true * init_f), requires_grad=(mode == "adam"))
    groups = [{"params": model.net.parameters(), "lr": lr}]
    if mode == "adam":
        groups.append({"params": [log_eps], "lr": eps_mult * lr})
    opt = torch.optim.Adam(groups)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    t_col = torch.linspace(1e-4, T, int(100 * T)).reshape(-1, 1)
    t_ic = torch.zeros(1, 1)
    t_obs = torch.tensor([[t] for t in obs_times])
    t_w = torch.linspace(0, T, int(100 * T)).reshape(-1, 1)
    wopt = torch.optim.Adam(model.net.parameters(), lr=1e-3)
    for _ in range(warm):
        wopt.zero_grad(); ((model.coeffs(t_w) - model.c0) ** 2).mean().backward(); wopt.step()

    def terms():
        tc = t_col.clone().requires_grad_(True)
        c = model.coeffs(tc)
        cols = [torch.autograd.grad(c[:, k], tc, grad_outputs=torch.ones_like(c[:, k]), create_graph=True, retain_graph=True)[0] for k in range(K + 1)]
        dcdt = torch.cat(cols, dim=1)
        u = model.field(c); ut = model.field(dcdt); uxx = model.field(c @ model.Dc2.T)
        rp, rm = model.bc_residual(c)
        L_bc = (rp ** 2).mean() + (rm ** 2).mean()
        L_ic = ((model.coeffs(t_ic) - model.c0) ** 2).mean()
        L_data = ((model.coeffs(t_obs) @ Tm_obs.T - u_obs_t) ** 2).mean()
        return ut, uxx, u, L_bc, L_ic, L_data

    t0 = time.time(); traj = []
    if phaseA > 0:
        optA = torch.optim.Adam(model.net.parameters(), lr=1e-3)
        for _ in range(phaseA):
            optA.zero_grad()
            ut, uxx, u, L_bc, L_ic, L_data = terms()
            (lam_dataA * L_data + lam_ic * L_ic + lam_bc * L_bc).backward(); optA.step()
    for ep in range(epochs):
        opt.zero_grad()
        ut, uxx, u, L_bc, L_ic, L_data = terms()
        if mode == "adam":
            e2 = torch.exp(log_eps) ** 2
        else:
            with torch.no_grad():
                wts = model.wq if model.grid == "cgl" else torch.ones(uxx.shape[-1])   # same x-weights as the PDE loss
                e2_ls = torch.clamp(((ut + f_ac(u)) * uxx * wts).sum() / ((uxx ** 2) * wts).sum(), min=1e-8, max=4.0)
                eps_val = torch.sqrt(e2_ls) if ep == 0 else torch.sqrt(0.9 * eps_val ** 2 + 0.1 * e2_ls)
            e2 = eps_val ** 2
        L_pde = model.pde_loss(ut - e2 * uxx + f_ac(u))
        loss = lam_pde * L_pde + lam_ic * L_ic + lam_bc * L_bc + lam_data * L_data
        loss.backward(); opt.step(); sched.step()
        if ep % 50 == 0:
            traj.append(float(torch.exp(log_eps)) if mode == "adam" else float(eps_val))
    eps_hat = float(torch.exp(log_eps)) if mode == "adam" else float(eps_val)
    return {"eps_true": eps_true, "eps_hat": eps_hat, "ratio_to_true": eps_hat / eps_true, "ratio_to_init": eps_hat / (eps_true * init_f),
            "rel_err_pct": 100 * abs(eps_hat - eps_true) / eps_true, "elapsed": time.time() - t0, "L_data": float(L_data), "L_pde": float(L_pde),
            "data_rms_misfit": float(np.sqrt(float(L_data))), "eps_traj_every200": traj[::4]}
