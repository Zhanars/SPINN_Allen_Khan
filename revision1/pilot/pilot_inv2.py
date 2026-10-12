"""Inverse-problem variants. The ablation showed eps_hat ~= 0.717*eps_init for the plain joint Adam scheme (an optimisation artefact).
Variants tested here:
  mode=adam    : original joint training of (theta, log eps)
  mode=varpro  : eps^2 is the closed-form least-squares minimiser of the PDE residual for the current network (variable projection)
  phaseA=N     : N epochs of data/IC/BC fitting (PDE weight 0, eps frozen) before the joint phase
usage: python pilot_inv2.py <ic> <eps_true> <out.json> key=value ..."""
import os, sys, json, time
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code")); sys.path.insert(0, HERE)
import numpy as np
import torch
torch.set_num_threads(1)
import pilot_ref as pr
import spectral_pinn as sp
from common import f_ac, T_FINAL, chebyshev_vandermonde

ic, eps_true, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
kw = dict(a.split("=") for a in sys.argv[4:])
g = lambda k, d, c=float: c(kw.get(k, d))
K = g("K", 32, int); epochs = g("epochs", 4000, int); lr = g("lr", 2e-3); lam_pde = g("lam_pde", 50); lam_ic = g("lam_ic", 10)
lam_bc = g("lam_bc", 1); lam_data = g("lam_data", 2); warm = g("warmup", 1000, int); seed = g("seed", 42, int)
eps_mult = g("eps_lr_mult", 10); eps_init_f = g("eps_init_f", 2.0); n_obs_x = g("n_obs_x", 40, int); noise = g("noise", 0.0)
mode = kw.get("mode", "adam"); lbfgs_n = g("lbfgs", 0, int); eps_b2 = g("eps_beta2", 0.999); phaseA = g("phaseA", 0, int); lam_dataA = g("lam_dataA", 100.0); nreal = g("noise_seed", 1, int)
obs_times = [float(s) for s in kw.get("obs_times", "0.2,0.4,0.6,0.8,1.0").split(",")]

u0 = pr.u0_orig if ic == "orig" else pr.u0_compat
sp.u0_fn = u0
x_obs = np.linspace(-1, 1, n_obs_x)
xn, ref = pr.reference(eps_true, u0, obs_times)
u_obs = np.stack([pr.bary_eval(xn, ref[t], x_obs) for t in obs_times])
if noise > 0:
    u_obs = u_obs + np.random.default_rng(nreal).normal(0, noise, u_obs.shape)
u_obs_t = torch.tensor(u_obs); Tm_obs = torch.tensor(chebyshev_vandermonde(x_obs, K))

torch.manual_seed(seed)
eps_val = torch.tensor(eps_true * eps_init_f)               # current eps (used by varpro)
log_eps = torch.tensor(np.log(eps_true * eps_init_f), requires_grad=(mode == "adam"))
model = sp.SpectralPINN(eps_true, K=K, seed=seed, grid=kw.get("grid", "uniform"), hard_bc=(kw.get("hard", "0") == "1"))
groups = [{"params": model.net.parameters(), "lr": lr}]
if mode == "adam":
    groups.append({"params": [log_eps], "lr": eps_mult * lr, "betas": (0.9, eps_b2)})
opt = torch.optim.Adam(groups)
sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
t_col = torch.linspace(1e-4, T_FINAL, 100).reshape(-1, 1)
t_ic = torch.zeros(1, 1)
t_obs = torch.tensor([[t] for t in obs_times])
t_w = torch.linspace(0, T_FINAL, 100).reshape(-1, 1)
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
if phaseA > 0:  # fit data with the PDE switched off
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
            num = ((ut + f_ac(u)) * uxx).sum(); den = (uxx ** 2).sum()
            e2_ls = torch.clamp(num / den, min=1e-8, max=4.0)
            # relax towards the LS value so that eps does not follow the transient of the network
            eps_val = torch.sqrt(e2_ls) if ep == 0 else torch.sqrt(0.9 * eps_val ** 2 + 0.1 * e2_ls)
        e2 = eps_val ** 2
    L_pde = ((ut - e2 * uxx + f_ac(u)) ** 2).mean()
    loss = lam_pde * L_pde + lam_ic * L_ic + lam_bc * L_bc + lam_data * L_data
    loss.backward(); opt.step(); sched.step()
    if ep % 50 == 0:
        traj.append(float(torch.exp(log_eps)) if mode == "adam" else float(eps_val))
if lbfgs_n > 0:
    if mode != "adam":
        raise SystemExit("lbfgs only with mode=adam")
    params = list(model.net.parameters()) + [log_eps]
    for p in params: p.requires_grad_(True)
    optL = torch.optim.LBFGS(params, lr=1.0, max_iter=lbfgs_n, max_eval=int(1.3 * lbfgs_n), tolerance_grad=1e-14, tolerance_change=1e-16, history_size=50, line_search_fn="strong_wolfe")
    def closure():
        optL.zero_grad()
        ut, uxx, u, L_bc, L_ic, L_data = terms()
        e2 = torch.exp(log_eps) ** 2
        L_pde = ((ut - e2 * uxx + f_ac(u)) ** 2).mean()
        loss = lam_pde * L_pde + lam_ic * L_ic + lam_bc * L_bc + lam_data * L_data
        loss.backward()
        return loss
    optL.step(closure)
    ut, uxx, u, L_bc, L_ic, L_data = terms()
    L_pde = ((ut - torch.exp(log_eps).detach() ** 2 * uxx + f_ac(u)) ** 2).mean()
eps_hat = float(torch.exp(log_eps)) if mode == "adam" else float(eps_val)
res = {"ic": ic, "eps_true": eps_true, "eps_hat": eps_hat, "ratio": eps_hat / eps_true, "rel_err_pct": 100 * abs(eps_hat - eps_true) / eps_true,
       "elapsed": time.time() - t0, "L_data": float(L_data), "L_pde": float(L_pde), "traj": traj[::4], "cfg": kw}
json.dump(res, open(out, "w"), indent=1)
print(ic, eps_true, kw, "eps_hat=%.5f ratio=%.3f L_data=%.2e L_pde=%.2e" % (eps_hat, eps_hat / eps_true, float(L_data), float(L_pde)))
