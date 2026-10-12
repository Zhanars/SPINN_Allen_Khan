"""Pilot for the inverse problem: is the ~38% eps bias caused by (a) the Neumann-incompatible IC, (b) the FD reference used to
generate the observations, or (c) something in the loss/optimisation?
usage: python pilot_inv.py <ic orig|compat> <data fd|spec> <eps_true> <out.json> [key=value ...]"""
import os, sys, json, time
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))
sys.path.insert(0, HERE)
import numpy as np
import torch
torch.set_num_threads(1)
import pilot_ref as pr
import common
import spectral_pinn as sp
from common import f_ac, T_FINAL, chebyshev_vandermonde

ic, data, eps_true, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
kw = dict(a.split("=") for a in sys.argv[5:])
K = int(kw.get("K", 32)); epochs = int(kw.get("epochs", 4000)); lr = float(kw.get("lr", 2e-3))
lam_pde = float(kw.get("lam_pde", 50)); lam_ic = float(kw.get("lam_ic", 10)); lam_bc = float(kw.get("lam_bc", 1)); lam_data = float(kw.get("lam_data", 2))
warm = int(kw.get("warmup", 1000)); seed = int(kw.get("seed", 42)); eps_mult = float(kw.get("eps_lr_mult", 10)); eps_init_f = float(kw.get("eps_init_f", 2.0))
n_obs_x = int(kw.get("n_obs_x", 40)); obs_times = [float(s) for s in kw.get("obs_times", "0.2,0.4,0.6,0.8,1.0").split(",")]
noise = float(kw.get("noise", 0.0))

u0 = pr.u0_orig if ic == "orig" else pr.u0_compat
sp.u0_fn = u0
x_obs = np.linspace(-1, 1, n_obs_x)
if data == "fd":
    assert ic == "orig"
    x_ref, snaps = common.reference_solution_fd(eps_true, save_times=obs_times)
    u_obs = np.stack([np.interp(x_obs, x_ref, snaps[t]) for t in obs_times])
else:
    xn, ref = pr.reference(eps_true, u0, obs_times)
    u_obs = np.stack([pr.bary_eval(xn, ref[t], x_obs) for t in obs_times])
if noise > 0:
    u_obs = u_obs + np.random.default_rng(1).normal(0, noise, u_obs.shape)
u_obs_t = torch.tensor(u_obs)
Tm_obs = torch.tensor(chebyshev_vandermonde(x_obs, K))

torch.manual_seed(seed)
log_eps = torch.tensor(np.log(eps_true * eps_init_f), requires_grad=True)
model = sp.SpectralPINN(eps_true, K=K, seed=seed)
opt = torch.optim.Adam([{"params": model.net.parameters(), "lr": lr}, {"params": [log_eps], "lr": eps_mult * lr}])
sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
t_col = torch.linspace(1e-4, T_FINAL, 100).reshape(-1, 1)
t_ic = torch.zeros(1, 1)
t_obs = torch.tensor([[t] for t in obs_times])
wopt = torch.optim.Adam(model.net.parameters(), lr=1e-3)
t_w = torch.linspace(0, T_FINAL, 100).reshape(-1, 1)
for _ in range(warm):
    wopt.zero_grad(); ((model.coeffs(t_w) - model.c0) ** 2).mean().backward(); wopt.step()
t0 = time.time(); hist = []
for ep in range(epochs):
    opt.zero_grad()
    e = torch.exp(log_eps)
    tc = t_col.clone().requires_grad_(True)
    c = model.coeffs(tc)
    cols = [torch.autograd.grad(c[:, k], tc, grad_outputs=torch.ones_like(c[:, k]), create_graph=True, retain_graph=True)[0] for k in range(K + 1)]
    dcdt = torch.cat(cols, dim=1)
    u = model.field(c); ut = model.field(dcdt); uxx = model.field(c @ model.Dc2.T)
    L_pde = ((ut - e ** 2 * uxx + f_ac(u)) ** 2).mean()
    L_ic = ((model.coeffs(t_ic) - model.c0) ** 2).mean()
    rp, rm = model.bc_residual(c)
    L_bc = (rp ** 2).mean() + (rm ** 2).mean()
    L_data = ((model.coeffs(t_obs) @ Tm_obs.T - u_obs_t) ** 2).mean()
    loss = lam_pde * L_pde + lam_ic * L_ic + lam_bc * L_bc + lam_data * L_data
    loss.backward(); opt.step(); sched.step()
    if ep % 50 == 0:
        hist.append(float(e))
eps_hat = float(torch.exp(log_eps))
res = {"ic": ic, "data": data, "eps_true": eps_true, "eps_hat": eps_hat, "rel_err_pct": 100 * abs(eps_hat - eps_true) / eps_true,
       "ratio": eps_hat / eps_true, "elapsed": time.time() - t0, "L_data": float(L_data), "L_pde": float(L_pde), "L_ic": float(L_ic), "L_bc": float(L_bc),
       "eps_traj": hist[::4], "cfg": kw}
json.dump(res, open(out, "w"), indent=1)
print(ic, data, eps_true, "eps_hat=%.5f ratio=%.3f  L_data=%.2e L_pde=%.2e" % (eps_hat, eps_hat / eps_true, float(L_data), float(L_pde)))
