"""Round-1 revision harness: unified training/evaluation of the 1D Allen-Cahn experiments.

* certified reference solutions (refsol.py: Chebyshev collocation + Radau, N-doubling certified),
* all errors on the same 200-point uniform grid (dx = 2/199), Linf on the same grid,
* maximum-bound-principle check (max |u| over a 1001 x (100 T + 1) space-time grid),
* energy monotonicity counted at several time resolutions, energy by Clenshaw-Curtis (exact for the polynomial integrand),
* a-posteriori energy balance diagnostics for the trained trajectory:
      dE/dt = -||u_t||^2 + <r, u_t> + eps^2 [u_x u_t]_{x=-1}^{x=1},   r = u_t - eps^2 u_xx + f(u)
  (identity checked against dE/dt = int (eps^2 u_x u_xt + f(u) u_t) dx, evaluated by the chain rule with dc/dt from autograd, and quantified through the integrated residual).
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import sys
import json
import time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import refsol
from common import (f_ac, F_ac, T_FINAL, chebyshev_vandermonde, cheb_coeff_derivative_matrix,
                    cheb_boundary_deriv_vectors)

torch.set_default_dtype(torch.float64)

XE = np.linspace(-1.0, 1.0, 200)
DXE = XE[1] - XE[0]
IC = {"orig": refsol.u0_orig, "compat": refsol.u0_compat}
REF_DIR = os.path.join(HERE, "..", "data_r1", "ref")
os.makedirs(REF_DIR, exist_ok=True)


def l2(a, b):
    return float(np.sqrt(np.sum((a - b) ** 2) * DXE))


def linf(a, b):
    return float(np.max(np.abs(a - b)))


def reference(eps, ic, times):
    """dict t -> u_ref on XE (certified spectral reference, cached on disk)."""
    times = [round(float(t), 6) for t in times]
    path = os.path.join(REF_DIR, f"{ic}_eps{eps:g}_T{max(times):g}_n{len(times)}.npz")
    if os.path.exists(path):
        z = np.load(path)
        return {t: z[f"t{i}"] for i, t in enumerate(times)}
    out = refsol.reference(eps, ic, max(times), times, XE)
    np.savez(path, **{f"t{i}": np.asarray(out[t]) for i, t in enumerate(times)})
    return {t: np.asarray(out[t]) for t in times}


# --------------------------------------------------------------------------------------------------
# quadrature / operators shared by the diagnostics
# --------------------------------------------------------------------------------------------------
class Quad:
    """Clenshaw-Curtis nodes of order Nq=4K+8 with Chebyshev Vandermonde and derivative blocks."""
    def __init__(self, K):
        self.K = K
        self.N = 4 * K + 8
        j = np.arange(self.N + 1)
        self.x = np.cos(np.pi * j / self.N)
        self.w = torch.tensor(refsol.clenshaw_curtis_weights(self.N))
        self.T = torch.tensor(chebyshev_vandermonde(self.x, K))                  # (Nq+1, K+1)
        Dc = cheb_coeff_derivative_matrix(K)
        self.Dc = torch.tensor(Dc)
        self.Dc2 = torch.tensor(Dc @ Dc)
        vp, vm = cheb_boundary_deriv_vectors(K)
        self.vp = torch.tensor(vp); self.vm = torch.tensor(vm)                    # u_x(+1), u_x(-1) from coefficients
        self.Tp = torch.tensor(chebyshev_vandermonde(np.array([1.0]), K))[0]      # u(+1)
        self.Tm = torch.tensor(chebyshev_vandermonde(np.array([-1.0]), K))[0]     # u(-1)


def energy_cc(q, c, eps):
    u = c @ q.T.t()
    ux = (c @ q.Dc.t()) @ q.T.t()
    return ((0.5 * eps ** 2 * ux ** 2 + F_ac(u)) * q.w).sum(dim=-1)


def count_increases(E, tol=1e-10):
    d = np.diff(np.asarray(E))
    return int(np.sum(d > tol)), float(d.max())


# --------------------------------------------------------------------------------------------------
# evaluation of a trained Spectral-PINN
# --------------------------------------------------------------------------------------------------
def eval_spectral(model, eps, ic, T=1.0, err_times=(0.2, 0.5, 0.8, 1.0), diag=True):
    K = model.K
    q = Quad(K)
    res = {"errors": {}}
    ref = reference(eps, ic, sorted(set(list(err_times) + [T])))
    Tm_e = torch.tensor(chebyshev_vandermonde(XE, K))
    with torch.no_grad():
        for t in err_times:
            if t > T + 1e-12:
                continue
            c = model.coeffs(torch.tensor([[float(t)]]))
            u = (c @ Tm_e.t()).numpy().ravel()
            res["errors"][str(t)] = {"l2": l2(u, ref[round(t, 6)]), "linf": linf(u, ref[round(t, 6)])}
        # maximum bound principle: max |u| on a fine space-time grid
        tg = torch.linspace(0, T, int(100 * T) + 1).reshape(-1, 1)
        c = model.coeffs(tg)
        xf = np.linspace(-1, 1, 1001)
        uf = (c @ torch.tensor(chebyshev_vandermonde(xf, K)).t()).numpy()
        res["max_abs_u"] = float(np.abs(uf).max())
        res["mbp_excess"] = float(max(0.0, np.abs(uf).max() - 1.0))
        # energy monotonicity at several resolutions
        res["energy"] = {}
        for tag, per_unit in (("dt0.02", 50), ("dt0.002", 500), ("dt0.0002", 5000)):
            nn = int(round(per_unit * T)) + 1
            tt = torch.linspace(0, T, nn).reshape(-1, 1)
            E = energy_cc(q, model.coeffs(tt), eps).numpy()
            cnt, mx = count_increases(E)
            res["energy"][tag] = {"n_points": nn, "n_increases": cnt, "max_step": mx, "E_first": float(E[0]), "E_last": float(E[-1])}
            if tag == "dt0.02":
                res["E_traj_dt0.02"] = E.tolist()
        # relative decrease margin: smallest decrement / typical decrement
    if diag:
        res["balance"] = balance_diagnostics(model, q, eps, T)
    return res


def balance_diagnostics(model, q, eps, T, nt=201):
    """Residual, boundary term and identity check along the trained trajectory (needs d c / d t by autograd)."""
    K = model.K
    t = torch.linspace(1e-6, T, nt).reshape(-1, 1).requires_grad_(True)
    c = model.coeffs(t)
    cols = []
    for k in range(K + 1):
        cols.append(torch.autograd.grad(c[:, k], t, grad_outputs=torch.ones_like(c[:, k]), retain_graph=True, create_graph=False)[0])
    cd = torch.cat(cols, dim=1).detach()
    c = c.detach()
    Tq = q.T.t()
    u = c @ Tq; ut = cd @ Tq
    ux = (c @ q.Dc.t()) @ Tq; uxx = (c @ q.Dc2.t()) @ Tq
    uxt = (cd @ q.Dc.t()) @ Tq
    r = ut - eps ** 2 * uxx + f_ac(u)
    w = q.w
    n_ut2 = (ut ** 2 * w).sum(-1); n_r2 = (r ** 2 * w).sum(-1); ip = (r * ut * w).sum(-1)
    # boundary term eps^2 [u_x u_t]_{-1}^{1}: u_x(+-1) from coefficients, u_t(+-1) from time derivative of coefficients
    ux_p = c @ q.vp; ux_m = c @ q.vm
    ut_p = cd @ q.Tp; ut_m = cd @ q.Tm
    bterm = eps ** 2 * (ux_p * ut_p - ux_m * ut_m)
    # dE/dt directly:  sum w (eps^2 u_x u_xt + f(u) u_t)
    dEdt = ((eps ** 2 * ux * uxt + f_ac(u) * ut) * w).sum(-1)
    rhs = -n_ut2 + ip + bterm
    tt = t.detach().numpy().ravel()
    scale = float(np.max(np.abs(dEdt.numpy())) + 1e-30)
    out = {"identity_max_abs_err": float(np.max(np.abs(dEdt.numpy() - rhs.numpy()))), "identity_scale": scale,
           "identity_rel_err": float(np.max(np.abs(dEdt.numpy() - rhs.numpy())) / scale)}
    out["int_half_r2"] = float(np.trapezoid(0.5 * n_r2.numpy(), tt))
    out["int_abs_bterm"] = float(np.trapezoid(np.abs(bterm.numpy()), tt))
    out["int_n_ut2"] = float(np.trapezoid(n_ut2.numpy(), tt))
    out["rms_residual_L2"] = float(np.sqrt(np.trapezoid(n_r2.numpy(), tt) / T))
    out["max_residual_L2"] = float(np.sqrt(n_r2.numpy().max()))
    out["max_abs_ux_wall"] = float(max(np.abs(ux_p.numpy()).max(), np.abs(ux_m.numpy()).max()))
    out["max_dEdt"] = float(dEdt.numpy().max())
    out["bound_E_increase"] = out["int_half_r2"] + out["int_abs_bterm"]
    E0 = float(energy_cc(q, c[:1], eps)[0]); E1 = float(energy_cc(q, c[-1:], eps)[0])
    out["E_first"] = E0; out["E_last"] = E1
    out["t_series"] = {"t": tt[::10].tolist(), "dEdt": dEdt.numpy()[::10].tolist(), "rhs": rhs.numpy()[::10].tolist(),
                       "n_ut2": n_ut2.numpy()[::10].tolist(), "n_r": np.sqrt(n_r2.numpy())[::10].tolist(), "bterm": bterm.numpy()[::10].tolist()}
    return out


# --------------------------------------------------------------------------------------------------
# vanilla PINN evaluation
# --------------------------------------------------------------------------------------------------
def eval_vanilla(net, eps, ic, T=1.0, err_times=(0.2, 0.5, 0.8, 1.0)):
    from vanilla_pinn import eval_on_grid
    ref = reference(eps, ic, sorted(set(list(err_times) + [T])))
    res = {"errors": {}}
    for t in err_times:
        u = eval_on_grid(net, XE, t)
        res["errors"][str(t)] = {"l2": l2(u, ref[round(t, 6)]), "linf": linf(u, ref[round(t, 6)])}
    xf = np.linspace(-1, 1, 1001)
    mx = 0.0
    for t in np.linspace(0, T, int(100 * T) + 1):
        mx = max(mx, float(np.abs(eval_on_grid(net, xf, float(t))).max()))
    res["max_abs_u"] = mx; res["mbp_excess"] = max(0.0, mx - 1.0)
    # energy along the trajectory (trapezoid on 1001 points, u_x by autograd)
    E = []
    tg = np.linspace(0, T, int(50 * T) + 1)
    x = torch.tensor(xf).reshape(-1, 1).requires_grad_(True)
    w = np.full(len(xf), xf[1] - xf[0]); w[0] *= 0.5; w[-1] *= 0.5
    for t in tg:
        xt = torch.cat([x, torch.full_like(x, float(t))], dim=1)
        u = net(xt)
        ux = torch.autograd.grad(u, x, grad_outputs=torch.ones_like(u))[0]
        e = (0.5 * eps ** 2 * ux ** 2 + F_ac(u)).detach().numpy().ravel()
        E.append(float((e * w).sum()))
    cnt, mx_step = count_increases(E)
    res["energy"] = {"dt0.02": {"n_points": len(E), "n_increases": cnt, "max_step": mx_step, "E_first": E[0], "E_last": E[-1]}}
    res["E_traj_dt0.02"] = E
    return res
