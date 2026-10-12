"""Checks supporting the code audit (part e).  Writes code_checks.json.
 1. coefficient-space operators of the repo (Dc, Dc2, BC vectors) vs numpy.polynomial.chebyshev
 2. BC-loss scale: L_bc of the IC fit c0 (u_x = -+2) vs K, relative to the PDE-loss scale
 3. wall slopes u_x(+-1,t) of the stored K=32 Spectral-PINN predictions (data/spectral_pinn.json)
 4. PDE residual of the repo's spectral operators evaluated on the trusted solution (sign/operator check)
 5. oracle test of the inverse-problem PDE-residual estimator of eps (no training): plug the K-projected
    trusted solution into the least-squares eps estimator
"""
import os, sys, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, r"D:\Videos\KazNPU\Symmetry_PINN_2026\code")
import refsol
import common  # repo (imports torch, sets float64)
from numpy.polynomial import chebyshev as C
DATA = r"D:\Videos\KazNPU\Symmetry_PINN_2026\data"
out = {}

# 1. operators
res = {}
for K in (8, 32):
    Dc = common.cheb_coeff_derivative_matrix(K)
    rng = np.random.default_rng(0)
    a = rng.standard_normal(K + 1)
    d1 = C.chebder(a)
    d1 = np.concatenate([d1, np.zeros(K + 1 - len(d1))])
    d2 = C.chebder(a, 2)
    d2 = np.concatenate([d2, np.zeros(K + 1 - len(d2))])
    res[f"K{K}_D_err"] = float(np.abs(Dc @ a - d1).max())
    res[f"K{K}_D2_err"] = float(np.abs(Dc @ (Dc @ a) - d2).max())
    vp, vm = common.cheb_boundary_deriv_vectors(K)
    res[f"K{K}_bc_err"] = float(max(abs(vp @ a - C.chebval(1.0, C.chebder(a))), abs(vm @ a - C.chebval(-1.0, C.chebder(a)))))
out["operators"] = res
print("operators:", res)

# 2. BC loss scale of the IC fit
x300 = np.linspace(-1, 1, 300)
bc = {}
for K in (16, 32, 64, 128, 256):
    c0 = common.fit_chebyshev_coeffs(common.u0_fn(x300), x300, K)
    vp, vm = common.cheb_boundary_deriv_vectors(K)
    s = max(1.0, K**2)
    r_plus, r_minus = c0 @ vp, c0 @ vm
    bc[K] = dict(ux_plus=float(r_plus), ux_minus=float(r_minus), L_bc_repo=float((r_plus / s)**2 + (r_minus / s)**2),
                 L_bc_unscaled=float(r_plus**2 + r_minus**2), scale_factor_K4=float(K**4))
    print("K", K, bc[K])
out["bc_scale"] = bc

# 3. wall slopes of stored predictions
xu = np.linspace(-1, 1, 200)
ws = {}
for name in ("spectral_pinn", "spectral_pinn_energy"):
    d = json.load(open(os.path.join(DATA, name + ".json")))
    ws[name] = {}
    for t in ("0.2", "0.5", "0.8", "1.0"):
        up = np.array(d["errors"][t]["u_pred"])
        c = C.chebfit(xu, up, 32)         # exact: the prediction is a degree-32 polynomial
        resid = np.abs(C.chebval(xu, c) - up).max()
        dc = C.chebder(c)
        ws[name][t] = dict(ux_plus1=float(C.chebval(1.0, dc)), ux_minus1=float(C.chebval(-1.0, dc)), fit_resid=float(resid))
    print(name, ws[name])
# trusted wall slope for t>0 is 0 by construction; slope of trusted K=32 projection (for comparison)
tr, xn, nodes = refsol.reference(0.01, "orig", 1.0, [0.2, 0.5, 0.8, 1.0], xu, N=256, return_nodes=True)
ws["trusted_K32_projection"] = {}
x300 = np.linspace(-1, 1, 300)
for t in (0.2, 0.5, 0.8, 1.0):
    u300 = refsol.bary_interp(xn, nodes[t], x300)
    c = C.chebfit(x300, u300, 32)
    dc = C.chebder(c)
    ws["trusted_K32_projection"][str(t)] = dict(ux_plus1=float(C.chebval(1.0, dc)), ux_minus1=float(C.chebval(-1.0, dc)))
print("trusted K32 projection", ws["trusted_K32_projection"])
out["wall_slopes"] = ws

# 4 + 5. operator check and oracle estimator of eps
orc = {}
for eps_true in (0.005, 0.01, 0.02, 0.05, 0.1):
    t_col = np.linspace(1e-4, 1.0, 100)
    h = 1e-5
    N = 512 if eps_true < 0.02 else 256
    _, xn, nodes = refsol.reference(eps_true, "orig", 1.0, list(t_col), x300, N=N, return_nodes=True)
    _, Dn = refsol.cheb_nodes_diff(N)
    D2n = Dn @ Dn
    for K in (16, 32, 64, 96):
        Dc = common.cheb_coeff_derivative_matrix(K)
        Dc2 = Dc @ Dc
        Tm = common.chebyshev_vandermonde(x300, K)
        num = den = 0.0
        res_true = 0.0
        resid_rows = []
        for t in t_col:
            u = nodes[t]
            ut = eps_true**2 * (D2n @ u) - (u**3 - u)           # exact PDE time derivative on nodes (t>0)
            u300 = refsol.bary_interp(xn, u, x300)
            ut300 = refsol.bary_interp(xn, ut, x300)
            c = np.linalg.lstsq(Tm, u300, rcond=None)[0]
            dcdt = np.linalg.lstsq(Tm, ut300, rcond=None)[0]
            uK, utK, uxxK = Tm @ c, Tm @ dcdt, Tm @ (Dc2 @ c)
            g = utK + (uK**3 - uK)                  # = eps^2 u_xx  (if exact)
            num += np.sum(g * uxxK)
            den += np.sum(uxxK**2)
            resid_rows.append(np.mean((g - eps_true**2 * uxxK)**2))
        eps2_hat = num / den
        orc[f"eps{eps_true}_K{K}"] = dict(eps_hat_over_eps_true=float(np.sqrt(max(eps2_hat, 0)) / eps_true),
                                          mean_sq_residual_at_true_eps=float(np.mean(resid_rows)))
    print("oracle eps_true", eps_true, {K: round(orc[f"eps{eps_true}_K{K}"]["eps_hat_over_eps_true"], 4) for K in (16, 32, 64, 96)},
          {K: f"{orc[f'eps{eps_true}_K{K}']['mean_sq_residual_at_true_eps']:.2e}" for K in (16, 32, 64, 96)}, flush=True)
out["oracle_eps_estimator"] = orc
json.dump(out, open(os.path.join(HERE, "code_checks.json"), "w"), indent=1)
