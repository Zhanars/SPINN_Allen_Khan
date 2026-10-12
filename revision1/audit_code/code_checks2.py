"""More checks for the code audit (part e).  Writes code_checks2.json.
 A. conditioning of the equispaced 300-point Chebyshev Vandermonde used for the IC fit and for the
    PDE collocation in x (spectral_pinn.py L65, L68, L80) as a function of K; IC-fit sup error between
    the grid points; coefficient size.
 B. location of the Spectral-PINN error (seed 42, eps=0.01, t=1): wall vs interior.
 C. statistics of the inverse-bias sweep (eps_hat/eps_true).
 D. data-fit/BC-conflict: IC representation error.
"""
import os, sys, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
from numpy.polynomial import chebyshev as C
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import refsol
DATA = r"D:\Videos\KazNPU\Symmetry_PINN_2026\data"
out = {}

# A. conditioning of the equispaced fit
x300 = np.linspace(-1, 1, 300)
xf = np.linspace(-1, 1, 40001)
A = {}
for K in (16, 32, 48, 64, 96, 128, 192, 256):
    Tm = C.chebvander(x300, K)
    sv = np.linalg.svd(Tm, compute_uv=False)
    c = np.linalg.lstsq(Tm, refsol.u0_orig(x300), rcond=None)[0]      # same as common.fit_chebyshev_coeffs (rcond=None)
    ufit = C.chebval(xf, c)
    e = ufit - refsol.u0_orig(xf)
    dc = C.chebder(c)
    A[K] = dict(cond_T300=float(sv[0] / sv[-1]), max_abs_coeff=float(np.abs(c).max()),
                IC_fit_Linf_fine=float(np.abs(e).max()), IC_fit_Linf_on_300=float(np.abs(C.chebval(x300, c) - refsol.u0_orig(x300)).max()),
                IC_fit_wall_slope=float(C.chebval(1.0, dc)))
    print("K", K, A[K], flush=True)
out["equispaced_fit"] = A
out["rule_of_thumb"] = "stable polynomial fitting on n equispaced points needs degree <~ 2*sqrt(n): n=300 -> ~35"

# B. error location (seed 42 stored predictions)
xu = np.linspace(-1, 1, 200)
d = json.load(open(os.path.join(DATA, "spectral_pinn.json")))
tr = refsol.reference(0.01, "orig", 1.0, [0.2, 0.5, 0.8, 1.0], xu, N=256)
B = {}
for t in (0.2, 0.5, 0.8, 1.0):
    up = np.array(d["errors"][str(t)]["u_pred"])
    e = up - tr[t]
    ae = np.abs(e)
    wall = np.abs(xu) > 0.95
    B[str(t)] = dict(Linf=float(ae.max()), argmax_x=float(xu[ae.argmax()]), Linf_wall_gt0p95=float(ae[wall].max()),
                     Linf_interior_le0p95=float(ae[~wall].max()),
                     L2_interior_le0p95=float(np.sqrt(np.sum(e[~wall]**2) * (xu[1] - xu[0]))),
                     L2_wall_gt0p95=float(np.sqrt(np.sum(e[wall]**2) * (xu[1] - xu[0]))))
    print("t", t, B[str(t)])
out["spectral_pinn_error_location"] = B
dv = json.load(open(os.path.join(DATA, "vanilla_pinn.json")))
Bv = {}
for t in (0.2, 0.5, 0.8, 1.0):
    up = np.array(dv["errors"][str(t)]["u_pred"])
    e = up - tr[t]
    ae = np.abs(e)
    Bv[str(t)] = dict(Linf=float(ae.max()), argmax_x=float(xu[ae.argmax()]), rel_L2=float(np.sqrt(np.sum(e**2)) / np.sqrt(np.sum(tr[t]**2))),
                      L2_trusted=float(np.sqrt(np.sum(tr[t]**2) * (xu[1] - xu[0]))),
                      L2_of_trivial_zero_prediction=float(np.sqrt(np.sum(tr[t]**2) * (xu[1] - xu[0]))))
    print("vanilla t", t, Bv[str(t)])
out["vanilla_error_location"] = Bv

# C. inverse sweep statistics
s = json.load(open(os.path.join(DATA, "inverse_bias_sweep.json")))["results"]
rat = {}
for r in s:
    rat.setdefault(r["K"], []).append(r["eps_hat"] / r["eps_true"])
C_ = {str(K): dict(mean=float(np.mean(v)), sd=float(np.std(v, ddof=1)), min=float(np.min(v)), max=float(np.max(v)),
                   mean_sq=float(np.mean(np.square(v)))) for K, v in rat.items()}
allr = [r["eps_hat"] / r["eps_true"] for r in s if r["K"] >= 32]
C_["K>=32 pooled"] = dict(mean=float(np.mean(allr)), sd=float(np.std(allr, ddof=1)), n=len(allr), sqrt2=float(np.sqrt(2)),
                          mean_sq=float(np.mean(np.square(allr))))
print("inverse sweep ratio stats", C_)
out["inverse_sweep_ratio_stats"] = C_
json.dump(out, open(os.path.join(HERE, "code_checks2.json"), "w"), indent=1)
