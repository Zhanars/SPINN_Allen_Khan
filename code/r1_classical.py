"""Classical baselines scored with the same protocol as the PINNs (200-point uniform grid, certified reference):
 - Chebyshev collocation + IMEX Euler (the submitted baseline), several (N, dt)
 - the reaction-only ODE u' = u - u^3 (diffusion neglected), closed form
Writes data_r1/classical.json"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import refsol, r1_core
from chebyshev_spectral import solve_chebyshev

out = {}
for ic in ("compat", "orig"):
    for eps in (0.01, 0.1):
        ref = r1_core.reference(eps, ic, [0.2, 0.5, 0.8, 1.0])
        u0 = r1_core.IC[ic]
        rec = {"cheb": []}
        for N, dt in [(40, 2e-3), (60, 1e-3), (80, 5e-4), (80, 1e-4), (120, 5e-5)]:
            ts = []
            for rep in range(3):
                t0 = time.time(); x, snaps = solve_chebyshev(eps, N=N, dt=dt, u0=u0, save_times=[1.0]); ts.append(time.time() - t0)
            u1 = refsol.bary_interp(x, snaps[1.0], r1_core.XE)
            rec["cheb"].append({"N": N, "dt": dt, "l2": r1_core.l2(u1, ref[1.0]), "linf": r1_core.linf(u1, ref[1.0]), "time_median": float(np.median(ts))})
        # reaction-only ODE
        xe = r1_core.XE; a = u0(xe)
        rec["ode"] = {}
        for t in (0.2, 0.5, 0.8, 1.0):
            u = a * np.exp(t) / np.sqrt(1 + a ** 2 * (np.exp(2 * t) - 1))
            rec["ode"][str(t)] = {"l2": r1_core.l2(u, ref[t]), "linf": r1_core.linf(u, ref[t])}
        out[f"{ic}_eps{eps:g}"] = rec
        print(ic, eps, [(c["N"], c["dt"], f"{c['l2']:.2e}", f"{c['linf']:.2e}", f"{c['time_median']:.2f}s") for c in rec["cheb"]], "ODE t=1:", rec["ode"]["1.0"], flush=True)
json.dump(out, open(os.path.join(HERE, "..", "data_r1", "classical.json"), "w"), indent=1)

# ---- B3: eps = 0.05, T = 5 (error at t = 5)
import time as _time
ref5 = r1_core.reference(0.05, "compat", [1.0, 2.0, 3.0, 4.0, 5.0])
u0 = r1_core.IC["compat"]
rec = {"cheb": []}
for N, dt in [(64, 1e-3), (96, 1e-3), (128, 5e-4), (128, 1e-4)]:
    ts = []
    for rep in range(2):
        t0 = _time.time(); x, snaps = solve_chebyshev(0.05, N=N, dt=dt, T=5.0, u0=u0, save_times=[5.0]); ts.append(_time.time() - t0)
    u5 = refsol.bary_interp(x, snaps[5.0], r1_core.XE)
    rec["cheb"].append({"N": N, "dt": dt, "l2": r1_core.l2(u5, ref5[5.0]), "linf": r1_core.linf(u5, ref5[5.0]), "time_median": float(np.median(ts))})
a = u0(r1_core.XE)
rec["ode"] = {str(t): {"l2": r1_core.l2(a * np.exp(t) / np.sqrt(1 + a ** 2 * (np.exp(2 * t) - 1)), ref5[t])} for t in (1.0, 5.0)}
out["compat_eps0.05_T5"] = rec
print("B3", rec["cheb"], rec["ode"], flush=True)
json.dump(out, open(os.path.join(HERE, "..", "data_r1", "classical.json"), "w"), indent=1)
