"""Classical baselines that complement data_r1/classical.json (same protocol as the PINN runs: 200-point grid, certified reference).

usage:  python r1_classical2.py b3        -> Chebyshev IMEX solver and reaction-only ODE on B3 (eps=0.05, T=5)
        python r1_classical2.py inverse   -> classical least-squares estimator of eps + identifiability (Fisher) numbers
        python r1_classical2.py floor     -> polynomial representation floor of the reference (degree K, unconstrained / Neumann basis)
        python r1_classical2.py timing    -> wall time of the Chebyshev IMEX solver on B1/B2 on this machine (run on an idle machine)
Each mode writes data_r1/classical2_<mode>.json.

Inverse: the data are generated exactly as in r1_inverse.py (certified reference at the 40 x 5 observation points, Gaussian noise
np.random.default_rng(noise_seed).normal(0, sigma, (5, 40))), so noise seeds 1..5 are the same realisations as the neural runs
(paired comparison); seeds 101..120 are extra replicates for the statistics.  Forward model: the Chebyshev/Radau solver of refsol.py
at a lower resolution (N=96, rtol 1e-10).  Estimator: nonlinear least squares in log(eps) started from 2 eps_true (the start of the
first version)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
import sys
import json
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import refsol

OUT_DIR = os.path.join(HERE, "..", "data_r1")
OBS_T = [0.2, 0.4, 0.6, 0.8, 1.0]
X_OBS = np.linspace(-1, 1, 40)


def save(name, d):
    json.dump(d, open(os.path.join(OUT_DIR, f"classical2_{name}.json"), "w"), indent=1)


def mode_b3():
    import r1_core
    from chebyshev_spectral import solve_chebyshev
    eps, T = 0.05, 5.0
    times = [0.2, 0.5, 0.8, 1.0, 2.0, 3.0, 5.0]          # same list (and cache file) as the B3 PINN runs
    ref = r1_core.reference(eps, "compat", times)
    rec = {"cheb": []}
    for N, dt in [(64, 5e-4), (80, 5e-4), (80, 1e-4), (128, 5e-4), (128, 1e-4)]:
        ts = []
        for rep in range(3):
            t0 = time.time(); x, snaps = solve_chebyshev(eps, N=N, dt=dt, T=T, u0=refsol.u0_compat, save_times=[1.0, T]); ts.append(time.time() - t0)
        e = {}
        for t in (1.0, T):
            u = refsol.bary_interp(x, snaps[t], r1_core.XE)
            e[str(t)] = {"l2": r1_core.l2(u, ref[t]), "linf": r1_core.linf(u, ref[t])}
        rec["cheb"].append({"N": N, "dt": dt, "errors": e, "time_median": float(np.median(ts))})
        print("B3 cheb", N, dt, e, "%.2fs" % np.median(ts), flush=True)
    a = refsol.u0_compat(r1_core.XE)
    rec["ode"] = {}
    for t in times:
        u = a * np.exp(t) / np.sqrt(1 + a ** 2 * (np.exp(2 * t) - 1))
        rec["ode"][str(t)] = {"l2": r1_core.l2(u, ref[t]), "linf": r1_core.linf(u, ref[t])}
    save("b3", rec)


def data_clean(eps_true):
    ref = refsol.reference(eps_true, "compat", max(OBS_T), [round(t, 6) for t in OBS_T], X_OBS)
    return np.stack([np.asarray(ref[round(t, 6)]) for t in OBS_T])


def model(eps, N=96, rtol=1e-10, atol=1e-12):
    ref = refsol.reference(float(eps), "compat", max(OBS_T), OBS_T, X_OBS, N=N, rtol=rtol, atol=atol)
    return np.stack([np.asarray(ref[t]) for t in OBS_T])


def mode_inverse():
    from scipy.optimize import least_squares
    out = {}
    for eps_true in (0.01, 0.05, 0.1):
        clean = data_clean(eps_true)
        h = 1e-3
        dlog = (model(eps_true * np.exp(h), rtol=1e-11, atol=1e-13) - model(eps_true * np.exp(-h), rtol=1e-11, atol=1e-13)) / (2 * h)
        S = float(np.linalg.norm(dlog))
        rec = {"S_norm_dlogeps": S, "S_rms": float(np.sqrt(np.mean(dlog ** 2))), "model_misfit_at_true": float(np.abs(model(eps_true) - clean).max()), "runs": {}}
        for sigma in (0.0, 0.01, 0.05):
            seeds = [0] if sigma == 0 else list(range(1, 6)) + list(range(101, 121))
            est = {}
            t0 = time.time()
            for ns in seeds:
                data = clean + (np.random.default_rng(ns).normal(0, sigma, clean.shape) if sigma > 0 else 0.0)
                fun = lambda p: (model(np.exp(p[0])) - data).ravel()
                r = least_squares(fun, [np.log(2 * eps_true)], method="trf", xtol=1e-12, ftol=1e-12, gtol=1e-12, bounds=([np.log(1e-4)], [np.log(2.0)]))
                est[str(ns)] = float(np.exp(r.x[0]))
            v = np.array(list(est.values()))
            rel = 100 * np.abs(v - eps_true) / eps_true
            rec["runs"][str(sigma)] = {"eps_hat": est, "rel_err_pct_mean": float(rel.mean()), "rel_rmse_pct": float(100 * np.sqrt(np.mean((v - eps_true) ** 2)) / eps_true),
                                       "rel_bias_pct": float(100 * (v.mean() - eps_true) / eps_true), "n": int(len(v)),
                                       "fisher_rel_se_pct": float(100 * sigma / S) if sigma > 0 else 0.0, "elapsed": time.time() - t0}
            print("inverse", eps_true, sigma, {k: round(x, 4) for k, x in rec["runs"][str(sigma)].items() if not isinstance(x, dict)}, flush=True)
        out[str(eps_true)] = rec
    save("inverse", out)


def mode_timing():
    from chebyshev_spectral import solve_chebyshev
    rec = {}
    for eps in (0.01, 0.1):
        for N, dt in [(80, 5e-4), (80, 1e-4)]:
            ts = []
            for rep in range(5):
                t0 = time.time(); solve_chebyshev(eps, N=N, dt=dt, T=1.0, u0=refsol.u0_compat, save_times=[1.0]); ts.append(time.time() - t0)
            rec[f"eps{eps:g}_N{N}_dt{dt:g}"] = {"times": ts, "median": float(np.median(ts))}
            print("timing", eps, N, dt, "%.3fs" % np.median(ts), flush=True)
    save("timing", rec)


def mode_floor():
    """Representation floor: smallest discrete L2 error (200-point grid) of a polynomial of degree <= K, unconstrained and in the exact Neumann
    basis phi_k = T_k - k^2/(k+2)^2 T_{k+2}, fitted by least squares to the certified reference on the evaluation grid."""
    xe = np.linspace(-1, 1, 200); dx = xe[1] - xe[0]
    rec = {}
    for name, eps, t in (("b1", 0.01, 1.0), ("b2", 0.1, 1.0), ("b3", 0.05, 5.0)):
        T = 5.0 if name == "b3" else 1.0
        ref = refsol.reference(eps, "compat", T, [t], xe)[t]
        for K in (16, 32, 64, 128):
            V = np.polynomial.chebyshev.chebvander(xe, K)
            c, *_ = np.linalg.lstsq(V, ref, rcond=None)
            e1 = float(np.sqrt(np.sum((V @ c - ref) ** 2) * dx))
            M = np.zeros((K + 1, K - 1))
            for k in range(K - 1):
                M[k, k] = 1.0; M[k + 2, k] = -k ** 2 / (k + 2) ** 2
            b, *_ = np.linalg.lstsq(V @ M, ref, rcond=None)
            e2 = float(np.sqrt(np.sum((V @ M @ b - ref) ** 2) * dx))
            rec[f"{name}.K{K}"] = e1; rec[f"{name}.K{K}.neumann"] = e2
            print(name, K, "%.2e %.2e" % (e1, e2), flush=True)
    save("floor", rec)


if __name__ == "__main__":
    {"b3": mode_b3, "inverse": mode_inverse, "timing": mode_timing, "floor": mode_floor}[sys.argv[1]]()
