"""Classical (non-neural) inverse baseline and identifiability diagnostics for eps.
Forward model: Chebyshev collocation (pilot_ref.reference).  Data: independent solver (N=256 Radau) + Gaussian noise.
Estimator: nonlinear least squares in log(eps) from eps_init = 2*eps_true (the same start as the PINN)."""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import least_squares
import pilot_ref as pr

obs_times = [0.2, 0.4, 0.6, 0.8, 1.0]
x_obs = np.linspace(-1, 1, 40)
out = {}
rng = np.random.default_rng(7)


def model(eps, u0, N=64):
    xn, ref = pr.reference(eps, u0, obs_times, N=N, rtol=1e-9, atol=1e-11)
    return np.concatenate([pr.bary_eval(xn, ref[t], x_obs) for t in obs_times])


for icname, u0 in (("compat", pr.u0_compat), ("orig", pr.u0_orig)):
    for eps_true in (0.01, 0.05, 0.1):
        xn, ref = pr.reference(eps_true, u0, obs_times, N=256)
        clean = np.concatenate([pr.bary_eval(xn, ref[t], x_obs) for t in obs_times])
        # sensitivity: d data / d log(eps)  (central difference) and Fisher standard error for sigma
        h = 1e-3
        dlog = (model(eps_true * np.exp(h), u0) - model(eps_true * np.exp(-h), u0)) / (2 * h)
        S = float(np.linalg.norm(dlog))
        rec = {"sensitivity_norm_dlogeps": S, "sensitivity_rms": float(np.sqrt(np.mean(dlog ** 2))), "runs": {}}
        for sigma in (0.0, 0.01, 0.05):
            ests = []
            for rep in range(1 if sigma == 0 else 20):
                data = clean + (rng.normal(0, sigma, clean.shape) if sigma > 0 else 0)
                f = lambda p: model(np.exp(p[0]), u0) - data
                r = least_squares(f, [np.log(2 * eps_true)], method="lm" if False else "trf", xtol=1e-12, ftol=1e-12, gtol=1e-12, x_scale=1.0)
                ests.append(float(np.exp(r.x[0])))
            ests = np.array(ests)
            rec["runs"][str(sigma)] = {"eps_hat_mean": float(ests.mean()), "eps_hat_sd": float(ests.std(ddof=1)) if len(ests) > 1 else 0.0,
                                       "rel_bias_pct": float(100 * (ests.mean() - eps_true) / eps_true),
                                       "rel_rmse_pct": float(100 * np.sqrt(np.mean((ests - eps_true) ** 2)) / eps_true),
                                       "fisher_rel_se_pct": float(100 * sigma / S) if sigma > 0 else 0.0}
        out[f"{icname}_{eps_true}"] = rec
        print(icname, eps_true, "S=%.3e" % S, {k: (round(v["rel_bias_pct"], 3), round(v["rel_rmse_pct"], 3)) for k, v in rec["runs"].items()}, flush=True)
json.dump(out, open("classical_inverse.json", "w"), indent=1)
