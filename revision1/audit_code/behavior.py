"""Part (c): qualitative/quantitative behaviour of the trusted solutions for both ICs;
eps-sensitivity; energy diagnostics; MBP; interface formation.  Writes behavior.json."""
import os, sys, json, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import refsol

XF = np.linspace(-1, 1, 20001)
out = {}


def zero_crossings(x, u):
    """positions of sign changes of u on x>0 (solution is even)."""
    m = x > 0
    xs, us = x[m], u[m]
    s = np.sign(us)
    idx = np.where(s[:-1] * s[1:] < 0)[0]
    return [float(xs[i] - us[i] * (xs[i + 1] - xs[i]) / (us[i + 1] - us[i])) for i in idx]


def energy_gl(eps, u_fun, up_fun, n=2000):
    xg, wg = np.polynomial.legendre.leggauss(n)
    return float(np.sum(wg * (0.5 * eps**2 * up_fun(xg)**2 + 0.25 * (u_fun(xg)**2 - 1)**2)))


def summarize(eps, ic, times, T):
    res, xn, nodes = refsol.reference(eps, ic, T, times, XF, return_nodes=True)
    _, D = refsol.cheb_nodes_diff(len(xn) - 1)
    rec = {}
    for t in times:
        u = res[t]
        if t == 0.0:
            ux = (refsol.u0_orig_prime if ic == "orig" else refsol.u0_compat_prime)(XF)
            E = energy_gl(eps, refsol._IC[ic], refsol.u0_orig_prime if ic == "orig" else refsol.u0_compat_prime)
            maxslope = float(np.max(np.abs(ux)))
        else:
            uxn = D @ nodes[t]
            ux = refsol.bary_interp(xn, uxn, XF)
            E = refsol.energy_cc(eps, xn, nodes[t])
            maxslope = float(np.max(np.abs(ux)))
        zc = zero_crossings(XF, u)
        rec[str(t)] = dict(
            max=float(u.max()), argmax_x=float(XF[np.argmax(u)]), min=float(u.min()),
            max_abs=float(np.abs(u).max()), L2=float(np.sqrt(np.trapezoid(u**2, XF))),
            E=E, max_abs_ux=maxslope, zero_crossings_x_gt0=zc,
            u_at=[float(refsol.bary_interp(xn, nodes[t], np.array([x0]))[0]) if t > 0 else float(refsol._IC[ic](np.array([x0]))[0]) for x0 in (0.0, 0.25, 0.5, 0.75, 1.0)],
            mbp_ok=bool(np.abs(u).max() <= 1.0 + 1e-12))
    return rec, res


sol = {}
t1 = [0.0, 0.1, 0.25, 0.5, 1.0]
for ic in ("orig", "compat"):
    for eps in (0.01, 0.02, 0.05, 0.1):
        rec, res = summarize(eps, ic, t1, 1.0)
        out[f"{ic}_eps{eps}_T1"] = rec
        sol[(ic, eps)] = res
        print(ic, eps, {t: (round(v['max'], 4), round(v['min'], 4), round(v['E'], 5), round(v['max_abs_ux'], 2)) for t, v in rec.items()}, flush=True)

# eps sensitivity (T<=1)
sens = {}
for ic in ("orig", "compat"):
    for (e1, e2) in ((0.01, 0.02), (0.05, 0.1), (0.01, 0.1)):
        for t in (0.1, 0.25, 0.5, 1.0):
            d = sol[(ic, e1)][t] - sol[(ic, e2)][t]
            sens[f"{ic}_{e1}_vs_{e2}_t{t}"] = dict(
                Linf=float(np.max(np.abs(d))), L2=float(np.sqrt(np.trapezoid(d**2, XF))),
                ref_max_abs=float(np.max(np.abs(sol[(ic, e1)][t]))))
            print(ic, e1, e2, t, f"Linf {sens[f'{ic}_{e1}_vs_{e2}_t{t}']['Linf']:.3e}  L2 {sens[f'{ic}_{e1}_vs_{e2}_t{t}']['L2']:.3e}", flush=True)
out["eps_sensitivity_T1"] = sens

# long time eps = 0.01
tl = [0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
for ic in ("orig", "compat"):
    rec, res = summarize(0.01, ic, tl, 10.0)
    out[f"{ic}_eps0.01_T10"] = rec
    sol[(ic, "long")] = res
    print(ic, "long", {t: (round(v['max'], 4), round(v['min'], 4), round(v['E'], 5), round(v['max_abs_ux'], 1), len(v['zero_crossings_x_gt0'])) for t, v in rec.items()}, flush=True)
# also eps 0.01 vs 0.02 at long times
sens_long = {}
for ic in ("orig", "compat"):
    for t in (1.0, 2.0, 3.0, 5.0, 10.0):
        res2 = refsol.reference(0.02, ic, 10.0, [1.0, 2.0, 3.0, 5.0, 10.0], XF)
        d = sol[(ic, "long")][t] - res2[t]
        sens_long[f"{ic}_0.01_vs_0.02_t{t}"] = dict(Linf=float(np.max(np.abs(d))), L2=float(np.sqrt(np.trapezoid(d**2, XF))))
        print("long sens", ic, t, sens_long[f"{ic}_0.01_vs_0.02_t{t}"], flush=True)
out["eps_sensitivity_long"] = sens_long
json.dump(out, open(os.path.join(HERE, "behavior.json"), "w"), indent=1)
