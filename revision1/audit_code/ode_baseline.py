"""How much does the diffusion term matter in the benchmark?  Error of the eps=0 pointwise ODE solution
u' = u - u^3 (no diffusion at all) against the trusted PDE solution on the 200-point PINN grid.
Writes ode_baseline.json."""
import os, sys, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import refsol
xu = np.linspace(-1, 1, 200); dx = xu[1] - xu[0]
def ode(u0, t):
    # exact solution of u' = u - u^3:  u(t) = u0 e^t / sqrt(1 + u0^2 (e^{2t} - 1))
    return u0 * np.exp(t) / np.sqrt(1.0 + u0**2 * (np.exp(2 * t) - 1.0))
out = {}
for ic in ("orig", "compat"):
    for eps in (0.001, 0.01, 0.05, 0.1):
        res = refsol.reference(eps, ic, 1.0, [0.1, 0.5, 1.0], xu)
        for t in (0.1, 0.5, 1.0):
            e = ode(refsol._IC[ic](xu), t) - res[t]
            out[f"{ic}_eps{eps}_t{t}"] = dict(L2_200=float(np.sqrt(np.sum(e**2) * dx)), Linf=float(np.abs(e).max()), argmax_x=float(xu[np.abs(e).argmax()]))
        print(ic, eps, {t: (f"{out[f'{ic}_eps{eps}_t{t}']['L2_200']:.2e}", f"{out[f'{ic}_eps{eps}_t{t}']['Linf']:.2e}", round(out[f'{ic}_eps{eps}_t{t}']['argmax_x'],2)) for t in (0.1, 0.5, 1.0)})
json.dump(out, open(os.path.join(HERE, "ode_baseline.json"), "w"), indent=1)

# ODE vs stored PINN errors at the PINN evaluation times (needs pinn_vs_trusted.json from pinn_vs_trusted.py)
pj = os.path.join(HERE, "pinn_vs_trusted.json")
if os.path.exists(pj):
    pv = json.load(open(pj))
    ts = [0.2, 0.5, 0.8, 1.0]
    tr = refsol.reference(0.01, "orig", 1.0, ts, xu)
    o2 = {}
    for t in ts:
        e = ode(refsol.u0_orig(xu), t) - tr[t]
        o2[str(t)] = dict(ode_L2=float(np.sqrt(np.sum(e**2) * dx)), ode_Linf=float(np.abs(e).max()),
                          pinn_L2=pv["spectral_pinn"][str(t)]["vs_trusted_l2"], pinn_Linf=pv["spectral_pinn"][str(t)]["vs_trusted_linf"],
                          vanilla_L2=pv["vanilla_pinn"][str(t)]["vs_trusted_l2"])
    json.dump(o2, open(os.path.join(HERE, "ode_vs_pinn.json"), "w"), indent=1)
    print("ode_vs_pinn.json written")
