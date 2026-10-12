"""Pilot: does the Neumann-incompatible IC matter for the Spectral-PINN? (K=32/64, orig vs compat IC).
usage: python pilot_sp.py <ic: orig|compat> <K> <seed> <epochs> <out.json>"""
import os, sys, json, time
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))
sys.path.insert(0, HERE)
import numpy as np
import torch
torch.set_num_threads(1)
import pilot_ref as pr
import spectral_pinn as sp

ic, K, seed, epochs, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
u0 = pr.u0_orig if ic == "orig" else pr.u0_compat
sp.u0_fn = u0                       # the model builds its IC coefficients through this name
eps = 0.01
times = [0.2, 0.5, 0.8, 1.0]
xn, ref = pr.reference(eps, u0, times)
xe = np.linspace(-1, 1, 200)
dx = xe[1] - xe[0]

model, hist, elapsed = sp.train_spectral_pinn(eps, K=K, epochs=epochs, energy_constrained=False, seed=seed, verbose=False)
res = {"ic": ic, "K": K, "seed": seed, "epochs": epochs, "elapsed": elapsed, "errors": {}}
for t in times:
    up = sp.eval_on_grid(model, t, xe)
    ur = pr.bary_eval(xn, ref[t], xe)
    res["errors"][str(t)] = {"l2": float(np.sqrt(np.sum((up - ur) ** 2) * dx)), "linf": float(np.max(np.abs(up - ur))),
                             "max_u": float(np.max(np.abs(up)))}
# final training losses
res["final"] = {k: hist[k][-1] for k in ("loss", "l_pde", "l_ic", "l_bc")}
# wall slope of the network at t=0 and t=1 (should be 0 for t>0)
with torch.no_grad():
    for t in (0.0, 1e-3, 0.5, 1.0):
        c = model.coeffs(torch.tensor([[t]]))
        res.setdefault("wall_slope", {})[str(t)] = [float(c @ model.v_plus * (K ** 2)), float(c @ model.v_minus * (K ** 2))]
json.dump(res, open(out, "w"), indent=1)
print(json.dumps({k: res[k] for k in ("ic", "K", "elapsed")}), res["errors"]["1.0"], res["final"])
