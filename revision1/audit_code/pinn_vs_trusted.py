"""(i) Recompute the saved single-seed PINN predictions (data/*.json, 200-pt grid) against
the trusted reference instead of the FD reference; (ii) K=32 Chebyshev representation floor
of the trusted solution; (iii) magnitude of reference-induced error in the stored numbers."""
import os, sys, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import refsol, compare_baselines as cb
DATA = r"D:\Videos\KazNPU\Symmetry_PINN_2026\data"
xu = np.linspace(-1, 1, 200); dxu = xu[1] - xu[0]
times = [0.2, 0.5, 0.8, 1.0]
tr = refsol.reference(0.01, "orig", 1.0, times, xu, N=256)
xf, snaps = cb.fd_solution(0.01, refsol.u0_orig, "current", T=1.0, save_times=times)
xg, snapg = cb.fd_solution(0.01, refsol.u0_orig, "ghost", T=1.0, save_times=times)
out = {}
for name in ("spectral_pinn", "spectral_pinn_energy", "vanilla_pinn"):
    d = json.load(open(os.path.join(DATA, name + ".json")))
    out[name] = {}
    for t in times:
        up = np.array(d["errors"][str(t)]["u_pred"])
        e_tr = up - tr[t]
        e_fd = up - np.interp(xu, xf, snaps[t])
        e_gh = up - np.interp(xu, xg, snapg[t])
        out[name][str(t)] = dict(
            stored_l2=d["errors"][str(t)]["l2"], stored_linf=d["errors"][str(t)]["linf"],
            recomputed_vs_FDcurrent_l2=float(np.sqrt(np.sum(e_fd**2) * dxu)),
            recomputed_vs_FDcurrent_linf=float(np.max(np.abs(e_fd))),
            vs_trusted_l2=float(np.sqrt(np.sum(e_tr**2) * dxu)),
            vs_trusted_linf=float(np.max(np.abs(e_tr))),
            vs_ghostFD_l2=float(np.sqrt(np.sum(e_gh**2) * dxu)),
            vs_ghostFD_linf=float(np.max(np.abs(e_gh))))
    print(name)
    for t in times:
        r = out[name][str(t)]
        print(f"  t={t}: stored L2 {r['stored_l2']:.4e} Linf {r['stored_linf']:.4e} | recomputed(FDcur) {r['recomputed_vs_FDcurrent_l2']:.4e} {r['recomputed_vs_FDcurrent_linf']:.4e} | vs trusted {r['vs_trusted_l2']:.4e} {r['vs_trusted_linf']:.4e}")
# K=32 floor: least squares fit on the 200-pt grid and on a 300-pt grid (as code), and truncation
fl = {}
for K in (16, 24, 32, 48, 64):
    fl[K] = {}
    for t in times + [0.1]:
        _, xn, nodes = refsol.reference(0.01, "orig", 1.0, [t], xu, N=256, return_nodes=True)
        xg300 = np.linspace(-1, 1, 300)
        uu = refsol.bary_interp(xn, nodes[t], xg300)
        Tm = np.polynomial.chebyshev.chebvander(xg300, K)
        c, *_ = np.linalg.lstsq(Tm, uu, rcond=None)
        ufit200 = np.polynomial.chebyshev.chebvander(xu, K) @ c
        e = ufit200 - refsol.bary_interp(xn, nodes[t], xu)
        fl[K][str(t)] = dict(L2_200=float(np.sqrt(np.sum(e**2) * dxu)), Linf_200=float(np.max(np.abs(e))))
    print("K=", K, {t: (f"{v['L2_200']:.2e}", f"{v['Linf_200']:.2e}") for t, v in fl[K].items()})
out["K_floor_orig_eps0.01"] = fl
json.dump(out, open(os.path.join(HERE, "pinn_vs_trusted.json"), "w"), indent=1)
