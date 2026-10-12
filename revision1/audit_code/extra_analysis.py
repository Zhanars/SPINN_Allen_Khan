"""Extra diagnostics for audit_code.md.  Usage: python extra_analysis.py [sens|long|energy|epssweep|all]
Writes extra_analysis_<part>.json."""
import os, sys, json, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import refsol, compare_baselines as cb
DATA = r"D:\Videos\KazNPU\Symmetry_PINN_2026\data"
XU = np.linspace(-1, 1, 200); DXU = XU[1] - XU[0]
XF = np.linspace(-1, 1, 20001)
out = {}
part = sys.argv[1] if len(sys.argv) > 1 else "all"


def l2_200(e):
    return float(np.sqrt(np.sum(e**2) * DXU))


_trc = {}


def _tr_cache(eps, N):
    key = (eps, N)
    if key not in _trc:
        _, xn, nodes = refsol.reference(eps, "orig", 1.0, [1.0], XU, N=N, return_nodes=True)
        _trc[key] = (xn, nodes[1.0])
    return _trc[key]


if part in ("all", "sens"):
    # 1. where does the eps-sensitivity live?
    s = {}
    for ic in ("orig", "compat"):
        for (e1, e2) in ((0.01, 0.02), (0.05, 0.1)):
            for t in (0.1, 1.0):
                r1 = refsol.reference(e1, ic, 1.0, [t], XF)[t]
                r2 = refsol.reference(e2, ic, 1.0, [t], XF)[t]
                d = np.abs(r1 - r2)
                inner = np.abs(XF) <= 0.9
                s[f"{ic}_{e1}_{e2}_t{t}"] = dict(Linf=float(d.max()), argmax_x=float(XF[d.argmax()]),
                                                 Linf_abs_x_le_0p9=float(d[inner].max()),
                                                 Linf_abs_x_le_0p5=float(d[np.abs(XF) <= 0.5].max()))
                print(ic, e1, e2, t, s[f"{ic}_{e1}_{e2}_t{t}"], flush=True)
    out["sens_where"] = s
    # 2. identifiability on the inverse-problem observation design
    xo = np.linspace(-1, 1, 40)
    tobs = [0.2, 0.4, 0.6, 0.8, 1.0]
    ident = {}
    for eps_true in (0.01, 0.1):
        base = refsol.reference(eps_true, "orig", 1.0, tobs, xo)
        # data actually used in the repo: current FD reference (np.interp to xo)
        xf, snaps = cb.fd_solution(eps_true, refsol.u0_orig, "current", T=1.0, save_times=tobs)
        data_err = np.concatenate([np.interp(xo, xf, snaps[t]) - base[t] for t in tobs])
        rec = dict(rms_repo_data_error_vs_trusted=float(np.sqrt(np.mean(data_err**2))),
                   linf_repo_data_error_vs_trusted=float(np.abs(data_err).max()), shifts={})
        for fac in (0.75, 0.95, 1.05, 1.25, 1.4, 2.0):
            tr = refsol.reference(eps_true * fac, "orig", 1.0, tobs, xo)
            d = np.concatenate([tr[t] - base[t] for t in tobs])
            rec["shifts"][str(fac)] = dict(rms=float(np.sqrt(np.mean(d**2))), linf=float(np.abs(d).max()))
        ident[str(eps_true)] = rec
        print("ident eps_true", eps_true, rec, flush=True)
    out["identifiability_orig"] = ident
    identc = {}
    for eps_true in (0.01, 0.1):
        base = refsol.reference(eps_true, "compat", 1.0, tobs, xo)
        rec = dict(shifts={})
        for fac in (0.75, 0.95, 1.05, 1.25, 1.4, 2.0):
            tr = refsol.reference(eps_true * fac, "compat", 1.0, tobs, xo)
            d = np.concatenate([tr[t] - base[t] for t in tobs])
            rec["shifts"][str(fac)] = dict(rms=float(np.sqrt(np.mean(d**2))), linf=float(np.abs(d).max()))
        identc[str(eps_true)] = rec
        print("ident(compat) eps_true", eps_true, rec, flush=True)
    out["identifiability_compat"] = identc

if part in ("all", "long"):
    # 3. long-time FD errors, eps = 0.01, repo parameters (Nx=500, dt=1e-4)
    tl = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
    lt = {}
    for ic in ("orig", "compat"):
        res, xn, nodes = refsol.reference(0.01, ic, 10.0, tl, XU, N=1024, return_nodes=True)
        tf = {t: (lambda xx, t=t: refsol.bary_interp(xn, nodes[t], xx)) for t in tl}
        t0 = time.time()
        xf, sc = cb.fd_solution(0.01, refsol._IC[ic], "current", T=10.0, save_times=tl)
        xg, sg = cb.fd_solution(0.01, refsol._IC[ic], "ghost", T=10.0, save_times=tl)
        print("FD long", ic, f"{time.time()-t0:.0f}s", flush=True)
        for t in tl:
            lt[f"{ic}_t{t}"] = dict(current=cb.fd_vs_trusted(xf, sc[t], tf[t]), ghost=cb.fd_vs_trusted(xg, sg[t], tf[t]))
            print(ic, t, {k: (round(v["L2_200"], 6), round(v["Linf_200"], 6)) for k, v in lt[f"{ic}_t{t}"].items()}, flush=True)
    out["long_time_FD_errors_eps0.01"] = lt
    # 4. K-representation floors at long times
    fl = {}
    res, xn, nodes = refsol.reference(0.01, "orig", 10.0, tl, XU, N=1024, return_nodes=True)
    x300 = np.linspace(-1, 1, 300)
    for K in (32, 64, 128, 256):
        Tm = np.polynomial.chebyshev.chebvander(x300, K)
        T200 = np.polynomial.chebyshev.chebvander(XU, K)
        for t in tl:
            u300 = refsol.bary_interp(xn, nodes[t], x300)
            c, *_ = np.linalg.lstsq(Tm, u300, rcond=None)
            uf = T200 @ c
            e = uf - res[t]
            fl[f"K{K}_t{t}"] = dict(L2_200=l2_200(e), Linf_200=float(np.abs(e).max()),
                                    max_abs_fit=float(np.abs(uf).max()))
        print("K", K, {t: (f"{fl[f'K{K}_t{t}']['L2_200']:.3e}", f"{fl[f'K{K}_t{t}']['max_abs_fit']:.3f}") for t in tl}, flush=True)
    out["K_floor_long_orig"] = fl
    d = json.load(open(os.path.join(DATA, "stress_time_sweep.json")))
    out["pinn_T10_reported_l2"] = {k: v["l2"] for k, v in d["results"][2]["errors"].items()}

if part in ("all", "energy"):
    # 5. energy: PINN E_traj vs trusted energy; quadrature error of trapezoid-300
    en = {}
    xg, wg = np.polynomial.legendre.leggauss(2000)
    E0 = float(np.sum(wg * (0.5 * 0.01**2 * refsol.u0_orig_prime(xg)**2 + 0.25 * (refsol.u0_orig(xg)**2 - 1)**2)))
    for name in ("spectral_pinn", "spectral_pinn_energy"):
        d = json.load(open(os.path.join(DATA, name + ".json")))
        tg = np.array(d["t_grid"])
        E = np.array(d["E_traj"])
        res, xn, nodes = refsol.reference(0.01, "orig", 1.0, [t for t in tg if t > 0], XU, N=256, return_nodes=True)
        Etr = np.array([refsol.energy_cc(0.01, xn, nodes[t]) if t > 0 else E0 for t in tg])
        dE_true = np.diff(Etr)
        dE_pinn = np.diff(E)
        en[name] = dict(max_abs_E_error=float(np.abs(E - Etr).max()), E_pinn_first=float(E[0]), E_pinn_last=float(E[-1]),
                        E_true_first=float(Etr[0]), E_true_last=float(Etr[-1]),
                        true_dissipation_per_step_min=float(-dE_true.max()), true_dissipation_per_step_max=float(-dE_true.min()),
                        pinn_dE_max=float(dE_pinn.max()), pinn_dE_min=float(dE_pinn.min()),
                        n_pinn_violations_tol1e10=int((dE_pinn > 1e-10).sum()), n_steps=int(len(dE_pinn)),
                        monotone_true=bool((dE_true < 0).all()))
        print(name, en[name], flush=True)
    qe = {}
    for ic in ("orig", "compat"):
        for eps in (0.01, 0.1):
            ts = [0.5, 1.0] if eps == 0.1 else [0.5, 1.0, 5.0, 10.0]
            res, xn, nodes = refsol.reference(eps, ic, max(ts), ts, XU, return_nodes=True)
            _, D = refsol.cheb_nodes_diff(len(xn) - 1)
            x300 = np.linspace(-1, 1, 300)
            for t in ts:
                u = refsol.bary_interp(xn, nodes[t], x300)
                ux = refsol.bary_interp(xn, D @ nodes[t], x300)
                Etrap = float(np.trapezoid(0.5 * eps**2 * ux**2 + 0.25 * (u**2 - 1)**2, x300))
                Ecc = refsol.energy_cc(eps, xn, nodes[t])
                qe[f"{ic}_eps{eps}_t{t}"] = dict(E_cc=Ecc, E_trap300=Etrap, abs_diff=abs(Etrap - Ecc), rel=abs(Etrap - Ecc) / Ecc)
                print("quad", ic, eps, t, qe[f"{ic}_eps{eps}_t{t}"], flush=True)
    en["trapezoid300_quadrature_error"] = qe
    out["energy"] = en
    mb = {}
    for name in ("spectral_pinn", "spectral_pinn_energy", "vanilla_pinn"):
        d = json.load(open(os.path.join(DATA, name + ".json")))
        mb[name] = {t: dict(max=float(np.max(v["u_pred"])), min=float(np.min(v["u_pred"]))) for t, v in d["errors"].items()}
    out["pinn_stored_range"] = mb
    print(mb)

if part in ("all", "epssweep"):
    # 6. reference quality in the eps sweep (stress_eps_sweep.py): FD ref with Nx=reference_nx_for_eps, dt=1e-4
    sys.path.insert(0, r"D:\Videos\KazNPU\Symmetry_PINN_2026\code")
    from stress_common import reference_nx_for_eps, chebyshev_n_for_eps
    es = {}
    for eps in (0.01, 0.005, 0.002, 0.001):
        conv = {}
        Ns = [512, 1024, 1536] if eps < 0.01 else [256, 384]
        sols = {}
        for N in Ns:
            info = {}
            t0 = time.time()
            sols[N] = refsol.reference(eps, "orig", 1.0, [1.0], XU, N=N, info=info)[1.0]
            conv[N] = dict(time_s=time.time() - t0, nsteps=info["nsteps"])
        for N in Ns[:-1]:
            conv[N]["Linf_vs_largest"] = float(np.abs(sols[N] - sols[Ns[-1]]).max())
        tr = sols[Ns[-1]]
        nx = reference_nx_for_eps(eps)
        xf, snaps = cb.fd_solution(eps, refsol.u0_orig, "current", Nx=nx, T=1.0, save_times=[1.0])
        u200 = np.interp(XU, xf, snaps[1.0])
        e = u200 - tr
        es[str(eps)] = dict(conv=conv, nx_ref=nx, FD_current_L2_200=l2_200(e), FD_current_Linf_200=float(np.abs(e).max()))
        print("eps sweep", eps, es[str(eps)], flush=True)
    out["eps_sweep_reference"] = es
    # corrected numbers for the Chebyshev baseline rows (Table 1 and the eps sweep)
    from scipy.linalg import lu_factor, lu_solve
    xfine = np.linspace(-1, 1, 20001)

    def cheb_run_lu(eps, N, dt, T=1.0):
        D, x = cb.repo_cheb_diff_matrix(N)
        D2 = D @ D
        A = np.eye(N + 1) - dt * eps**2 * D2
        A[0, :] = D[0, :]
        A[-1, :] = D[-1, :]
        lu = lu_factor(A)
        u = refsol.u0_orig(x).copy()
        for n in range(int(round(T / dt))):
            rhs = u - dt * (u**3 - u)
            rhs[0] = 0.0
            rhs[-1] = 0.0
            u = lu_solve(lu, rhs)
        return x, u

    cb_rows = {}
    for (eps, N, dt) in [(0.01, 40, 2e-3), (0.01, 60, 1e-3), (0.01, 80, 5e-4), (0.01, 80, 5e-5),
                         (0.01, 400, 1e-4), (0.005, 400, 5e-5), (0.002, 400, 2e-5), (0.001, 400, 1e-5)]:
        t0 = time.time()
        x, u = cheb_run_lu(eps, N, dt)
        el = time.time() - t0
        Nt = 512 if eps >= 0.01 else 1536
        tr_nodes = refsol.reference(eps, "orig", 1.0, [1.0], xfine, N=Nt)[1.0]
        m = cb.cheb_vs_trusted(x, u, lambda xx: refsol.bary_interp(*_tr_cache(eps, Nt), xx))
        cb_rows[f"eps{eps}_N{N}_dt{dt:g}"] = dict(m, elapsed_s=el)
        print("cheb row", eps, N, dt, {k: f"{v:.3e}" for k, v in m.items()}, flush=True)
    out["chebyshev_baseline_corrected"] = cb_rows
json.dump(out, open(os.path.join(HERE, f"extra_analysis_{part}.json"), "w"), indent=1)
