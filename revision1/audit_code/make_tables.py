"""Builds markdown tables (tables_*.md fragments, printed to stdout and saved in tables.json-free form)
from the JSON outputs of the audit scripts.  Used to assemble audit_code.md without transcription errors."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda n: json.load(open(os.path.join(HERE, n)))


def e(v, nd=2):
    return f"{v:.{nd}e}"


def table(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


frag = {}

# ---------------- certification T<=1 ----------------
c = L("cert_reference.json")
rows = []
for key, rec in c.items():
    ic, eps = key.split("_")[0], key.split("_")[1][3:]
    mx = lambda d: max(d.values())
    Nst = {int(k): mx(v) for k, v in rec["N_study"].items()}
    tol = {k: mx(v["diff"]) for k, v in rec["tol_study"].items()}
    rows.append([ic, eps, e(Nst[64], 1), e(Nst[96], 1), e(Nst[128], 1), e(Nst[192], 1), e(Nst[256], 1), e(Nst[384], 1),
                 e(tol["1e-08"], 1), e(tol["1e-10"], 1), e(tol["default(1e-11)"], 1),
                 e(mx(rec["method_study"]["BDF"]["diff"]), 1), e(mx(rec["method_study"]["LSODA"]["diff"]), 1),
                 e(mx(rec["half_vs_full_N128"]), 1)])
frag["cert_T1"] = table(["IC", "eps", "N=64", "N=96", "N=128", "N=192", "N=256", "N=384",
                         "rtol 1e-8", "rtol 1e-10", "rtol 1e-11 (default)", "BDF", "LSODA", "half vs full"], rows)

# ---------------- independent check ----------------
ind = L("independent_check.json")
rows = []
for key, rec in ind.items():
    ic, eps = key.split("_")[0], key.split("_")[1][3:]
    for t, d in rec.items():
        rows.append([ic, eps, t, e(d["1000"], 1), e(d["2000"], 1), e(d["4000"], 1), e(d["8000"], 1),
                     e(d["Richardson(2000,4000)"], 1), e(d["Richardson(4000,8000)"], 1)])
frag["indep"] = table(["IC", "eps", "t", "FD Nx=1000", "2000", "4000", "8000", "Richardson(2000,4000)", "Richardson(4000,8000)"], rows)

# ---------------- long-time certification ----------------
lc = L("long_cert.json")
lc2 = L("long_cert2.json")
rows = []
for ic in ("orig", "compat"):
    for N in ("256", "384", "512", "768"):
        d = lc[ic][N]
        rows.append([ic, f"N={N} vs N=1024"] + [e(d[t], 1) for t in ("1.0", "2.0", "3.0", "5.0", "7.0", "10.0")])
    for tag, name in (("N1024_default", "N=1024 vs N=1536 (rtol 1e-12)"), ("N1280_default", "N=1280 vs N=1536"),
                      ("N1024_rtol1e-9", "N=1024, rtol 1e-9 vs N=1536"), ("N1024_BDF_1e-10", "N=1024, BDF vs N=1536")):
        d = lc2[ic][tag]
        rows.append([ic, name] + [e(d[t], 1) for t in ("1.0", "2.0", "3.0", "5.0", "7.0", "10.0")])
frag["cert_long"] = table(["IC", "comparison (eps=0.01)", "t=1", "t=2", "t=3", "t=5", "t=7", "t=10"], rows)

# ---------------- baselines vs trusted ----------------
cb = L("compare_baselines.json")
rows = []
for key, rec in cb.items():
    ic, eps = key.split("_")[0], key.split("_")[1][3:]
    for t in ("0.1", "0.5", "1.0"):
        r = rec[t]
        cells = [ic, eps, t]
        for name in ("(i)_FD_current", "(ii)_FD_ghost", "(iv)_Cheb_N80_dt5e-4"):
            m = r[name]
            cells += [e(m["L2_200"]), e(m["L2_cont"]), e(m["Linf_200"])]
        rows.append(cells)
frag["baselines"] = table(["IC", "eps", "t", "(i) L2_200", "(i) L2_cont", "(i) Linf", "(ii) L2_200", "(ii) L2_cont", "(ii) Linf",
                           "(iv) L2_200", "(iv) L2_cont", "(iv) Linf"], rows)

rows = []
for key, rec in cb.items():
    ic, eps = key.split("_")[0], key.split("_")[1][3:]
    for t in ("0.1", "0.5", "1.0"):
        r = rec[t]
        rows.append([ic, eps, t, e(r["linear_interp_floor_Nx500"]["L2_200"]), e(r["linear_interp_floor_Nx500"]["Linf_200"]),
                     e(r["(i)_FD_current"]["Linf_nodes"]), e(r["diag_current_dt2.5e-5"]["Linf_nodes"]), e(r["diag_current_Nx2000_dt1e-4"]["Linf_nodes"]),
                     e(r["(ii)_FD_ghost"]["Linf_nodes"]), e(r["diag_ghost_dt2.5e-5"]["Linf_nodes"]), e(r["diag_ghost_Nx2000_dt1e-4"]["Linf_nodes"]),
                     e(r["(i)-(ii)"]["Linf_200"])])
frag["decomp"] = table(["IC", "eps", "t", "lin-interp floor L2_200", "lin-interp floor Linf",
                        "(i) nodal Linf Nx500,dt1e-4", "(i) dt 2.5e-5", "(i) Nx 2000",
                        "(ii) nodal Linf Nx500,dt1e-4", "(ii) dt 2.5e-5", "(ii) Nx 2000", "(i) vs (ii) Linf_200"], rows)

rows = []
for key, rec in cb.items():
    ic, eps = key.split("_")[0], key.split("_")[1][3:]
    for t in ("0.1", "0.5", "1.0"):
        r = rec[t]
        rows.append([ic, eps, t, e(r["repo_style_cheb_vs_FDcurrent"]["L2"]), e(r["repo_style_cheb_vs_FDcurrent"]["Linf"]),
                     e(r["(iv)_Cheb_N80_dt5e-4"]["L2_nodes_CC"]), e(r["(iv)_Cheb_N80_dt5e-4"]["L2_cont"]), e(r["(iv)_Cheb_N80_dt5e-4"]["Linf_nodes"])])
frag["cheb_repo_vs_corrected"] = table(["IC", "eps", "t", "repo-style L2 (vs FD ref, dx_ref weights)", "repo-style Linf (vs FD ref)",
                                        "corrected L2 (Clenshaw-Curtis, vs trusted)", "corrected L2_cont (barycentric)", "corrected Linf (nodes, vs trusted)"], rows)

# ---------------- long time FD ----------------
ex = L("extra_analysis_long.json")
rows = []
for k, v in ex["long_time_FD_errors_eps0.01"].items():
    ic, t = k.split("_t")
    rows.append([ic, t, e(v["current"]["L2_200"]), e(v["current"]["Linf_200"]), e(v["ghost"]["L2_200"]), e(v["ghost"]["Linf_200"])])
frag["fd_long"] = table(["IC", "t", "(i) L2_200", "(i) Linf", "(ii) L2_200", "(ii) Linf"], rows)
rows = []
fl = ex["K_floor_long_orig"]
for K in (32, 64, 128, 256):
    rows.append([K] + [f"{fl[f'K{K}_t{t}']['L2_200']:.2e} / {fl[f'K{K}_t{t}']['max_abs_fit']:.3f}" for t in ("1.0", "2.0", "3.0", "5.0", "7.0", "10.0")])
frag["K_floor_long"] = table(["K (LS fit on 300 equispaced pts, as the code does)", "t=1", "t=2", "t=3", "t=5", "t=7", "t=10"], rows)

# ---------------- eps sweep + cheb rows ----------------
es = L("extra_analysis_epssweep.json")
rows = []
rep = L("../../data/stress_eps_sweep.json") if os.path.exists(os.path.join(HERE, "../../data/stress_eps_sweep.json")) else None
for eps, rec in es["eps_sweep_reference"].items():
    rows.append([eps, rec["nx_ref"], e(rec["FD_current_L2_200"]), e(rec["FD_current_Linf_200"])])
frag["eps_sweep_ref"] = table(["eps", "nx_ref (repo)", "FD-current ref L2_200 vs trusted", "FD-current ref Linf vs trusted"], rows)
repo_cheb = {0.01: (1.17e-4, 4.26e-4), 0.005: (3.4e-5, 2.04e-4), 0.002: (7e-6, 7.9e-5), 0.001: (3e-6, 3.9e-5)}
rows = []
tab1 = {(0.01, 40, "0.002"): (7.578e-5, 5.500e-4), (0.01, 60, "0.001"): (7.780e-5, 5.528e-4), (0.01, 80, "0.0005"): (8.511e-5, 5.563e-4)}
for key, v in es["chebyshev_baseline_corrected"].items():
    parts = key.split("_")
    eps = float(parts[0][3:])
    N = int(parts[1][1:])
    dt = parts[2][2:]
    rep_l2 = rep_linf = "-"
    if N in (40, 60, 80) and dt in ("0.002", "0.001", "0.0005"):
        rep_l2, rep_linf = tab1[(eps, N, dt)]
        rep_l2, rep_linf = e(rep_l2, 2), e(rep_linf, 2)
    elif N == 400:
        sw = [r for r in json.load(open(r"D:\Videos\KazNPU\Symmetry_PINN_2026\data\stress_eps_sweep.json"))["results"]
              if r["method"] == "chebyshev" and abs(r["eps"] - eps) < 1e-12][0]
        rep_l2, rep_linf = e(sw["l2"], 2), e(sw["linf"], 2)
    rows.append([eps, N, dt, rep_l2, rep_linf, e(v["L2_cont"]), e(v["Linf_200"]), e(v["L2_200"])])
frag["cheb_rows"] = table(["eps", "N", "dt", "stored L2 (repo)", "stored Linf (repo)", "corrected L2_cont", "corrected Linf_200", "corrected L2_200"], rows)

# ---------------- PINN recomputed ----------------
pv = L("pinn_vs_trusted.json")
rows = []
for name in ("spectral_pinn", "spectral_pinn_energy", "vanilla_pinn"):
    for t in ("0.2", "0.5", "0.8", "1.0"):
        r = pv[name][t]
        rows.append([name, t, e(r["stored_l2"], 3), e(r["stored_linf"], 3), e(r["vs_trusted_l2"], 3), e(r["vs_trusted_linf"], 3),
                     f"{100*(r['vs_trusted_l2']/r['stored_l2']-1):+.1f}%", f"{100*(r['vs_trusted_linf']/r['stored_linf']-1):+.1f}%"])
frag["pinn"] = table(["run (seed 42)", "t", "stored L2 (vs FD ref)", "stored Linf", "L2 vs trusted", "Linf vs trusted", "dL2", "dLinf"], rows)
rows = []
for K in ("16", "24", "32", "48", "64"):
    rows.append([K] + [f"{pv['K_floor_orig_eps0.01'][K][t]['L2_200']:.1e} / {pv['K_floor_orig_eps0.01'][K][t]['Linf_200']:.1e}" for t in ("0.1", "0.2", "0.5", "0.8", "1.0")])
frag["Kfloor"] = table(["K", "t=0.1", "t=0.2", "t=0.5", "t=0.8", "t=1.0"], rows)

# ---------------- behaviour ----------------
b = L("behavior.json")
rows = []
for ic in ("orig", "compat"):
    rec = b[f"{ic}_eps0.01_T10"]
    for t, v in rec.items():
        zc = ", ".join(f"{z:.4f}" for z in v["zero_crossings_x_gt0"])
        rows.append([ic, t, f"{v['max']:.4f}", f"{v['min']:.4f}", f"{v['max_abs_ux']:.1f}", f"{v['E']:.5f}", zc, f"{v['u_at'][0]:.4f}"])
frag["behav_long"] = table(["IC", "t", "max u", "min u", "max abs u_x", "E(t)", "zero x>0", "u(0,t)"], rows)
rows = []
for ic in ("orig", "compat"):
    for eps in ("0.01", "0.02", "0.05", "0.1"):
        rec = b[f"{ic}_eps{eps}_T1"]
        for t in ("0.0", "0.5", "1.0"):
            v = rec[t]
            rows.append([ic, eps, t, f"{v['max']:.4f}", f"{v['min']:.4f}", f"{v['max_abs_ux']:.2f}", f"{v['E']:.5f}", ", ".join(f"{z:.4f}" for z in v["zero_crossings_x_gt0"])])
frag["behav_eps"] = table(["IC", "eps", "t", "max u", "min u", "max abs u_x", "E", "zero x>0"], rows)
rows = []
sv = b["eps_sensitivity_T1"]
for ic in ("orig", "compat"):
    for pair in ("0.01_vs_0.02", "0.05_vs_0.1", "0.01_vs_0.1"):
        rows.append([ic, pair.replace("_vs_", " vs ")] + [f"{sv[f'{ic}_{pair}_t{t}']['Linf']:.2e} / {sv[f'{ic}_{pair}_t{t}']['L2']:.2e}" for t in ("0.1", "0.25", "0.5", "1.0")])
frag["sens"] = table(["IC", "eps pair", "t=0.1 (Linf / L2)", "t=0.25", "t=0.5", "t=1.0"], rows)
sl = b["eps_sensitivity_long"]
rows = []
for ic in ("orig", "compat"):
    rows.append([ic] + [f"{sl[f'{ic}_0.01_vs_0.02_t{t}']['Linf']:.2e} / {sl[f'{ic}_0.01_vs_0.02_t{t}']['L2']:.2e}" for t in ("1.0", "2.0", "3.0", "5.0", "10.0")])
frag["sens_long"] = table(["IC (eps 0.01 vs 0.02)", "t=1", "t=2", "t=3", "t=5", "t=10"], rows)
sw = L("extra_analysis_sens.json")["sens_where"]
rows = []
for k, v in sw.items():
    ic, e1, e2, t = k.split("_")
    rows.append([ic, f"{e1} vs {e2}", t[1:], e(v["Linf"]), f"{v['argmax_x']:.3f}", e(v["Linf_abs_x_le_0p9"]), e(v["Linf_abs_x_le_0p5"])])
frag["sens_where"] = table(["IC", "eps pair", "t", "Linf", "argmax x", "Linf on |x|<=0.9", "Linf on |x|<=0.5"], rows)
ob = L("ode_baseline.json")
rows = []
for ic in ("orig", "compat"):
    for eps in ("0.001", "0.01", "0.05", "0.1"):
        rows.append([ic, eps] + [f"{ob[f'{ic}_eps{eps}_t{t}']['L2_200']:.2e} / {ob[f'{ic}_eps{eps}_t{t}']['Linf']:.2e}" for t in ("0.1", "0.5", "1.0")])
frag["ode"] = table(["IC", "eps", "t=0.1 (L2_200 / Linf)", "t=0.5", "t=1.0"], rows)
idn = L("extra_analysis_sens.json")
rows = []
for eps in ("0.01", "0.1"):
    r = idn["identifiability_orig"][eps]
    rows.append(["orig", eps, e(r["rms_repo_data_error_vs_trusted"], 2)] + [e(r["shifts"][f]["rms"], 2) for f in ("0.95", "1.05", "1.25", "1.4", "2.0")])
for eps in ("0.01", "0.1"):
    r = idn["identifiability_compat"][eps]
    rows.append(["compat", eps, "(n/a)"] + [e(r["shifts"][f]["rms"], 2) for f in ("0.95", "1.05", "1.25", "1.4", "2.0")])
frag["ident"] = table(["IC", "eps_true", "RMS error of repo's FD-generated data vs trusted", "RMS data change for eps x0.95", "x1.05", "x1.25", "x1.4", "x2.0"], rows)

# ---------------- energy ----------------
en = L("extra_analysis_energy.json")["energy"]
rows = []
for name in ("spectral_pinn", "spectral_pinn_energy"):
    v = en[name]
    rows.append([name, f"{v['E_pinn_first']:.5f}", f"{v['E_true_first']:.5f}", f"{v['E_pinn_last']:.5f}", f"{v['E_true_last']:.5f}", e(v["max_abs_E_error"]),
                 f"{v['true_dissipation_per_step_min']:.3e} .. {v['true_dissipation_per_step_max']:.3e}", v["n_pinn_violations_tol1e10"]])
frag["energy"] = table(["run", "E_pinn(0)", "E_true(0)", "E_pinn(1)", "E_true(1)", "max abs E error over 51 pts", "true dissipation per 0.02 step", "violations counted"], rows)
rows = []
for k, v in en["trapezoid300_quadrature_error"].items():
    rows.append([k, f"{v['E_cc']:.8f}", f"{v['E_trap300']:.8f}", e(v["rel"], 1)])
frag["quad"] = table(["case", "E (Clenshaw-Curtis)", "E (trapezoid, 300 uniform pts)", "rel. diff"], rows)

# ---------------- code checks ----------------
cc = L("code_checks.json")
rows = []
for K, v in cc["bc_scale"].items():
    rows.append([K, f"{v['ux_plus']:.4g}", e(v["L_bc_repo"], 2), e(1.0 / float(v["scale_factor_K4"]), 2)])
frag["bcscale"] = table(["K", "u_x(+1) of IC fit c0", "L_bc(c0) as coded (scaled by 1/K^4)", "scale factor 1/K^4"], rows)
c2 = L("code_checks2.json")
rows = []
for K, v in c2["equispaced_fit"].items():
    rows.append([K, e(v["cond_T300"], 1), e(v["IC_fit_Linf_on_300"], 1), e(v["IC_fit_Linf_fine"], 1), f"{v['IC_fit_wall_slope']:.4g}"])
frag["eqfit"] = table(["K", "cond(Vandermonde 300 equispaced pts)", "IC-fit error at the 300 pts", "IC-fit sup error (40001 pts)", "wall slope of IC fit (true -2)"], rows)
rows = []
for k, v in cc["oracle_eps_estimator"].items():
    pass
orc = cc["oracle_eps_estimator"]
for eps in ("0.005", "0.01", "0.02", "0.05", "0.1"):
    rows.append([eps] + [f"{orc[f'eps{eps}_K{K}']['eps_hat_over_eps_true']:.3f}" for K in (16, 32, 64, 96)])
frag["oracle"] = table(["eps_true", "K=16", "K=32", "K=64", "K=96"], rows)
ws = cc["wall_slopes"]
rows = []
for t in ("0.2", "0.5", "0.8", "1.0"):
    rows.append([t, f"{ws['spectral_pinn'][t]['ux_plus1']:.3f}", f"{ws['spectral_pinn'][t]['ux_minus1']:.3f}",
                 f"{ws['trusted_K32_projection'][t]['ux_plus1']:.3f}"])
frag["wallslope"] = table(["t", "u_x(+1), Spectral-PINN seed 42", "u_x(-1)", "u_x(+1) of K=32 projection of trusted solution (true value 0)"], rows)

json.dump(frag, open(os.path.join(HERE, "tables_fragments.json"), "w"), indent=1)
for k, v in frag.items():
    print("\n#### " + k + "\n" + v)
