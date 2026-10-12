"""Recompute every JSON-derivable number quoted in manuscript.tex.
Run:  python -I compute_all.py > compute_all_output.txt
Reads only ../../data/*.json (read-only)."""
import json, os, sys
import numpy as np

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data")
J = lambda n: json.load(open(os.path.join(D, n)))


def ms(a, ddof=1):
    a = np.asarray(a, float)
    return a.mean(), a.std(ddof=ddof)


def hdr(s):
    print("\n" + "=" * 8, s, "=" * 8)


# ---------------------------------------------------------------- Table 1
hdr("TABLE 1 (forward_multiseed.json, chebyshev_baseline.json)")
fm = J("forward_multiseed.json")["results"]
by = {}
for r in fm:
    by.setdefault(r["method"], []).append(r)
for m, rs in by.items():
    l2 = [r["l2"] for r in rs]
    li = [r["linf"] for r in rs]
    el = [r["elapsed"] for r in rs]
    print(f"{m:16s} n={len(rs)} seeds={[r['seed'] for r in rs]}")
    for ddof in (1, 0):
        a, b = ms(l2, ddof); c, d = ms(li, ddof); e, f = ms(el, ddof)
        print(f"   ddof={ddof}: L2 {a:.4e} +- {b:.3e} | Linf {c:.4e} +- {d:.3e} | time {e:.2f} +- {f:.2f}")
    print("   per-seed L2:", ["%.4e" % x for x in l2])
    print("   SD/mean (L2, ddof=1): %.3f" % (np.std(l2, ddof=1) / np.mean(l2)))

cb = J("chebyshev_baseline.json")["results"]
print("Chebyshev baseline rows:")
for r in cb:
    print("  N=%d dt=%g L2=%.4e Linf=%.4e t=%.3fs" % (r["N"], r["dt"], r["l2"], r["linf"], r["elapsed"]))

# ratios 49 vs 47
hdr("RATIO 49 vs 47 (vanilla / spectral L2)")
van = {r["seed"]: r for r in by["vanilla"]}
spe = {r["seed"]: r for r in by["spectral"]}
spe_e = {r["seed"]: r for r in by["spectral_energy"]}
vm = np.mean([r["l2"] for r in van.values()])
sm = np.mean([r["l2"] for r in spe.values()])
sem = np.mean([r["l2"] for r in spe_e.values()])
print("ratio of 5-seed means  vanilla/spectral        = %.3f" % (vm / sm))
print("ratio of 5-seed means  vanilla/spectral_energy = %.3f" % (vm / sem))
print("ratio of rounded table values 6.74e-2/1.38e-3  = %.3f" % (6.74e-2 / 1.38e-3))
pr = [van[s]["l2"] / spe[s]["l2"] for s in van]
print("paired per-seed ratios (same seed index)       =", ["%.2f" % x for x in pr], " mean %.2f  median %.2f  geo-mean %.2f" %
      (np.mean(pr), np.median(pr), np.exp(np.mean(np.log(pr)))))
pre = [van[s]["l2"] / spe_e[s]["l2"] for s in van]
print("paired ratios vanilla/spectral_energy          =", ["%.2f" % x for x in pre], " mean %.2f" % np.mean(pre))
print("seed 42 only (multiseed json)                  = %.3f (energy variant %.3f)" % (van[42]["l2"] / spe[42]["l2"], van[42]["l2"] / spe_e[42]["l2"]))
v1 = J("vanilla_pinn.json"); s1 = J("spectral_pinn.json"); e1 = J("spectral_pinn_energy.json")
print("single-run files (seed 42): vanilla_pinn.json t=1 L2 %.4e / spectral_pinn.json %.4e = %.3f ; / spectral_pinn_energy.json %.4e = %.3f" % (
    v1["errors"]["1.0"]["l2"], s1["errors"]["1.0"]["l2"], v1["errors"]["1.0"]["l2"] / s1["errors"]["1.0"]["l2"],
    e1["errors"]["1.0"]["l2"], v1["errors"]["1.0"]["l2"] / e1["errors"]["1.0"]["l2"]))
# bootstrap CI for ratio of means
rng = np.random.default_rng(0)
vl = np.array([van[s]["l2"] for s in van]); sl = np.array([spe[s]["l2"] for s in van])
bs = []
for _ in range(20000):
    i = rng.integers(0, 5, 5); j = rng.integers(0, 5, 5)
    bs.append(vl[i].mean() / sl[j].mean())
print("bootstrap (independent resampling) ratio-of-means 95%% CI: [%.1f, %.1f]" % tuple(np.percentile(bs, [2.5, 97.5])))
# delta method CI with SD of means
a, sa = vl.mean(), vl.std(ddof=1) / np.sqrt(5); b, sb = sl.mean(), sl.std(ddof=1) / np.sqrt(5)
r0 = a / b; se = r0 * np.sqrt((sa / a) ** 2 + (sb / b) ** 2)
print("delta-method ratio %.2f +- %.2f (1 SE)" % (r0, se))
# Linf ratio
print("Linf ratio of means vanilla/spectral = %.2f" % (np.mean([r['linf'] for r in by['vanilla']]) / np.mean([r['linf'] for r in by['spectral']])))

hdr("TIME RATIO 5.4 / 40x degradation")
tv = np.mean([r["elapsed"] for r in by["vanilla"]]); ts = np.mean([r["elapsed"] for r in by["spectral"]])
te = np.mean([r["elapsed"] for r in by["spectral_energy"]]); t64 = np.mean([r["elapsed"] for r in by["spectral_K64"]])
print("mean times  vanilla %.1f  spectral %.1f  spectral_energy %.1f  K64 %.1f" % (tv, ts, te, t64))
print("ratio vanilla/spectral = %.3f ; vanilla/spectral_energy = %.3f ; table-rounded 919/169 = %.3f ; 919/191=%.3f" % (tv / ts, tv / te, 919 / 169, 919 / 191))
print("per-seed ratio of elapsed:", ["%.2f" % (van[s]["elapsed"] / spe[s]["elapsed"]) for s in van])
print("single-run files: vanilla %.1fs spectral %.1fs energy %.1fs -> ratio %.2f (energy %.2f)" % (v1["elapsed"], s1["elapsed"], e1["elapsed"], v1["elapsed"] / s1["elapsed"], v1["elapsed"] / e1["elapsed"]))
print("vanilla epochs 8000 ; spectral 6000 (+1000 warm-up NOT included in elapsed, t0 set after warm-up)")
k64 = np.array([r["l2"] for r in by["spectral_K64"]])
print("K64/K32 ratio of means = %.2f ; paired per seed =" % (k64.mean() / sl.mean()), ["%.1f" % (k64[i] / sl[i]) for i in range(5)])
sl_ord = [spe[s]["l2"] for s in van]
print("K=64 per-seed L2:", ["%.4e" % x for x in k64], " min K64 %.3e vs max K32 %.3e ; SD/mean = %.3f" % (k64.min(), max(sl_ord), k64.std(ddof=1) / k64.mean()))

# ---------------------------------------------------------------- energy
hdr("ENERGY VIOLATIONS (single-run files, seed 42)")
for nm, d in (("spectral_pinn.json (NO penalty)", s1), ("spectral_pinn_energy.json (penalty)", e1)):
    E = np.array(d["E_traj"])
    print(nm, "n_energy_violations=%d, n_points=%d, max consecutive diff %.3e, E0=%.6f E1=%.6f" % (d["n_energy_violations"], len(E), np.diff(E).max(), E[0], E[-1]))
    print("   recomputed count of increases > 1e-10: %d ; count of any increase > 0: %d ; min decrement %.3e" % (int((np.diff(E) > 1e-10).sum()), int((np.diff(E) > 0).sum()), -np.diff(E).max()))
Ea = np.array(s1["E_traj"]); Eb = np.array(e1["E_traj"])
print("max |E_pen - E_nopen| over the 51 points = %.3e  (relative %.3e)" % (np.abs(Ea - Eb).max(), (np.abs(Ea - Eb) / Ea).max()))
eh = np.array(e1["energy_history"])
print("energy_history (L_E sampled every 50 epochs, 120 values): nonzero %d, max %.3e, sum %.3e" % ((eh > 0).sum(), eh.max(), eh.sum()))
print("   -> penalty term lam_E*L_E max = %.3e (lam_E=50 in code)" % (50 * eh.max()))
print("   loss_history final: no-penalty %.4e, penalty %.4e" % (s1["loss_history"][-1], e1["loss_history"][-1]))

# ---------------------------------------------------------------- inverse
hdr("TABLE 2 INVERSE SWEEP (inverse_bias_sweep.json)")
iv = J("inverse_bias_sweep.json")["results"]
eps_l = sorted(set(r["eps_true"] for r in iv)); K_l = sorted(set(r["K"] for r in iv))
print("n=%d, seeds=%s, eps_init/eps_true=%s" % (len(iv), set(r["seed"] for r in iv), set(round(r["eps_init"] / r["eps_true"], 6) for r in iv)))
allv = np.array([r["rel_err_pct"] for r in iv])
for ddof in (1, 0):
    print("overall (n=20) mean %.3f  SD(ddof=%d) %.3f" % (allv.mean(), ddof, allv.std(ddof=ddof)))
print("signed relative error (eps_hat-eps_true)/eps_true*100: min %.2f max %.2f (all overestimates: %s)" % (
    min((r['eps_hat'] - r['eps_true']) / r['eps_true'] * 100 for r in iv), max((r['eps_hat'] - r['eps_true']) / r['eps_true'] * 100 for r in iv),
    all(r['eps_hat'] > r['eps_true'] for r in iv)))
print("table:")
print("eps\\K " + " ".join("%7d" % k for k in K_l) + "   row mean+-SD")
for e in eps_l:
    row = [next(r["rel_err_pct"] for r in iv if r["eps_true"] == e and r["K"] == k) for k in K_l]
    print("%6g " % e + " ".join("%7.2f" % x for x in row) + "   %.1f +- %.1f" % (np.mean(row), np.std(row, ddof=1)))
for k in K_l:
    col = [r["rel_err_pct"] for r in iv if r["K"] == k]
    print("K=%3d col mean %.2f SD(ddof=1) %.2f SD(ddof=0) %.2f" % (k, np.mean(col), np.std(col, ddof=1), np.std(col, ddof=0)))
mask = np.array([not (r["eps_true"] == 0.01 and r["K"] == 16) for r in iv])
print("excluding the 12.6%% outlier: mean %.2f SD %.2f (n=19)" % (allv[mask].mean(), allv[mask].std(ddof=1)))
print("K=16 col without outlier: mean %.2f" % np.mean([r['rel_err_pct'] for r in iv if r['K'] == 16 and not (r['eps_true'] == 0.01)]))
print("Trend vs K at fixed eps (rel err list K=16,32,64,96):")
for e in eps_l:
    row = [next(r["rel_err_pct"] for r in iv if r["eps_true"] == e and r["K"] == k) for k in K_l]
    d = np.diff(row)
    print("  eps=%g:" % e, ["%.1f" % x for x in row], "decreases (K->next K):", int((d < 0).sum()), "of 3 ; net change K16->K96 %+.1f" % (row[-1] - row[0]),
          "; Spearman-like sign(corr K,err) %+.2f" % np.corrcoef(np.log(K_l), row)[0, 1])
print("Trend vs eps at fixed K:")
for k in K_l:
    col = [next(r["rel_err_pct"] for r in iv if r["eps_true"] == e and r["K"] == k) for e in eps_l]
    print("  K=%d:" % k, ["%.1f" % x for x in col], "corr(log eps, err) %+.2f" % np.corrcoef(np.log(eps_l), col)[0, 1])
print("Compare noise-experiment numbers quoted in text (eps_hat=0.01402,0.01402,0.01400; 40.2/40.2/40.0%) with the sweep at eps=0.01:")
for r in iv:
    if r["eps_true"] == 0.01:
        print("   K=%d eps_hat=%.6f rel_err=%.3f%%" % (r["K"], r["eps_hat"], r["rel_err_pct"]))
print("   sweep eps_true=0.05,K=32 -> %.2f%% (text: single run at eps_true=0.05 gave 14.7%%)" % next(r["rel_err_pct"] for r in iv if r["eps_true"] == 0.05 and r["K"] == 32))
print("   eps_hat/eps_true ratios:", sorted(set(round(r["eps_hat"] / r["eps_true"], 3) for r in iv)))
print("   (all runs start at eps_init = 2*eps_true; 19 of 20 final estimates lie at 1.30-1.45 x truth, i.e. between truth and init; 1 run (eps=0.01,K=16) ends at 0.874 x truth = UNDER-estimate)")
sg=[(r["eps_hat"]-r["eps_true"])/r["eps_true"]*100 for r in iv]
print("   signed mean rel err = %.2f%%, signed SD(ddof=1) = %.2f; #overestimates=%d, #underestimates=%d" % (np.mean(sg), np.std(sg,ddof=1), sum(x>0 for x in sg), sum(x<0 for x in sg)))

# ---------------------------------------------------------------- stress K
hdr("TABLE 3 STRESS K (stress_k_sweep.json)")
sk = J("stress_k_sweep.json")["results"]
for r in sk:
    print("K=%3d %-16s warm %5d + ep %5d  L2 %.4e  Linf %.4e  t=%.0fs (%.1f min)" % (r["K"], r["tag"], r["warmup"], r["epochs"], r["l2"], r["linf"], r["elapsed"], r["elapsed"] / 60))
k256a = next(r for r in sk if r["tag"] == "K256_1x"); k256b = next(r for r in sk if r["tag"] == "K256_8x_matched")
print("K=256: extra training time of 8x vs 1x = %.0f s = %.1f min ; total 8x = %.1f min ; 1x = %.1f min ; L2 ratio 1x/8x = %.3f ; Linf %.2f vs %.2f" % (
    k256b["elapsed"] - k256a["elapsed"], (k256b["elapsed"] - k256a["elapsed"]) / 60, k256b["elapsed"] / 60, k256a["elapsed"] / 60, k256a["l2"] / k256b["l2"], k256a["linf"], k256b["linf"]))
print("matched-budget improvement factors: K64 %.2fx, K128 %.2fx, K256 %.2fx" % (
    next(r for r in sk if r['tag'] == 'K64_1x')['l2'] / next(r for r in sk if r['tag'] == 'K64_2x_matched')['l2'],
    next(r for r in sk if r['tag'] == 'K128_1x')['l2'] / next(r for r in sk if r['tag'] == 'K128_4x_matched')['l2'],
    k256a['l2'] / k256b['l2']))
print("K64_1x equals forward_multiseed seed-42 K64 L2? %.6e vs %.6e" % (next(r for r in sk if r['tag'] == 'K64_1x')['l2'], next(r for r in by['spectral_K64'] if r['seed'] == 42)['l2']))
print("matched L2 relative to K=32 mean 1.378e-3: K64 %.1fx  K128 %.1fx  K256 %.0fx" % tuple(next(r for r in sk if r['tag'] == t)['l2'] / sl.mean() for t in ('K64_2x_matched', 'K128_4x_matched', 'K256_8x_matched')))
print("baseline_K32 stored in json:", J("stress_k_sweep.json")["baseline_K32"], " (manually typed 1.38e-3/3.32e-3 = Table 1 rounded means)")

# ---------------------------------------------------------------- stress eps
hdr("TABLE 4 STRESS EPS (stress_eps_sweep.json)")
se = J("stress_eps_sweep.json")["results"]
for r in se:
    print("%-13s eps=%-6g nx_ref=%5d N=%s dt=%s  L2 %.4e Linf %.4e t=%.0fs" % (r["method"], r["eps"], r["nx_ref"], r.get("N", "-"), r.get("dt", "-"), r["l2"], r["linf"], r["elapsed"]))
for e in (0.01, 0.005, 0.002, 0.001):
    r = next(r for r in se if r["method"] == "chebyshev" and r["eps"] == e)
    print("  eps=%g: Chebyshev code-norm weight sum = (N+1)*dx_ref = %d*%.3e = %.3f (true measure 2) -> norm factor sqrt(%.3f/2)=%.3f" % (
        e, r["N"] + 1, 2.0 / r["nx_ref"], (r["N"] + 1) * 2.0 / r["nx_ref"], (r["N"] + 1) * 2.0 / r["nx_ref"], np.sqrt((r["N"] + 1) * 2.0 / r["nx_ref"] / 2)))

# ---------------------------------------------------------------- stress T
hdr("TABLE 5 STRESS T (stress_time_sweep.json)")
st = J("stress_time_sweep.json")["results"]
for r in st:
    errs = r["errors"]
    tmax = max(float(k) for k in errs)
    print("T=%g epochs=%d nt_col=%d n_ckpt(penalty)=%d  eval points=%d  violations=%d  L2(t=1)=%.4e  L2(t=T)=%.4e  Linf(t=T)=%.3f  E_first=%.4f E_last=%.4f elapsed=%.0fs" % (
        r["T"], r["epochs"], r["nt_col"], r["n_ckpt"], r["n_eval_points"], r["n_energy_violations"], errs["1.0"]["l2"], errs[str(tmax)]["l2"], errs[str(tmax)]["linf"], r["E_first"], r["E_last"], r["elapsed"]))
    print("   eval spacing = T/(n_eval-1) = %.3f ; penalty checkpoint spacing = T/n_ckpt = %.3f ; evaluation is %dx denser than the penalty checkpoints" % (r["T"] / (r["n_eval_points"] - 1), r["T"] / r["n_ckpt"], round(r["n_ckpt"] * 2 / r["T"] / (r["n_ckpt"] / r["T"]))))
    print("   L2 vs t:", {k: round(v["l2"], 3) for k, v in errs.items()})
print("T=10: first t with L2>=0.95:", next(k for k, v in st[2]["errors"].items() if v["l2"] >= 0.95))
print("T=10: L2 for t>=5: min %.3f max %.3f" % (min(v['l2'] for k, v in st[2]['errors'].items() if float(k) >= 5), max(v['l2'] for k, v in st[2]['errors'].items() if float(k) >= 5)))
print("T=1 stress run equals multiseed spectral_energy seed 42? %.6e vs %.6e" % (st[0]["errors"]["1.0"]["l2"], spe_e[42]["l2"]))

# ---------------------------------------------------------------- 2D
hdr("TABLE 6 (experiment_2d.json + symmetry_2d_multiseed.json)")
e2 = J("experiment_2d.json")
for t, v in e2["violations"].items():
    print("t=%s Fourier S=%.4e  vanilla(seed42 single run) S=%.4e" % (t, v["fourier_spectral"], v["vanilla"]))
print("times: fourier %.1fs vanilla %.1fs (not quoted in the manuscript)" % (e2["fourier_time"], e2["vanilla_time"]))
m2 = J("symmetry_2d_multiseed.json")["results"]
for t in ("0.0", "0.2", "0.5", "1.0"):
    a = np.array([r["violations"][t] for r in m2])
    print("t=%s vanilla 5 seeds: mean %.4e SD(ddof=1) %.4e (ddof=0 %.4e)  -> (%.2f +- %.2f)e%d   per-seed %s" % (
        t, a.mean(), a.std(ddof=1), a.std(ddof=0), a.mean() / 10 ** np.floor(np.log10(a.mean())), a.std(ddof=1) / 10 ** np.floor(np.log10(a.mean())), np.floor(np.log10(a.mean())), ["%.2e" % x for x in a]))
for r in m2:
    v = r["violations"]
    ts = ["0.0", "0.2", "0.5", "1.0"]
    mono = all(v[ts[i + 1]] > v[ts[i]] for i in range(3))
    print("  seed %4d monotone increasing: %s ; S(1)/S(0) = %.2f" % (r["seed"], mono, v["1.0"] / v["0.0"]))
print("growth of mean S: S(1)/S(0) = %.2f  (log10 = %.2f decades)" % (np.mean([r['violations']['1.0'] for r in m2]) / np.mean([r['violations']['0.0'] for r in m2]),
      np.log10(np.mean([r['violations']['1.0'] for r in m2]) / np.mean([r['violations']['0.0'] for r in m2]))))
print("Fourier S range in Table 6: %.2e .. %.2e" % (min(v['fourier_spectral'] for v in e2['violations'].values()), max(v['fourier_spectral'] for v in e2['violations'].values())))
print("vanilla seed-42 in multiseed json vs experiment_2d.json: %.6e vs %.6e" % (m2[0]["violations"]["0.0"], e2["violations"]["0.0"]["vanilla"]))

hdr("TABLE 7 (stress_2d_sweep.json) and the sentence 'S(t) remains at 1e-31 -- 1e-31.4'")
s2 = J("stress_2d_sweep.json")["results"]
allS = []
for r in s2:
    print("M=%2d eps=%.1f epochs=%d n_grid=%d  S(0)=%.3e S(0.5)=%.3e S(1)=%.3e   elapsed %.0fs" % (
        r["M"], r["eps"], r["epochs"], r["n_grid"], r["violations"]["0.0"], r["violations"]["0.5"], r["violations"]["1.0"], r["elapsed"]))
    allS += list(r["violations"].values())
S1 = [r["violations"]["1.0"] for r in s2]
print("S(1) range over Table 7: %.3e .. %.3e  =  10^%.3f .. 10^%.3f" % (min(S1), max(S1), np.log10(min(S1)), np.log10(max(S1))))
print("All S values (t=0,0.5,1) in stress_2d_sweep.json: %.3e .. %.3e = 10^%.3f .. 10^%.3f" % (min(allS), max(allS), np.log10(min(allS)), np.log10(max(allS))))
print("Text says 10^-31 .. 10^-31.4  =  %.2e .. %.2e  -> actual min %.2e is %.1fx ABOVE the quoted lower end 1e-31, actual max %.2e is %.1fx above 1e-31 and %.1fx above 10^-31.4" % (
    1e-31, 10 ** -31.4, min(allS), min(allS) / 1e-31, max(allS), max(allS) / 1e-31, max(allS) / 10 ** -31.4))
print("Text in Table-6 paragraph says ~1e-32 .. 1e-31 ; actual %.2e .. %.2e" % (min(v['fourier_spectral'] for v in e2['violations'].values()), max(v['fourier_spectral'] for v in e2['violations'].values())))
print("Table 7 baseline row S(1)=%.3e vs Table 6 Fourier S(1)=%.3e (different re-run, 1.27e-31 vs 1.26e-31 -> Table 6 prints 1.3e-31)" % (s2[0]["violations"]["1.0"], e2["violations"]["1.0"]["fourier_spectral"]))
print("sqrt(S) at S=3.74e-31: %.2e (RMS asymmetry) ; machine eps = %.2e" % (np.sqrt(3.74e-31), np.finfo(float).eps))
