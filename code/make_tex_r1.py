"""Generate numbers_r1.tex (\\N{key} macros) and the table fragments of the revised manuscript from data_r1.
Every experimental number in the manuscript comes from here (no manual transcription).  Missing data never stop the build:
the corresponding macro prints [n/a: key] (bold) and the table row shows '--'.  The list of missing keys is printed at the end.
usage: python make_tex_r1.py"""
import os
import sys
import json
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r1_tables as T

BS = chr(92)
DATA = os.path.join(HERE, "..", "data_r1")
OUT = os.path.join(HERE, "..", "manuscript", "r1_src")
D0 = os.path.join(HERE, "..", "data")
NUM, MISSING = {}, []
NA = BS + "textbf{--}"


# ----------------------------------------------------------------------------- helpers
def runs(pat):
    return [r for r in T.load(pat) if "errors" in r]


def inv_runs(pat):
    return [r for r in T.load(pat) if "eps_hat" in r]


def arr(rs, f):
    if len(rs) == 0:
        raise KeyError('empty group')
    return np.array([f(r) for r in rs], float)


def l2_at(t):
    return lambda r: r["errors"][t]["l2"]


def linf_at(t):
    return lambda r: r["errors"][t]["linf"]


def sci(x, d=2):
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    if round(m, d) >= 10:
        m /= 10
        e += 1
    return f"{m:.{d}f}{BS}times10^{{{e}}}"


def pm(v, d=2):
    """math-mode string: (m +- s) x 10^e (mean of the values, sample SD) or m x 10^e for one value"""
    v = np.asarray(v, float)
    if len(v) == 0:
        return NA
    m = v.mean()
    if len(v) == 1:
        return BS + "ensuremath{" + sci(m, d) + "}"
    s = v.std(ddof=1)
    e = int(np.floor(np.log10(abs(m)))) if m != 0 else 0
    if m != 0 and round(m / 10 ** e, d) >= 10:
        e += 1
    return BS + "ensuremath{" + f"({m / 10 ** e:.{d}f}{BS}pm{s / 10 ** e:.{d}f}){BS}times10^{{{e}}}" + "}"


def mathsci(x, d=2):
    return BS + "ensuremath{" + sci(x, d) + "}"


def put(key, val):
    NUM[key] = val


def put_or_miss(key, thunk):
    try:
        put(key, thunk())
    except Exception as ex:                       # missing data -> visible marker, never a crash
        MISSING.append(key)


def boot_ci(a, b, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    r = [rng.choice(a, len(a)).mean() / rng.choice(b, len(b)).mean() for _ in range(n)]
    return np.percentile(r, [2.5, 97.5])


def paired(a_runs, b_runs, f):
    """per-seed ratios f(a)/f(b) for seeds present in both lists"""
    sa = {r["job"]["seed"]: f(r) for r in a_runs}
    sb = {r["job"]["seed"]: f(r) for r in b_runs}
    return np.array([sa[s] / sb[s] for s in sorted(set(sa) & set(sb))])


def nparams(K, hard):
    out = K - 1 if hard else K + 1
    return 64 + 64 + 64 * 64 + 64 + 64 * out + out


VANILLA_PARAMS = 2 * 64 + 64 + 2 * (64 * 64 + 64) + 64 + 1

# ----------------------------------------------------------------------------- data groups
G = {
    "van": runs("fwd_compat_vanilla_s*"), "soft32": runs("fwd_compat_K32_s*"), "soft32E": runs("fwd_compat_K32_E50_s*"),
    "soft64": runs("fwd_compat_K64_s*"), "soft128": runs("var_K128_uniform_soft"),
    "h32": runs("fh_K32_s*"), "h32E": runs("fh_K32_E50_s*"), "h64": runs("fh_K64_s*"), "h128": runs("fh_K128_s*"),
    "lbvan": runs("lb_vanilla_s*"), "lbh32": runs("lb_hard_K32_s*"), "lbh64": runs("lb_hard_K64_s*"),
    "b2van": runs("b2_vanilla_s*"), "b2hard": runs("b2_hard_K32_s*"), "b2soft": runs("b2_soft_K32_s*"),
    "b2lbvan": runs("lb_b2_vanilla_s*"), "b2lbh": runs("lb_b2_hard_K32_s*"),
    "b3van": runs("b3_vanilla_s*"), "b3h32": runs("b3_K32_E0_s*"), "b3h32E": runs("b3_K32_E50_s*"), "b3h64": runs("b3_K64_E0_s*"), "b3h64E": runs("b3_K64_E50_s*"),
    "fix32": runs("vs_K32_bcfix_s*"), "fix64": runs("vs_K64_bcfix_s*"), "fix128": runs("vs_K128_bcfix_s*"),
    "k64x2": runs("kh_K64_x2"), "k128x2": runs("kh_K128_x2"), "k256": runs("kh_K256_x1"),
    "o32": runs("fwd_orig_K32_s*"), "o64": runs("fwd_orig_K64_s*"),
    "v64cgl_soft": runs("var_K64_cgl_soft"), "v64u_hard": runs("var_K64_uniform_hard"), "v64cgl_hard": runs("var_K64_cgl_hard"),
    "v32cgl_soft": runs("var_K32_cgl_soft"),
}
timing = {k: runs(f"timing_{k}_s42") for k in ("vanilla", "soft_K32", "hard_K32", "hard_K64")}
CL = json.load(open(os.path.join(DATA, "classical.json")))


def cheb(ic, eps, N, dt, key="l2"):
    for c in CL[f"{ic}_eps{eps}"]["cheb"] if f"{ic}_eps{eps}" in CL else []:
        if c["N"] == N and abs(c["dt"] - dt) < 1e-12:
            return c[key]
    raise KeyError


def cheb3(N, dt, key="l2"):
    for c in CL["compat_eps0.05_T5"]["cheb"]:
        if c["N"] == N and abs(c["dt"] - dt) < 1e-12:
            return c[key]
    raise KeyError


def tstr(rs, tkey, which="train_time"):
    """time string: clean timing run if available, else mean of the main runs marked with a dagger"""
    t = timing.get(tkey)
    if t:
        return f"{t[0][which]:.0f}"
    if rs:
        return f"{np.mean(arr(rs, lambda r: r[which])):.0f}" + "$^\\dagger$"
    return NA


# ----------------------------------------------------------------------------- numbers
for ic in ("compat", "orig"):
    for eps in ("0.01", "0.1"):
        put_or_miss(f"ode.{ic}.{eps}.l2", lambda ic=ic, eps=eps: mathsci(CL[f"{ic}_eps{eps}"]["ode"]["1.0"]["l2"]))

put_or_miss("b1.cheb80.l2", lambda: mathsci(cheb("compat", "0.01", 80, 5e-4)))
put_or_miss("b1.cheb80.dt1e-4.l2", lambda: mathsci(cheb("compat", "0.01", 80, 1e-4)))
put_or_miss("b1.cheb80.time", lambda: f"{cheb('compat', '0.01', 80, 5e-4, 'time_median'):.1f}")
put_or_miss("b2.cheb80.l2", lambda: mathsci(cheb("compat", "0.1", 80, 5e-4)))
put_or_miss("b3.cheb.l2", lambda: mathsci(cheb3(96, 1e-3)))

l2 = l2_at("1.0")
put_or_miss("b1.fh_K32.l2", lambda: pm(arr(G["h32"], l2)))
put_or_miss("b1.van.l2", lambda: pm(arr(G["van"], l2)))
put_or_miss("b1.lb_van.l2", lambda: pm(arr(G["lbvan"], l2)))
put_or_miss("b1.lb_hard_K32.l2", lambda: pm(arr(G["lbh32"], l2)))
put_or_miss("b1.ratio.van_fhK32", lambda: f"{arr(G['van'], l2).mean() / arr(G['h32'], l2).mean():.0f}")
put_or_miss("b1.ratio.van_fhK32.ci", lambda: "[%.0f, %.0f]" % tuple(boot_ci(arr(G["van"], l2), arr(G["h32"], l2))))
put_or_miss("b1.ratio.van_fhK32.range", lambda: "%.0f--%.0f" % (paired(G["van"], G["h32"], l2).min(), paired(G["van"], G["h32"], l2).max()))
put_or_miss("b1.ratio.van_fhK32.linf", lambda: f"{arr(G['van'], linf_at('1.0')).mean() / arr(G['h32'], linf_at('1.0')).mean():.0f}")
put_or_miss("b1.ratio.lbvan_fhK32", lambda: f"{arr(G['lbvan'], l2).mean() / arr(G['h32'], l2).mean():.0f}")
put_or_miss("b1.ratio.fh_K32_cheb80", lambda: f"{arr(G['h32'], l2).mean() / cheb('compat', '0.01', 80, 5e-4):.0f}")
put_or_miss("b1.ratio.ode_fh_K32", lambda: f"{CL['compat_eps0.01']['ode']['1.0']['l2'] / arr(G['h32'], l2).mean():.1f}")
put_or_miss("b1.speed.fh_K32_cheb80", lambda: f"{np.mean(arr(G['h32'], lambda r: r['train_time'])) / cheb('compat', '0.01', 80, 5e-4, 'time_median'):.0f}")
put_or_miss("b1.speed.van_fhK32", lambda: f"{np.mean(arr(G['van'], lambda r: r['train_time'])) / np.mean(arr(G['h32'], lambda r: r['train_time'])):.1f}")
put_or_miss("b2.hard_K32.l2", lambda: pm(arr(G["b2hard"], l2)))
put_or_miss("b2.van.l2", lambda: pm(arr(G["b2van"], l2)))
put_or_miss("b2.ratio.ode_hard", lambda: f"{CL['compat_eps0.1']['ode']['1.0']['l2'] / arr(G['b2hard'], l2).mean():.0f}")
l25 = l2_at("5.0")
put_or_miss("b3.K32_E0.l2", lambda: pm(arr(G["b3h32"], l25)))
put_or_miss("b3.K64_E0.l2", lambda: pm(arr(G["b3h64"], l25)))
put_or_miss("b3.van.l2", lambda: pm(arr(G["b3van"], l25)))

# K scaling
put_or_miss("kk.soft64.l2", lambda: pm(arr(G["soft64"], l2)))
put_or_miss("kk.hard64.l2", lambda: pm(arr(G["h64"], l2)))
put_or_miss("kk.soft128.l2", lambda: pm(arr(G["soft128"], l2)))
put_or_miss("kk.hard128.l2", lambda: pm(arr(G["h128"], l2)))
put_or_miss("kk.fix64.l2", lambda: pm(arr(G["fix64"], l2)))
put_or_miss("kk.hard32.l2", lambda: pm(arr(G["h32"], l2)))
put_or_miss("kk.soft64.wall", lambda: f"{max(r['balance']['max_abs_ux_wall'] for r in G['soft64']):.1f}")
put_or_miss("kk.cgl_soft_gain", lambda: f"{G['soft64'][0]['errors']['1.0']['l2'] / G['v64cgl_soft'][0]['errors']['1.0']['l2']:.0f}")

# energy
allh = G["h32"] + G["h32E"]
put_or_miss("en.n_runs", lambda: str(len(allh)) if allh else (_ for _ in ()).throw(KeyError))
put_or_miss("en.min_step", lambda: mathsci(min(-r["energy"]["dt0.02"]["max_step"] for r in allh)))
put_or_miss("en.max_diff", lambda: mathsci(max(np.max(np.abs(np.array(a["E_traj_dt0.02"]) - np.array(b["E_traj_dt0.02"]))) for a in G["h32"] for b in G["h32E"] if a["job"]["seed"] == b["job"]["seed"])))
put_or_miss("en.nonzero", lambda: "%d" % max(r["penalty"]["n_epochs_nonzero"] for r in G["h32E"]))
put_or_miss("en.van_n", lambda: str(len(G["van"])) if G["van"] else (_ for _ in ()).throw(KeyError))
put_or_miss("en.van_viol", lambda: str(sum(r["energy"]["dt0.02"]["n_increases"] for r in G["van"])))

# maximum bound principle
allspec = G["soft32"] + G["soft32E"] + G["soft64"] + G["soft128"] + G["h32"] + G["h32E"] + G["h64"] + G["h128"]
put_or_miss("mbp.max_b1", lambda: "%.4f" % max(r["max_abs_u"] for r in allspec))
put_or_miss("mbp.max_hard", lambda: "%.4f" % max(r["max_abs_u"] for r in G["h32"] + G["h32E"] + G["h64"] + G["h128"]))
put_or_miss("mbp.max_van", lambda: "%.4f" % max(r["max_abs_u"] for r in G["van"]))
put_or_miss("mbp.excess", lambda: mathsci(max(r["max_abs_u"] for r in allspec) - 1.0, 1))
put_or_miss("mbp.excess_all", lambda: "%.1f\\%%" % (100 * (max(r["max_abs_u"] for r in allspec + G["van"]) - 1.0)))

# balance
put_or_miss("bal.id_err", lambda: mathsci(max(r["balance"]["identity_rel_err"] for r in allspec if "balance" in r), 1))
put_or_miss("bal.int_r_hard", lambda: mathsci(np.mean(arr(G["h32E"], lambda r: r["balance"]["int_half_r2"]))))
put_or_miss("bal.dE_hard", lambda: "%.3f" % np.mean(arr(G["h32E"], lambda r: r["balance"]["E_first"] - r["balance"]["E_last"])))
put_or_miss("bal.frac_hard", lambda: "%.3f\\%%" % (100 * np.mean(arr(G["h32E"], lambda r: r["balance"]["int_half_r2"])) / np.mean(arr(G["h32E"], lambda r: r["balance"]["E_first"] - r["balance"]["E_last"]))))
put_or_miss("bal.int_b_soft", lambda: mathsci(np.mean(arr(G["soft64"], lambda r: r["balance"]["int_abs_bterm"]))))

# identifiability (classical inverse)
ID = json.load(open(os.path.join(DATA, "ident.json")))
put("ident.S001", "%.3f" % ID["compat_0.01"]["sensitivity_norm_dlogeps"])
put("ident.S005", "%.2f" % ID["compat_0.05"]["sensitivity_norm_dlogeps"])
put("ident.S01", "%.2f" % ID["compat_0.1"]["sensitivity_norm_dlogeps"])
put("ident.se001", "%.0f\\%%" % ID["compat_0.01"]["runs"]["0.01"]["fisher_rel_se_pct"])
put("ident.cl001", "%.0f\\%%" % ID["compat_0.01"]["runs"]["0.01"]["rel_rmse_pct"])
put("ident.cl005", "%.1f\\%%" % ID["compat_0.05"]["runs"]["0.01"]["rel_rmse_pct"])

# inverse
adam = inv_runs("inv_adam_*")
if adam:
    rr = arr(adam, lambda r: r["ratio_to_init"])
    put("inv.ratio_mean", "%.2f" % rr.mean())
    put("inv.ratio_range", "%.2f--%.2f" % (rr.min(), rr.max()))
    put("inv.steps", "%.0f" % (abs(np.log(rr.mean())) / 0.02))
else:
    MISSING.append("inv.ratio_*")
for e, key in ((0.05, "005"), (0.1, "01")):
    put_or_miss(f"inv.vp_err_{key}", lambda e=e: "%.1f\\%%" % inv_runs(f"inv_vp_e{e}_ld2000")[0]["rel_err_pct"])

# 2D (unchanged experiments)
e2 = json.load(open(os.path.join(D0, "experiment_2d.json")))["violations"]
m2 = json.load(open(os.path.join(D0, "symmetry_2d_multiseed.json")))["results"]
s2 = json.load(open(os.path.join(D0, "stress_2d_sweep.json")))["results"]
F = [v["fourier_spectral"] for v in e2.values()]
put("s2d.fmin", mathsci(min(F))); put("s2d.fmax", mathsci(max(F)))
gr = [r["violations"]["1.0"] / r["violations"]["0.0"] for r in m2]
put("s2d.growth", "%.1f" % (np.mean([r["violations"]["1.0"] for r in m2]) / np.mean([r["violations"]["0.0"] for r in m2])))
put("s2d.growth_range", "%.1f--%.1f" % (min(gr), max(gr)))
S1 = [r["violations"]["1.0"] for r in s2]
put("s2d.smin", mathsci(min(S1))); put("s2d.smax", mathsci(max(S1)))


# ----------------------------------------------------------------------------- tables
def table(label, caption, colspec, header, body, extra="", size="small"):
    t = [BS + "begin{table}[H]", BS + "caption{" + caption + BS + "label{" + label + "}}", BS + size, BS + "begin{tabularx}{" + BS + "textwidth}{" + colspec + "}", BS + "toprule", header, BS + "midrule"]
    t += body + [BS + "bottomrule", BS + "end{tabularx}", extra, BS + "end{table}"]
    return "\n".join(t) + "\n"


def row(*cells):
    return " & ".join(cells) + BS + BS


def maxu(rs):
    return "%.4f" % max(r["max_abs_u"] for r in rs) if rs else NA


def viol(rs, k="dt0.02"):
    return "%d/%d" % (sum(r["energy"][k]["n_increases"] for r in rs), len(rs) * 50) if rs else NA


def meth_row(name, rs, par, tkey, tfinal="1.0"):
    if not rs:
        return row(name, par, NA, NA, NA, NA, NA, NA)
    f = lambda r: r["errors"][tfinal]
    return row(name, par, pm(arr(rs, lambda r: f(r)["l2"])), pm(arr(rs, lambda r: f(r)["linf"])), str(len(rs)), tstr(rs, tkey), maxu(rs), viol(rs))


def ode_row(ic, eps):
    return row("Reaction-only ODE ($\\varepsilon=0$)", "0", mathsci(CL[f"{ic}_eps{eps}"]["ode"]["1.0"]["l2"]), mathsci(CL[f"{ic}_eps{eps}"]["ode"]["1.0"]["l2"] * 0 + 0) if False else NA, "1", "--", "--", "--")


def cheb_row(ic, eps, N, dt):
    try:
        return row(f"Chebyshev collocation, $N={N}$, $\\delta t={dt:g}$", str(N + 1), mathsci(cheb(ic, eps, N, dt)), mathsci(cheb(ic, eps, N, dt, "linf")), "1", f"{cheb(ic, eps, N, dt, 'time_median'):.2f}", "--", "--")
    except KeyError:
        return row(f"Chebyshev collocation, $N={N}$, $\\delta t={dt:g}$", NA, NA, NA, NA, NA, NA, NA)


hdr = row(BS + "textbf{Method}", BS + "textbf{Params}", "$L^2$ error", "$L^\\infty$ error", "$n$", "Time (s)", "$\\max|u|$", "Energy incr.")
body = [BS + "multicolumn{8}{l}{" + BS + "textit{B1: $\\varepsilon=0.01$, $T=1$ (error at $t=1$)}}" + BS + BS]
body += [cheb_row("compat", "0.01", 80, 5e-4), cheb_row("compat", "0.01", 80, 1e-4), ode_row("compat", "0.01"),
         meth_row("Vanilla PINN", G["van"], "8577", "vanilla"),
         meth_row("Vanilla PINN + L-BFGS", G["lbvan"], "8577", "vanilla"),
         meth_row("Spectral-PINN, as submitted, $K=32$", G["soft32"], str(nparams(32, False)), "soft_K32"),
         meth_row("\\quad + energy penalty ($\\lambda_E=50$)", G["soft32E"], str(nparams(32, False)), "soft_K32"),
         meth_row("Spectral-PINN, as submitted, $K=64$", G["soft64"], str(nparams(64, False)), "soft_K32"),
         meth_row("Spectral-PINN, revised, $K=32$", G["h32"], str(nparams(32, True)), "hard_K32"),
         meth_row("\\quad + energy penalty ($\\lambda_E=50$)", G["h32E"], str(nparams(32, True)), "hard_K32"),
         meth_row("\\quad + L-BFGS", G["lbh32"], str(nparams(32, True)), "hard_K32"),
         meth_row("Spectral-PINN, revised, $K=64$", G["h64"], str(nparams(64, True)), "hard_K64"),
         meth_row("Spectral-PINN, revised, $K=128$", G["h128"], str(nparams(128, True)), "hard_K64")]
body += [BS + "midrule", BS + "multicolumn{8}{l}{" + BS + "textit{B2: $\\varepsilon=0.1$, $T=1$}}" + BS + BS,
         cheb_row("compat", "0.1", 80, 5e-4), ode_row("compat", "0.1"),
         meth_row("Vanilla PINN", G["b2van"], "8577", "vanilla"), meth_row("Vanilla PINN + L-BFGS", G["b2lbvan"], "8577", "vanilla"),
         meth_row("Spectral-PINN, as submitted, $K=32$", G["b2soft"], str(nparams(32, False)), "soft_K32"),
         meth_row("Spectral-PINN, revised, $K=32$", G["b2hard"], str(nparams(32, True)), "hard_K32"),
         meth_row("\\quad + L-BFGS", G["b2lbh"], str(nparams(32, True)), "hard_K32")]
b3rows = [BS + "midrule", BS + "multicolumn{8}{l}{" + BS + "textit{B3: $\\varepsilon=0.05$, $T=5$ (error at $t=5$)}}" + BS + BS]
try:
    b3rows.append(row("Chebyshev collocation, $N=96$, $\\delta t=10^{-3}$", "97", mathsci(cheb3(96, 1e-3)), mathsci(cheb3(96, 1e-3, "linf")), "1", f"{cheb3(96, 1e-3, 'time_median'):.2f}", "--", "--"))
except KeyError:
    pass
b3rows += [row("Reaction-only ODE ($\\varepsilon=0$)", "0", mathsci(CL["compat_eps0.05_T5"]["ode"]["5.0"]["l2"]), NA, "1", "--", "--", "--")]
b3rows += [meth_row("Vanilla PINN", G["b3van"], "8577", "vanilla", "5.0"),
           meth_row("Spectral-PINN, revised, $K=32$", G["b3h32"], str(nparams(32, True)), "hard_K32", "5.0"),
           meth_row("\\quad + energy penalty", G["b3h32E"], str(nparams(32, True)), "hard_K32", "5.0"),
           meth_row("Spectral-PINN, revised, $K=64$", G["b3h64"], str(nparams(64, True)), "hard_K64", "5.0"),
           meth_row("\\quad + energy penalty", G["b3h64E"], str(nparams(64, True)), "hard_K64", "5.0")]
body += b3rows
cap = ("\\rev{Forward problem on the three benchmarks. Errors against the certified reference on the common 200-point grid, mean $\\pm$ SD over $n$ seeds; ``Energy incr.'' = number of increases of the energy over the intervals of a grid of spacing $0.02$ (total over seeds, $50T$ intervals per run); $\\max|u|$ is the largest value over all runs of the row. "
       "Time is the training time of one run on an idle single core (rows with $^\\dagger$: mean over the main runs, which were executed concurrently and are not a clean timing); the classical time is the median of three solves. The vanilla PINN has 8000 Adam epochs, the Spectral-PINN 6000 (+1000 warm-up) for $T=1$.}")
open(os.path.join(OUT, "tab_forward.tex"), "w", encoding="utf-8").write(table("tab:forward", cap, "Xccccccc", hdr, body, size="footnotesize"))

# ---- K table
hdrK = row(BS + "textbf{Variant}", "$K$", "Budget", "$L^2$ error", "max $|u_x(\\pm1)|$", "$n$", "Time (s)")


def krow(name, K, bud, rs):
    if not rs:
        return row(name, str(K), bud, NA, NA, "0", NA)
    wall = "%.2f" % max(r["balance"]["max_abs_ux_wall"] for r in rs) if "balance" in rs[0] else NA
    return row(name, str(K), bud, pm(arr(rs, l2)), wall, str(len(rs)), f"{np.mean(arr(rs, lambda r: r['train_time'])):.0f}$^\\dagger$")


bodyK = [krow("soft penalty, equispaced grid (as submitted)", 32, "$1\\times$", G["soft32"]), krow("", 64, "$1\\times$", G["soft64"]), krow("", 128, "$1\\times$", G["soft128"]),
         krow("soft penalty, Gauss--Lobatto grid", 32, "$1\\times$", G["v32cgl_soft"]), krow("", 64, "$1\\times$", G["v64cgl_soft"]),
         krow("soft penalty rescaled ($\\lambda_{\\rm BC}=K^4$)", 32, "$1\\times$", G["fix32"]), krow("", 64, "$1\\times$", G["fix64"]), krow("", 128, "$1\\times$", G["fix128"]),
         BS + "midrule",
         krow("exact Neumann basis, equispaced grid", 64, "$1\\times$", G["v64u_hard"]),
         krow("exact Neumann basis, Gauss--Lobatto grid (revised)", 32, "$1\\times$", G["h32"]), krow("", 64, "$1\\times$", G["h64"]), krow("", 128, "$1\\times$", G["h128"]),
         krow("", 64, "$2\\times$", G["k64x2"]), krow("", 128, "$2\\times$", G["k128x2"]), krow("", 256, "$1\\times$", G["k256"])]
capK = ("\\rev{Effect of the mode number $K$ for different treatments of the Neumann condition and of the collocation grid (B1, $L^2$ error at $t=1$). Budget is the number of epochs relative to $6000+1000$ warm-up epochs. Single-seed rows have $n=1$ (SD not defined). The wall slope is the maximum of $|u_x(\\pm1,t)|$ over the trained trajectory. Times marked $^\\dagger$ come from concurrently executed runs.}")
open(os.path.join(OUT, "tab_K.tex"), "w", encoding="utf-8").write(table("tab:K", capK, "p{4.8cm}cccccc", hdrK, bodyK))

# ---- energy table
hdrE = row("Setting", "$n$", "$L^2$ error at $t=1$", "increases (0.02 / 0.002 / $2\\times10^{-4}$)", "penalty $>0$ in epochs")
bodyE = []
for lam, pat in ((0, "fh_K32_s*"), (1, "ah_lamE1_s*"), (10, "ah_lamE10_s*"), (50, "fh_K32_E50_s*"), (500, "ah_lamE500_s*"), (5000, "ah_lamE5000_s*")):
    rs = runs(pat)
    if rs:
        inc = " / ".join(str(sum(r["energy"][k]["n_increases"] for r in rs)) for k in ("dt0.02", "dt0.002", "dt0.0002"))
        nz = "%d" % max(r["penalty"]["n_epochs_nonzero"] for r in rs)
        bodyE.append(row("$\\lambda_E=" + str(lam) + "$", str(len(rs)), pm(arr(rs, l2)), inc, nz))
    else:
        bodyE.append(row("$\\lambda_E=" + str(lam) + "$", "0", NA, NA, NA))
bodyE.append(BS + "midrule")
for name, p0, p1 in (("300 epochs, warm-up 100", "rh_short_E0_s*", "rh_short_E50_s*"), ("no warm-up", "rh_nowarm_E0_s*", "rh_nowarm_E50_s*"), ("learning rate $10^{-4}$, 1500 epochs", "rh_pde0_E0_s*", "rh_pde0_E50_s*")):
    for lam, pat in ((0, p0), (50, p1)):
        rs = runs(pat)
        if rs:
            inc = " / ".join(str(sum(r["energy"][k]["n_increases"] for r in rs)) for k in ("dt0.02", "dt0.002", "dt0.0002"))
            bodyE.append(row(f"{name}, $\\lambda_E={lam}$", str(len(rs)), pm(arr(rs, l2)), inc, "%d" % max(r["penalty"]["n_epochs_nonzero"] for r in rs)))
        else:
            bodyE.append(row(f"{name}, $\\lambda_E={lam}$", "0", NA, NA, NA))
capE = ("\\rev{Energy penalty weight and regime (revised Spectral-PINN, $K=32$, B1). ``increases'' is the total number of energy increases (more than $10^{-10}$) over all seeds on grids of spacing $0.02$, $0.002$, $2\\times10^{-4}$; ``penalty $>0$'' is the largest number, over the seeds, of training epochs at which $\\mathcal L_E>0$ (out of the total number of epochs).}")
open(os.path.join(OUT, "tab_energy.tex"), "w", encoding="utf-8").write(table("tab:energy", capE, "p{4.6cm}cccc", hdrE, bodyE))

# ---- balance table
hdrB = row(BS + "textbf{Network}", "$n$", "max $|u_x(\\pm1)|$", "$\\frac12\\!\\int\\!\\|r\\|^2$", "$\\int|\\beta|$", "$E(0)-E(T)$", "id. error")
bodyB = []
for name, rs in (("revised $K=32$", G["h32"]), ("revised $K=32$, $\\lambda_E=50$", G["h32E"]), ("revised $K=64$", G["h64"]), ("revised $K=128$", G["h128"]),
                 ("as submitted $K=64$", G["soft64"]), ("revised, B2", G["b2hard"]), ("revised $K=64$, B3", G["b3h64"])):
    if rs and "balance" in rs[0]:
        b = lambda k: arr(rs, lambda r: r["balance"][k])
        bodyB.append(row(name, str(len(rs)), "%.1e" % b("max_abs_ux_wall").max(), mathsci(b("int_half_r2").mean()), mathsci(b("int_abs_bterm").mean()) if b("int_abs_bterm").mean() > 0 else "0",
                         "%.3f" % (b("E_first") - b("E_last")).mean(), mathsci(b("identity_rel_err").max(), 1)))
    else:
        bodyB.append(row(name, "0", NA, NA, NA, NA, NA))
capB = ("\\rev{A-posteriori energy balance of the trained trajectories (Section~\\ref{sec:balance}): wall slope, time-integrated residual $\\frac12\\int_0^T\\|r\\|^2dt$, time-integrated boundary defect, actual energy decrease, and relative error of the identity~\\eqref{eq:balance-main} checked against the autograd derivative of $E(t)$ (maximum over seeds). Means over seeds.}")
open(os.path.join(OUT, "tab_balance.tex"), "w", encoding="utf-8").write(table("tab:balance", capB, "Xcccccc", hdrB, bodyB))

# ---- inverse tables
def invrow(name, rs, f=lambda r: r["ratio_to_true"]):
    return row(name, str(len(rs)), "%.2f" % np.mean(arr(rs, f)) if rs else NA)


hdrI = row("$\\varepsilon_{\\rm init}/\\varepsilon_{\\rm true}$", "variant", "$\\varepsilon_{\\rm true}$", "$\\hat\\varepsilon/\\varepsilon_{\\rm true}$", "$\\hat\\varepsilon/\\varepsilon_{\\rm init}$", "data misfit (RMS)")
bodyI = []
for e in (0.05, 0.1):
    for init in (0.5, 1.0, 1.5, 2.0, 4.0):
        for var in ("soft", "hard"):
            rs = inv_runs(f"inv_adam_{var}_e{e}_i{init}")
            if rs:
                r = rs[0]
                bodyI.append(row(f"{init:g}", "as submitted" if var == "soft" else "revised", f"{e:g}", "%.2f" % r["ratio_to_true"], "%.2f" % r["ratio_to_init"], "%.1e" % r["data_rms_misfit"]))
            else:
                bodyI.append(row(f"{init:g}", var, f"{e:g}", NA, NA, NA))
capI1 = ("\\rev{Joint Adam training of the network and $\\log\\varepsilon$ (as in the first version): recovered $\\hat\\varepsilon$ for different starting values $\\varepsilon_{\\rm init}$ (noise-free data, 5 snapshots, 40 points, $\\lambda_{\\rm data}=2$, one run each). The ratio to the starting value is nearly constant, and data and truth have little influence.}")
t1 = table("tab:inverse1", capI1, "ccccXc".replace("X", "c"), hdrI, bodyI)

hdrJ = row("Setting", "$n$", "$\\varepsilon_{\\rm true}$", "rel. error of $\\hat\\varepsilon$ (\\%)", "data misfit (RMS)")
bodyJ = []
for e in (0.05, 0.1):
    for ld in (20, 200, 2000, 20000):
        rs = inv_runs(f"inv_vp_e{e}_ld{ld}")
        bodyJ.append(row(f"variable projection, $\\lambda_{{\\rm data}}={ld}$", str(len(rs)), f"{e:g}", "%.1f" % rs[0]["rel_err_pct"] if rs else NA, "%.1e" % rs[0]["data_rms_misfit"] if rs else NA))
bodyJ.append(BS + "midrule")
for tag, name in (("vp_e0.05_nx10", "10 observation points"), ("vp_e0.05_nx160", "160 observation points"), ("vp_e0.05_nt3", "3 snapshots"), ("vp_e0.05_nt20", "20 snapshots"),
                  ("vp_e0.05_lpde5", "$\\lambda_{\\rm PDE}=5$"), ("vp_e0.05_lpde500", "$\\lambda_{\\rm PDE}=500$"), ("vp_e0.05_nophaseA", "no data-fitting phase"), ("vp_e0.05_K64", "$K=64$")):
    rs = inv_runs("inv_" + tag)
    bodyJ.append(row(name + ", $\\lambda_{\\rm data}=2000$", str(len(rs)), "0.05", "%.1f" % rs[0]["rel_err_pct"] if rs else NA, "%.1e" % rs[0]["data_rms_misfit"] if rs else NA))
capI2 = ("\\rev{Variable-projection estimator: noise-free data, effect of the data weight (upper part) and of the observation design and loss weights (lower part, $\\varepsilon_{\\rm true}=0.05$).}")
t2 = table("tab:inverse2", capI2, "Xcccc", hdrJ, bodyJ)

hdrN = row("$\\varepsilon_{\\rm true}$", "$\\sigma$", "neural (variable projection): mean $|$rel. err.$|$ (\\%)", "classical LS: RMSE (\\%)", "Cram\\'er--Rao (\\%)")
bodyN = []
for e in (0.05, 0.1):
    for nz in (0.01, 0.05):
        rs = inv_runs(f"inv_vpnoise_e{e}_n{nz}_r*")
        cl = ID[f"compat_{e}"]["runs"][str(nz)]
        bodyN.append(row(f"{e:g}", f"{nz:g}", "%.1f (n=%d)" % (np.mean(arr(rs, lambda r: r["rel_err_pct"])), len(rs)) if rs else NA, "%.1f" % cl["rel_rmse_pct"], "%.1f" % cl["fisher_rel_se_pct"]))
for nz in (0.0, 0.01):
    rs = inv_runs(f"inv_vp_e0.01_n{nz}")
    cl = ID["compat_0.01"]["runs"][str(nz)]
    bodyN.append(row("0.01", f"{nz:g}", "%.1f" % rs[0]["rel_err_pct"] if rs else NA, "%.1f" % cl["rel_rmse_pct"], "%.1f" % cl["fisher_rel_se_pct"]))
capI3 = ("\\rev{Inverse problem with noisy data: neural variable-projection estimator (5 noise realizations; $\\lambda_{\\rm data}=2000$) against the classical least-squares estimator with the Chebyshev solver (20 realizations) and the Cram\\'er--Rao bound $\\sigma/\\|\\partial_{\\log\\varepsilon}\\mathbf d\\|$. At $\\varepsilon=0.01$ the problem is not identifiable for $\\sigma\\ge0.01$.}")
t3 = table("tab:inverse3", capI3, "ccXcc", hdrN, bodyN)
open(os.path.join(OUT, "tab_inverse.tex"), "w", encoding="utf-8").write(t1 + "\n" + t2 + "\n" + t3)

# ---- 2D tables (data unchanged from the first version)
bodyS = []
vv = {t: [r["violations"][t] for r in m2] for t in ("0.0", "0.2", "0.5", "1.0")}
for t in ("0.0", "0.2", "0.5", "1.0"):
    bodyS.append(row(t, mathsci(e2[t]["fourier_spectral"]), pm(vv[t])))
capS = ("Reflection-symmetry violation $S(t)$ (Section~\\ref{sec:2d}) for the Fourier-Spectral-PINN (single run; architecturally exact regardless of seed) versus a vanilla PINN with the same periodic input encoding (mean $\\pm$ SD over 5 seeds). Data of the first version, unchanged.")
tS = table("tab:2d", capS, "ccc", row("$t$", "Fourier-Spectral-PINN", "Vanilla PINN (5 seeds)"), bodyS)
bodyS2 = [row(str(r["M"]), f"{r['eps']:g}", str(r["epochs"]), mathsci(r["violations"]["1.0"])) for r in s2]
capS2 = ("Fourier-Spectral-PINN reflection-symmetry violation $S(1)$ at larger $M$ and smaller $\\varepsilon$ (data of the first version, unchanged).")
tS2 = table("tab:stress-2d", capS2, "cccc", row("$M$", "$\\varepsilon$", "epochs", "$S(1)$"), bodyS2)
open(os.path.join(OUT, "tab_2d.tex"), "w", encoding="utf-8").write(tS + "\n" + tS2)

# ---- appendix: original initial condition
bodyO = []
for name, rs, rc in (("as submitted, $K=32$", G["o32"], G["soft32"]), ("as submitted, $K=64$", G["o64"], G["soft64"])):
    wall = lambda g: ("%.2f" % max(r["balance"]["max_abs_ux_wall"] for r in g)) if g and "balance" in g[0] else NA
    bodyO.append(row(name, str(len(rs)), pm(arr(rs, l2)) if rs else NA, wall(rs), str(len(rc)), pm(arr(rc, l2)) if rc else NA, wall(rc)))
capO = (r"\rev{Original initial condition $x^2\cos(\pi x)$ versus the compatible condition~\eqref{eq:IC}: $L^2$ error at $t=1$ ($\varepsilon=0.01$) of the as-submitted Spectral-PINN, mean $\pm$ SD over $n$ seeds, and the largest wall slope $|u_x(\pm1,t)|$.}")
open(os.path.join(OUT, "tab_orig.tex"), "w", encoding="utf-8").write(table("tab:orig", capO, "Xcccccc", row("Variant", "$n$", "original: $L^2$", "wall slope", "$n$", "compatible: $L^2$", "wall slope"), bodyO))

# ----------------------------------------------------------------------------- numbers file
lines = ["% generated by make_tex_r1.py -- do not edit", BS + "makeatletter",
         BS + "newcommand{" + BS + "N}[1]{" + BS + "@ifundefined{num@#1}{" + BS + "textbf{[n/a: " + BS + "detokenize{#1}]}}{" + BS + "csname num@#1" + BS + "endcsname}}", BS + "makeatother"]
for k, v in sorted(NUM.items()):
    lines.append(BS + "expandafter" + BS + "def" + BS + "csname num@" + k + BS + "endcsname{" + v + "}")
open(os.path.join(OUT, "numbers_r1.tex"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
json.dump(NUM, open(os.path.join(DATA, "numbers.json"), "w"), indent=1)
print(len(NUM), "numbers written;", len(MISSING), "missing:", sorted(set(MISSING)))
