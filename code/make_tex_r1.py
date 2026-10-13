"""Generate numbers_r1.tex (\\N{key} macros) and the tables tab_*.tex of the revised manuscript from data_r1.
Every result of the training runs and classical baselines of revision 1 quoted in the manuscript is produced here (no manual transcription).

usage: python make_tex_r1.py
  -> manuscript/r1_src/numbers_r1.tex, manuscript/r1_src/tab_{forward,K,energy,balance,inverse,2d,orig}.tex, data_r1/numbers.json
  prints the keys used in manuscript/r1_src/*.tex that could not be generated (missing data) and a short data report."""
import os
import re
import sys
import json
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r1_tables as T

BS = chr(92)
ROOT = os.path.join(HERE, "..")
OUT_DIR = os.path.join(ROOT, "manuscript", "r1_src")
os.makedirs(OUT_DIR, exist_ok=True)
RUNS = T.RUNS
NUM = {}
REPORT = []
SEEDS = [42, 123, 7, 2024, 99]


# ------------------------------------------------------------------------------------------------ formatting
def _exp(x):
    return int(np.floor(np.log10(abs(x))))


def tex_sci(x, d=2):
    """3.45e-4 -> 3.45\\times10^{-4} (no $)"""
    if x == 0:
        return "0"
    e = _exp(x)
    m = x / 10 ** e
    if round(abs(m), d) >= 10:
        m /= 10; e += 1
    if e == 0:
        return f"{m:.{d}f}"
    return f"{m:.{d}f}{BS}times10^{{{e}}}"


def sci_pm(mean, sd, d=2):
    """(m +- s) x 10^e with the exponent of the mean (no $)."""
    if mean == 0:
        return "0"
    e = _exp(mean)
    if round(abs(mean) / 10 ** e, d) >= 10:
        e += 1
    if e == 0:
        return f"{mean:.{d}f}{BS}pm{sd:.{d}f}"
    return f"({mean / 10 ** e:.{d}f}{BS}pm{sd / 10 ** e:.{d}f}){BS}times10^{{{e}}}"


def num(x, sig=2):
    """plain number with `sig` significant digits (no exponent for 1e-3 <= |x| < 1e5)"""
    if x == 0:
        return "0"
    e = _exp(x)
    if e < -3 or e >= 5:
        return tex_sci(x, sig - 1)
    dec = max(0, sig - 1 - e)
    if dec == 0:
        x = round(x, sig - 1 - e)          # 756 -> 760 for two significant digits
    s = f"{round(x, dec):.{dec}f}"
    if abs(x) >= 1e4:
        s = s[:-3] + BS + "," + s[-3:]
    return s


def pct(x, sig=2):
    return num(x, sig) + BS + "%"


def m(s):
    return BS + "ensuremath{" + s + "}"


def stat(v):
    v = np.asarray(v, float)
    return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0, "n": int(len(v)), "min": float(v.min()), "max": float(v.max())}


def f_stat(v, d=2):
    s = stat(v)
    return m(sci_pm(s["mean"], s["sd"], d)) if s["n"] > 1 else m(tex_sci(s["mean"], d))


def f_sci(x, d=2):
    return m(tex_sci(x, d))


def put(key, val):
    NUM[key] = val


# ------------------------------------------------------------------------------------------------ data access
_cache = {}


def runs(pat):
    if pat not in _cache:
        out = []
        for f in sorted(glob.glob(os.path.join(RUNS, pat + ".json"))):
            r = json.load(open(f))
            if "errors" in r or "eps_hat" in r:
                out.append(r)
        _cache[pat] = out
    return _cache[pat]


def run1(tag):
    rr = runs(tag)
    return rr[0] if rr else None


def seed(r):
    return r["job"].get("seed", 42)


def l2(rr, t="1.0", norm="l2"):
    return [r["errors"][t][norm] for r in rr]


def by_seed(rr):
    return {seed(r): r for r in rr}


def jload(name):
    p = os.path.join(ROOT, name)
    return json.load(open(p)) if os.path.exists(p) else None


def need(*xs):
    return all(x for x in xs)


CL = jload("data_r1/classical.json")
CL2_B3 = jload("data_r1/classical2_b3.json")
CL2_INV = jload("data_r1/classical2_inverse.json")
CL2_TIME = jload("data_r1/classical2_timing.json")
CL2_FLOOR = jload("data_r1/classical2_floor.json")
QUAD = jload("data_r1/check_quad.json")
if QUAD:
    NUM["quad.max_rel"] = BS + "ensuremath{" + tex_sci(max(v["rel_diff"] for v in QUAD.values()), 0) + "}"


def cheb(ic, eps, N, dt):
    for c in CL[f"{ic}_eps{eps}"]["cheb"]:
        if c["N"] == N and abs(c["dt"] - dt) < 1e-12:
            return c
    return None


def cheb_time(eps, N, dt):
    """wall time of the classical solver measured on the server (same machine as the PINN timing runs) if available."""
    if CL2_TIME:
        k = f"eps{eps:g}_N{N}_dt{dt:g}"
        if k in CL2_TIME:
            return CL2_TIME[k]["median"], "server"
    return None, None


def boot_ratio(a, b, n=20000, seed_=0):
    rng = np.random.default_rng(seed_)
    a = np.asarray(a); b = np.asarray(b)
    r = [rng.choice(a, len(a)).mean() / rng.choice(b, len(b)).mean() for _ in range(n)]
    return float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))


def ttime(tag):
    r = run1(tag)
    return r["train_time"] if r else None


# ================================================================================================ classical / ODE numbers
for ic in ("compat", "orig"):
    for eps in ("0.01", "0.1"):
        rec = CL[f"{ic}_eps{eps}"]
        put(f"ode.{ic}.{eps}.l2", f_sci(rec["ode"]["1.0"]["l2"]))
        for c in rec["cheb"]:
            put(f"cheb.{ic}.{eps}.N{c['N']}.dt{c['dt']:g}.l2", f_sci(c["l2"]))


# ================================================================================================ B1
fh32, fh32E, fh64, fh128 = runs("fh_K32_s*"), runs("fh_K32_E50_s*"), runs("fh_K64_s*"), runs("fh_K128_s*")
van, lbvan, lbh32, lbh64 = runs("fwd_compat_vanilla_s*"), runs("lb_vanilla_s*"), runs("lb_hard_K32_s*"), runs("lb_hard_K64_s*")
so32, so32E, so64 = runs("fwd_compat_K32_s*"), runs("fwd_compat_K32_E50_s*"), runs("fwd_compat_K64_s*")
c80 = cheb("compat", "0.01", 80, 5e-4)
ode1 = CL["compat_eps0.01"]["ode"]["1.0"]["l2"]

if fh32:
    put("b1.fh_K32.l2", f_stat(l2(fh32)))
    put("b1.ratio.ode_fh_K32", num(ode1 / np.mean(l2(fh32))))
    put("b1.ratio.fh_K32_cheb80", num(np.mean(l2(fh32)) / c80["l2"]))
if van:
    put("b1.van.l2", f_stat(l2(van)))
if fh32 and van:
    a, b = l2(van), l2(fh32)
    put("b1.ratio.van_fhK32", num(np.mean(a) / np.mean(b)))
    lo, hi = boot_ratio(a, b)
    put("b1.ratio.van_fhK32.ci", f"[{num(lo)}, {num(hi)}]")
    bv, bf = by_seed(van), by_seed(fh32)
    pr = [bv[s]["errors"]["1.0"]["l2"] / bf[s]["errors"]["1.0"]["l2"] for s in bv if s in bf]
    put("b1.ratio.van_fhK32.range", f"{num(min(pr))}--{num(max(pr))}")
    put("b1.ratio.van_fhK32.linf", num(np.mean(l2(van, norm="linf")) / np.mean(l2(fh32, norm="linf"))))
    REPORT.append(f"B1 vanilla/fh_K32 ratio of means {np.mean(a) / np.mean(b):.1f}, CI [{lo:.1f},{hi:.1f}], per seed {[round(x, 1) for x in pr]}")
if lbvan:
    put("b1.lb_van.l2", f_stat(l2(lbvan)))
    put("b1.lb_van.n", str(len(lbvan)))
    if fh32:
        put("b1.ratio.lbvan_fhK32", num(np.mean(l2(lbvan)) / np.mean(l2(fh32))))
if lbh32:
    put("b1.lb_hard_K32.l2", f_stat(l2(lbh32)))
    if lbvan:
        put("b1.ratio.lb", num(np.mean(l2(lbvan)) / np.mean(l2(lbh32))))
if so32:
    put("b1.so_K32.l2", f_stat(l2(so32)))
if so32E:
    put("b1.so_K32E.l2", f_stat(l2(so32E)))
if fh32 and CL2_FLOOR:
    put("b1.floor_ratio.K32", num(np.mean(l2(fh32)) / CL2_FLOOR["b1.K32.neumann"]))
if lbh64:
    put("b1.lb_hard_K64.l2", f_stat(l2(lbh64)))
put("b1.cheb80.l2", f_sci(c80["l2"]))
put("b1.cheb80.dt1e-4.l2", f_sci(cheb("compat", "0.01", 80, 1e-4)["l2"]))
tc80, where = cheb_time(0.01, 80, 5e-4)
if tc80:
    put("b1.cheb80.time", num(tc80))
    if ttime("timing_hard_K32_s42"):
        put("b1.speed.fh_K32_cheb80", num(ttime("timing_hard_K32_s42") / tc80))
if ttime("timing_vanilla_s42") and ttime("timing_hard_K32_s42"):
    put("b1.speed.van_fhK32", num(ttime("timing_vanilla_s42") / ttime("timing_hard_K32_s42")))
for tag in ("timing_vanilla_s42", "timing_soft_K32_s42", "timing_hard_K32_s42", "timing_hard_K64_s42"):
    if ttime(tag):
        put(f"time.{tag}", num(ttime(tag), 3))

# ================================================================================================ B2
b2h, b2v, b2s = runs("b2_hard_K32_s*"), runs("b2_vanilla_s*"), runs("b2_soft_K32_s*")
lbb2v, lbb2h = runs("lb_b2_vanilla_s*"), runs("lb_b2_hard_K32_s*")
ode2 = CL["compat_eps0.1"]["ode"]["1.0"]["l2"]
put("b2.cheb80.l2", f_sci(cheb("compat", "0.1", 80, 5e-4)["l2"]))
if b2h:
    put("b2.hard_K32.l2", f_stat(l2(b2h)))
    put("b2.ratio.ode_hard", num(ode2 / np.mean(l2(b2h))))
if b2v:
    put("b2.van.l2", f_stat(l2(b2v)))
if b2s:
    put("b2.soft_K32.l2", f_stat(l2(b2s)))
    put("b2.soft.maxu", m(f"{max(r['max_abs_u'] for r in b2s):.2f}"))
    if "balance" in b2s[0]:
        put("b2.soft.wall", num(max(r["balance"]["max_abs_ux_wall"] for r in b2s)))
if b2h:
    put("b2.ratio.hard_cheb80", num(np.mean(l2(b2h)) / cheb("compat", "0.1", 80, 5e-4)["l2"]))
if so32 and fh32:
    put("b1.ratio.so_fh_K32", num(np.mean(l2(so32)) / np.mean(l2(fh32))))
if b2h and b2v:
    put("b2.ratio.van_hard", num(np.mean(l2(b2v)) / np.mean(l2(b2h))))
if lbb2v:
    put("b2.lb_van.l2", f_stat(l2(lbb2v)))
if lbb2h:
    put("b2.lb_hard_K32.l2", f_stat(l2(lbb2h)))
    if lbb2v:
        put("b2.ratio.lb_hard_van", num(np.mean(l2(lbb2h)) / np.mean(l2(lbb2v))))

# ================================================================================================ B3 (errors at t = 5)
b3K32, b3K32E = runs("b3_K32_E0_s*"), runs("b3_K32_E50_s*")
b3K64, b3K64E, b3v = runs("b3_K64_E0_s*"), runs("b3_K64_E50_s*"), runs("b3_vanilla_s*")
if b3K32:
    put("b3.K32_E0.l2", f_stat(l2(b3K32, "5.0")))
if b3K64:
    put("b3.K64_E0.l2", f_stat(l2(b3K64, "5.0")))
if b3K64E:
    put("b3.K64_E50.l2", f_stat(l2(b3K64E, "5.0")))
if b3v:
    put("b3.van.l2", f_stat(l2(b3v, "5.0")))
if CL2_B3:
    def b3cheb(N, dt):
        return [c for c in CL2_B3["cheb"] if c["N"] == N and abs(c["dt"] - dt) < 1e-12][0]
    put("b3.cheb.l2", f_sci(b3cheb(80, 5e-4)["errors"]["5.0"]["l2"]))
    put("b3.cheb128.l2", f_sci(b3cheb(128, 1e-4)["errors"]["5.0"]["l2"]))
    put("b3.ode.l2", f_sci(CL2_B3["ode"]["5.0"]["l2"]))
if CL2_FLOOR:
    for k, v in CL2_FLOOR.items():
        put(f"floor.{k}", f_sci(v, 1))

# ================================================================================================ mode number K
vs64, vs32, vs128 = runs("vs_K64_bcfix_s*"), runs("vs_K32_bcfix_s*"), runs("vs_K128_bcfix_s*")
so128 = runs("var_K128_uniform_soft")
if so64:
    put("kk.soft64.wall", num(max(r["balance"]["max_abs_ux_wall"] for r in so64)))
    put("kk.soft64.l2", f_stat(l2(so64)))
if fh64:
    put("kk.hard64.l2", f_stat(l2(fh64)))
if fh128:
    put("kk.hard128.l2", f_stat(l2(fh128)))
if so128:
    put("kk.soft128.l2", f_stat(l2(so128)))
if vs64:
    put("kk.fix64.l2", f_stat(l2(vs64)))
if vs128:
    put("kk.fix128.l2", f_stat(l2(vs128)))
v64cs, s64_42 = run1("var_K64_cgl_soft"), by_seed(so64).get(42)
if v64cs and s64_42:
    put("kk.cgl_soft_gain", num(s64_42["errors"]["1.0"]["l2"] / v64cs["errors"]["1.0"]["l2"]))
kh = {k: run1(k) for k in ("kh_K64_x2", "kh_K128_x2", "kh_K256_x1")}
for k, r in kh.items():
    if r:
        put(f"kk.{k}.l2", f_sci(r["errors"]["1.0"]["l2"]))
if fh128 and 42 in by_seed(fh128):
    put("kk.hard128_s42", f_sci(by_seed(fh128)[42]["errors"]["1.0"]["l2"]))
_rep = []
for t_tag, o_tag in (("timing_hard_K32_s42", "fh_K32_s42"), ("timing_vanilla_s42", "fwd_compat_vanilla_s42"), ("timing_hard_K64_s42", "fh_K64_s42"),
                     ("timing_soft_K32_s42", "fwd_compat_K32_s42")):
    a, b = run1(t_tag), run1(o_tag)
    if a and b:
        _rep.append((o_tag, abs(a["errors"]["1.0"]["l2"] / b["errors"]["1.0"]["l2"] - 1)))
if _rep:
    ws = set(open(os.path.join(ROOT, "data_r1", "runs_workstation.txt")).read().split())
    cross = [d for tag, d in _rep if tag in ws]
    same = [d for tag, d in _rep if tag not in ws]
    if cross:
        put("proto.repro", pct(100 * max(cross), 1))
    if same:
        put("proto.repro_same", pct(100 * max(same), 1) if max(same) > 1e-12 else "0")
    REPORT.append(f"repeat runs (timing group vs original): {[(t, round(100 * d, 3)) for t, d in _rep]} %")
if so128 and "balance" in so128[0]:
    put("kk.soft128.wall", num(max(r["balance"]["max_abs_ux_wall"] for r in so128)))
    put("kk.soft128.maxu", m(f"{max(r['max_abs_u'] for r in so128):.2f}"))
if fh32:
    put("kk.hard32.l2", f_stat(l2(fh32)))
for tag, key in (("var_K64_uniform_hard", "kk.hard_uniform64"), ("var_K64_cgl_soft", "kk.soft_cgl64")):
    r = run1(tag)
    if r:
        put(key, f_sci(r["errors"]["1.0"]["l2"]))
if s64_42:
    put("kk.soft_uniform64_s42", f_sci(s64_42["errors"]["1.0"]["l2"]))
if fh64 and 42 in by_seed(fh64):
    put("kk.hard_cgl64_s42", f_sci(by_seed(fh64)[42]["errors"]["1.0"]["l2"]))
if vs64:
    put("kk.fix64.rms", f_sci(np.mean([r["balance"]["rms_residual_L2"] for r in vs64]), 1))
    put("kk.fix64.wall", num(max(r["balance"]["max_abs_ux_wall"] for r in vs64)))
    put("kk.fix64.maxu", m(f"{max(r['max_abs_u'] for r in vs64):.2f}"))
if fh64:
    put("kk.hard64.rms", f_sci(np.mean([r["balance"]["rms_residual_L2"] for r in fh64]), 1))
if fh32 and fh64:
    put("kk.ratio64_32", num(np.mean(l2(fh64)) / np.mean(l2(fh32))))
if fh64 and fh128:
    put("kk.ratio128_64", num(np.mean(l2(fh128)) / np.mean(l2(fh64))))
if so64 and fh64:
    put("kk.ratio_soft_hard64", num(np.mean(l2(so64)) / np.mean(l2(fh64))))

# ================================================================================================ energy
B1_SPEC = ["fh_K32_s*", "fh_K32_E50_s*", "ah_lamE*_s*", "fh_K64_s*", "fh_K128_s*", "fwd_compat_K32_s*", "fwd_compat_K32_E50_s*",
           "fwd_compat_K64_s*", "vs_K*_bcfix_s*", "lb_hard_K*_s*", "var_K*", "kh_K*", "reg_K64_E50_s*"]
b1spec = []
for p in B1_SPEC:
    b1spec += runs(p)
if b1spec:
    # every B1 Spectral-PINN with the compatible IC: B1_SPEC + poor-trajectory regimes + spectral timing repeats
    b1en = b1spec + runs("rh_*_s*") + [r for r in runs("timing_*_s*") if r["job"]["method"] == "spectral"]
    inc = {k: sum(1 for r in b1en if r["energy"][k]["n_increases"] > 0) for k in ("dt0.02", "dt0.002", "dt0.0002")}
    put("en.n_runs", str(len(b1en)))
    put("en.n_runs_viol", str(max(inc.values())))
    put("en.min_step_all", f_sci(min(-r["energy"]["dt0.02"]["max_step"] for r in b1spec), 1))
    put("en.max_l2", num(max(r["errors"]["1.0"]["l2"] for r in b1spec), 2))
    if fh32 + fh32E:
        put("en.min_step", f_sci(min(-r["energy"]["dt0.02"]["max_step"] for r in fh32 + fh32E), 1))
    REPORT.append(f"energy: {len(b1spec)} B1 spectral runs (default regime); runs with increases per resolution {inc}; "
                  f"offenders: {[r['job']['tag'] for r in b1spec if any(r['energy'][k]['n_increases'] for k in r['energy'])]}")
if fh32 and fh32E:
    a, b = by_seed(fh32), by_seed(fh32E)
    md = max(float(np.max(np.abs(np.array(a[s]["E_traj_dt0.02"]) - np.array(b[s]["E_traj_dt0.02"])))) for s in a if s in b)
    put("en.max_diff", f_sci(md, 1))
    nz = sum(r["penalty"]["n_epochs_nonzero"] for r in fh32E)
    tot = sum(r["epochs"] for r in fh32E)
    put("en.nonzero", f"{nz} of the {num(tot, 5)}")
_allspec = [json.load(open(f)) for f in glob.glob(os.path.join(RUNS, "*.json"))]
_allspec = [r for r in _allspec if r["job"]["method"] == "spectral" and "energy" in r and not r["job"]["tag"].startswith("smoke")]
if _allspec:
    put("en.n_all_spec", str(len(_allspec)))
    put("en.n_all_spec_viol", str(sum(1 for r in _allspec if any(v["n_increases"] for v in r["energy"].values()))))
if b3v:
    put("en.b3van_inc", str(b3v[0]["energy"]["dt0.02"]["n_increases"]))
    put("en.b3van_int", str(b3v[0]["energy"]["dt0.02"]["n_points"] - 1))
if van:
    put("en.van_viol", str(sum(r["energy"]["dt0.02"]["n_increases"] for r in van + lbvan + runs("timing_vanilla_s*"))))
    put("en.van_n", str(len(van + lbvan + runs("timing_vanilla_s*"))))
put("reg.short_epochs", "300")
_lam = [fh32, runs("ah_lamE1_s*"), runs("ah_lamE10_s*"), fh32E, runs("ah_lamE500_s*"), runs("ah_lamE5000_s*")]
_lam = [rr for rr in _lam if rr]
if _lam:
    means = [np.mean(l2(rr)) for rr in _lam]
    put("en.lam.min", f_sci(min(means))); put("en.lam.max", f_sci(max(means)))
    put("en.lam.sd", f_sci(stat(l2(fh32))["sd"], 1) if fh32 else "")
    pen = [r for rr in _lam for r in rr if r["penalty"]["lam_E"] > 0]
    put("en.lam.max_nonzero", str(max(r["penalty"]["n_epochs_nonzero"] for r in pen)))
    put("en.lam.maxpen", f_sci(max(r["penalty"]["max_lamE_LE"] for r in pen), 1))
    put("en.lam.nruns", str(sum(len(rr) for rr in _lam)))
for base in ("short", "nowarm", "pde0"):
    for lam in (0, 50):
        rr = runs(f"rh_{base}_E{lam}_s*")
        if rr:
            put(f"reg.{base}.E{lam}", f_stat(l2(rr)))
            put(f"reg.{base}.E{lam}.maxnz", str(max(r["penalty"]["n_epochs_nonzero"] for r in rr)))
            put(f"reg.{base}.E{lam}.maxpen", f_sci(max(r["penalty"]["max_lamE_LE"] for r in rr), 1) if max(r["penalty"]["max_lamE_LE"] for r in rr) > 0 else "0")
_nw0, _nw50 = by_seed(runs("rh_nowarm_E0_s*")), by_seed(runs("rh_nowarm_E50_s*"))
if _nw0 and _nw50:
    _fac = sorted(_nw50[s]["errors"]["1.0"]["l2"] / _nw0[s]["errors"]["1.0"]["l2"] for s in _nw0 if s in _nw50)
    put("reg.nowarm.factors", ", ".join(num(x) for x in _fac[:-1]) + " and " + num(_fac[-1]))
_r64 = runs("reg_K64_E50_s*")
if _r64 and so64:
    sd_ = set(seed(r) for r in _r64)
    put("reg.soft64.E50", f_stat(l2(_r64)))
    put("reg.soft64.E0", f_stat(l2([r for r in so64 if seed(r) in sd_])))
    nz = [r["penalty"]["n_epochs_nonzero"] for r in _r64]
    put("reg.soft64.nz", f"{min(nz)}--{max(nz)}")
_b3E = runs("b3_K64_E50_s*") + runs("b3_K32_E50_s*")
if _b3E:
    nz = [r["penalty"]["n_epochs_nonzero"] for r in _b3E]
    put("b3.E50.nz", f"{min(nz)}--{max(nz)}")
for k, rr in (("K32", b3K32), ("K64", b3K64)):
    if rr:
        put(f"b3.{k}.t02", f_stat(l2(rr, "0.2")))
        put(f"b3.{k}.t1", f_stat(l2(rr, "1.0")))
_reg_all = [r for p in ("rh_*_s*", "reg_K64_E50_s*") for r in runs(p)] + [r for r in so64 if seed(r) in {seed(q) for q in runs("reg_K64_E50_s*")}]
if _reg_all:
    put("reg.n_runs", str(len(_reg_all)))
    put("reg.n_viol", str(sum(r["energy"][k]["n_increases"] for r in _reg_all for k in r["energy"])))

# ================================================================================================ maximum bound
b1van_all = van + lbvan + runs("timing_vanilla_s*")
if b1spec:
    put("mbp.max_b1", m(f"{max(r['max_abs_u'] for r in b1spec):.4f}"))
    # "revised" = exact Neumann basis AND Gauss-Lobatto grid; regimes reported separately
    hard = [r for r in b1spec + runs("timing_hard_K*_s*") if r["job"].get("hard_bc") and r["job"].get("grid") == "cgl"]
    if hard:
        put("mbp.max_hard", m(f"{max(r['max_abs_u'] for r in hard):.4f}"))
        put("mbp.excess", f_sci(max(r["max_abs_u"] for r in hard) - 1, 1))
    hreg = [r for r in runs("rh_*_s*") if r["job"].get("hard_bc") and r["job"].get("grid") == "cgl"]
    if hreg:
        put("mbp.n_hard_reg", str(len(hreg)))
        put("mbp.max_hard_reg", m(f"{max(r['max_abs_u'] for r in hreg):.4f}"))
        put("mbp.excess_reg", f_sci(max(r["max_abs_u"] for r in hreg) - 1, 1))
if b1van_all:
    put("mbp.max_van", m(f"{max(r['max_abs_u'] for r in b1van_all):.4f}"))
if b1spec and b1van_all:
    put("mbp.excess_all", f_sci(max(r["max_abs_u"] for r in b1spec + b1van_all) - 1, 1))
    allb1 = b1spec + b1van_all
    put("mbp.n_b1", str(len(allb1)))
    put("mbp.n_over", str(sum(1 for r in allb1 if r["max_abs_u"] > 1.0)))
    hard = [r for r in b1spec + runs("timing_hard_K*_s*") if r["job"].get("hard_bc") and r["job"].get("grid") == "cgl"]
    put("mbp.n_hard", str(len(hard)))
    put("mbp.n_hard_over", str(sum(1 for r in hard if r["max_abs_u"] > 1.0)))
    put("mbp.n_van", str(len(b1van_all)))
    put("mbp.n_van_over", str(sum(1 for r in b1van_all if r["max_abs_u"] > 1.0)))
    put("mbp.van_excess", f_sci(max(r["max_abs_u"] for r in b1van_all) - 1, 1))
    REPORT.append("MBP: B1 runs with max|u|>1: %d of %d; revised %d of %d; vanilla %d of %d" % (
        sum(1 for r in allb1 if r["max_abs_u"] > 1), len(allb1), sum(1 for r in hard if r["max_abs_u"] > 1), len(hard),
        sum(1 for r in b1van_all if r["max_abs_u"] > 1), len(b1van_all)))
_zr = os.path.join(ROOT, "data_r1", "ref", "compat_eps0.01_T1_n4.npz")
if os.path.exists(_zr):
    _z = np.load(_zr)
    put("mbp.ref_max_t1", m(f"{np.abs(_z['t3']).max():.4f}"))

# ================================================================================================ energy balance
allbal = [r for r in glob.glob(os.path.join(RUNS, "*.json"))]
allbal = [json.load(open(f)) for f in allbal]
allbal = [r for r in allbal if "balance" in r and not r["job"]["tag"].startswith("smoke")]
if allbal:
    put("bal.id_err", f_sci(max(r["balance"]["identity_rel_err"] for r in allbal), 0))
    put("bal.n", str(len(allbal)))
if fh32:
    ir = np.mean([r["balance"]["int_half_r2"] for r in fh32])
    dE = np.mean([r["balance"]["E_first"] - r["balance"]["E_last"] for r in fh32])
    put("bal.int_r_hard", f_sci(ir, 1))
    put("bal.dE_hard", f_sci(dE, 1))
    put("bal.frac_hard", f_sci(ir / dE, 0))
    put("bal.int_r_hard_max", f_sci(max(r["balance"]["int_half_r2"] for r in fh32), 1))
    put("bal.frac_hard_max", f_sci(max(r["balance"]["int_half_r2"] / (r["balance"]["E_first"] - r["balance"]["E_last"]) for r in fh32), 0))
for tag, rr in (("soft32", so32), ("soft64", so64)):
    if rr:
        put(f"bal.int_b_{tag}", f_sci(np.mean([r["balance"]["int_abs_bterm"] for r in rr]), 1))
if so64:
    put("bal.int_b_soft", NUM["bal.int_b_soft64"])
    put("bal.int_r_soft64", f_sci(np.mean([r["balance"]["int_half_r2"] for r in so64]), 1))

# ================================================================================================ inverse
sweep = []
for var in ("soft", "hard"):
    for e in (0.05, 0.1):
        for i in (0.5, 1.0, 1.5, 2.0, 4.0):
            r = run1(f"inv_adam_{var}_e{e}_i{i}")
            if r:
                sweep.append(r)
if sweep:
    ri = [r["ratio_to_init"] for r in sweep]
    put("inv.ratio_mean", num(np.mean(ri)))
    put("inv.ratio_range", f"{num(min(ri))}--{num(max(ri))}")
    put("inv.n_sweep", str(len(sweep)))
    steps = []
    for r in sweep:
        tr = np.log(np.array(r["eps_traj_every200"]))
        fin = np.log(r["eps_hat"])
        idx = [i for i in range(len(tr)) if np.all(np.abs(tr[i:] - fin) < 0.01)]
        if idx:
            steps.append(200 * idx[0])
    if steps:
        put("inv.steps", num(np.median(steps), 2))
        REPORT.append(f"inverse joint Adam: epochs until log eps within 0.01 of its final value: {sorted(steps)}")
    # share of the total displacement of log(eps) reached after the first 200 epochs
    fr = []
    for r in sweep:
        tr = np.log(np.array(r["eps_traj_every200"]))
        tot = np.log(r["eps_hat"]) - np.log(r["eps_true"] * r["job"].get("eps_init_f", 2.0))
        if abs(tot) > 1e-3:
            fr.append((tr[1] - np.log(r["eps_true"] * r["job"].get("eps_init_f", 2.0))) / tot)
    if fr:
        put("inv.frac200", pct(100 * float(np.median(fr))))
        put("inv.frac200_min", pct(100 * float(np.min(fr))))
for tag, key in (("inv_adam_hard_e0.05_i2.0", "base"), ("inv_ab_mult1", "mult1"), ("inv_ab_mult100", "mult100"), ("inv_ab_ld2000", "ld2000"),
                 ("inv_ab_ld2000_i1.0", "ld2000_i1"), ("inv_ab_ep16000", "ep16000"), ("inv_ab_nx10", "nx10"), ("inv_ab_nx160", "nx160"),
                 ("inv_ab_nt20", "nt20"), ("inv_ab_K16", "K16"), ("inv_ab_K64", "K64"), ("inv_ab_lpde5", "lpde5"), ("inv_ab_lpde500", "lpde500")):
    r = run1(tag)
    if r:
        put(f"inv.ab.{key}.rt", num(r["ratio_to_true"], 3))
        put(f"inv.ab.{key}.ri", num(r["ratio_to_init"], 3))
        put(f"inv.ab.{key}.disp", num(abs(np.log(r["ratio_to_init"])), 2))
_r = run1("inv_adam_hard_e0.05_i1.0")
if _r:
    put("inv.i1.rt", num(_r["ratio_to_true"], 2))
_abl = [run1(t) for t in ("inv_adam_hard_e0.05_i2.0", "inv_ab_nx10", "inv_ab_nx160", "inv_ab_nt20", "inv_ab_K16", "inv_ab_K64", "inv_ab_lpde5", "inv_ab_lpde500")]
_abl = [r for r in _abl if r]
if _abl:
    put("inv.ab.design_range", f"{num(min(r['ratio_to_true'] for r in _abl), 3)}--{num(max(r['ratio_to_true'] for r in _abl), 3)}")
for e in (0.05, 0.1):
    r = run1(f"inv_vp_e{e}_ld2000")
    if r:
        put("inv.vp_err_" + {0.05: "005", 0.1: "01"}[e], pct(r["rel_err_pct"]))
    for ld in (20, 200, 20000):
        r = run1(f"inv_vp_e{e}_ld{ld}")
        if r:
            put(f"inv.vp.e{e:g}.ld{ld}", pct(r["rel_err_pct"]))
    for nz in (0.01, 0.05):
        rr = [run1(f"inv_vpnoise_e{e}_n{nz}_r{s}") for s in range(1, 6)]
        rr = [r for r in rr if r]
        if rr:
            v = [r["rel_err_pct"] for r in rr]
            put(f"inv.vpn.e{e:g}.n{nz:g}", pct(np.mean(v)))
            put(f"inv.vpn.e{e:g}.n{nz:g}.bias", pct(100 * (np.mean([r["ratio_to_true"] for r in rr]) - 1)))
            put(f"inv.vpn.e{e:g}.n{nz:g}.maxratio", num(max(r["ratio_to_true"] for r in rr), 2))
        if CL2_INV:
            c = CL2_INV[str(e)]["runs"][str(nz)]
            put(f"inv.cl.e{e:g}.n{nz:g}", pct(c["rel_err_pct_mean"]))
            put(f"inv.cl.e{e:g}.n{nz:g}.rmse", pct(c["rel_rmse_pct"]))
            put(f"inv.fisher.e{e:g}.n{nz:g}", pct(c["fisher_rel_se_pct"]))
            p5 = [100 * abs(c["eps_hat"][str(s)] - e) / e for s in range(1, 6)]
            put(f"inv.cl5.e{e:g}.n{nz:g}", pct(np.mean(p5)))
for tag in ("nx10", "nx160", "nt3", "nt20", "lpde5", "lpde500", "nophaseA", "K64"):
    r = run1(f"inv_vp_e0.05_{tag}")
    if r:
        put(f"inv.vp.{tag}", pct(r["rel_err_pct"]))
        put(f"inv.vp.{tag}.rt", num(r["ratio_to_true"], 3))
for nz in (0.0, 0.01):
    r = run1(f"inv_vp_e0.01_n{nz}")
    if r:
        put(f"inv.vp001.n{nz:g}", pct(r["rel_err_pct"]))
        put(f"inv.vp001.n{nz:g}.rt", num(r["ratio_to_true"], 3))
if CL2_INV:
    put("ident.S001", f_sci(CL2_INV["0.01"]["S_norm_dlogeps"], 1))
    put("ident.S005", f_sci(CL2_INV["0.05"]["S_norm_dlogeps"], 1))
    put("ident.S01", f_sci(CL2_INV["0.1"]["S_norm_dlogeps"], 1))
    put("ident.se001", pct(CL2_INV["0.01"]["runs"]["0.01"]["fisher_rel_se_pct"]))
    put("ident.se005", pct(CL2_INV["0.05"]["runs"]["0.01"]["fisher_rel_se_pct"]))
    put("ident.cl001", pct(CL2_INV["0.01"]["runs"]["0.01"]["rel_rmse_pct"]))
    put("ident.cl005", pct(CL2_INV["0.05"]["runs"]["0.01"]["rel_rmse_pct"]))
    put("ident.n_rep", str(CL2_INV["0.05"]["runs"]["0.01"]["n"]))

# ================================================================================================ original initial condition (Appendix)
_ro32, _ro64 = runs("fwd_orig_K32_s*"), runs("fwd_orig_K64_s*")
for key, ro, rc in (("K32", _ro32, so32), ("K64", _ro64, so64)):
    if ro:
        sd = set(seed(r) for r in ro)
        rcs = [r for r in rc if seed(r) in sd]
        put(f"orig.{key}", f_stat(l2(ro)))
        if rcs:
            put(f"orig.{key}c", f_stat(l2(rcs)))
            if "balance" in rcs[0]:
                put(f"orig.{key}c.wall", num(max(r["balance"]["max_abs_ux_wall"] for r in rcs)))
        if "balance" in ro[0]:
            put(f"orig.{key}.wall", num(max(r["balance"]["max_abs_ux_wall"] for r in ro)))
if _ro32 or _ro64:
    put("orig.viol", str(sum(r["energy"]["dt0.02"]["n_increases"] for r in _ro32 + _ro64)))

# ================================================================================================ 2D (unchanged data of the first version)
ex2 = jload("data/experiment_2d.json")
ms2 = jload("data/symmetry_2d_multiseed.json")
st2 = jload("data/stress_2d_sweep.json")
fv = [ex2["violations"][t]["fourier_spectral"] for t in ex2["violations"]]
put("s2d.fmin", f_sci(min(fv), 1)); put("s2d.fmax", f_sci(max(fv), 1))
S0 = np.array([r["violations"]["0.0"] for r in ms2["results"]]); S1 = np.array([r["violations"]["1.0"] for r in ms2["results"]])
put("s2d.growth", num(S1.mean() / S0.mean()))
put("s2d.growth_range", f"{num((S1 / S0).min())}--{num((S1 / S0).max())}")
mono = all(all(np.diff([r["violations"][t] for t in ("0.0", "0.2", "0.5", "1.0")]) > 0) for r in ms2["results"])
REPORT.append(f"2D vanilla S(t) monotonically increasing in every seed: {mono}")
s1 = [r["violations"]["1.0"] for r in st2["results"]]
put("s2d.smin", f_sci(min(s1), 2)); put("s2d.smax", f_sci(max(s1), 2))
put("s2d.rms_max", f_sci(np.sqrt(max(s1)), 0))


# ================================================================================================ tables
def wide(caption, label, colspec, header, rows, note=None):
    out = [BS + "begin{table}[H]", BS + "begin{adjustwidth}{-" + BS + "extralength}{0cm}", BS + "tablesize{" + BS + "footnotesize}",
           BS + "caption{" + caption + BS + "label{" + label + "}}", BS + "newcolumntype{C}{>{" + BS + "centering" + BS + "arraybackslash}X}",
           BS + "begin{tabularx}{" + BS + "fulllength}{" + colspec + "}", BS + "toprule", header + BS + BS, BS + "midrule"]
    out += rows
    out += [BS + "bottomrule", BS + "end{tabularx}", BS + "end{adjustwidth}"]
    if note:
        out += [BS + "noindent{" + BS + "footnotesize{}" + note + "}"]
    out += [BS + "end{table}", ""]
    return "\n".join(out)


def narrow(caption, label, colspec, header, rows):
    out = [BS + "begin{table}[H]", BS + "tablesize{" + BS + "footnotesize}", BS + "caption{" + caption + BS + "label{" + label + "}}",
           BS + "newcolumntype{C}{>{" + BS + "centering" + BS + "arraybackslash}X}", BS + "begin{tabularx}{" + BS + "textwidth}{" + colspec + "}",
           BS + "toprule", header + BS + BS, BS + "midrule"]
    out += rows
    out += [BS + "bottomrule", BS + "end{tabularx}", BS + "end{table}", ""]
    return "\n".join(out)


def row(*cells):
    return " & ".join(cells) + BS + BS


def group_row(text, ncol):
    return BS + "multicolumn{" + str(ncol) + "}{l}{" + BS + "emph{" + text + "}}" + BS + BS


DASH = "---"


def incr_cell(rr, key="dt0.02"):
    tot = sum(r["energy"][key]["n_increases"] for r in rr)
    k = sum(1 for r in rr if r["energy"][key]["n_increases"] > 0)
    return "0" if tot == 0 else f"{tot} ({k}/{len(rr)})"


def pinn_row(label, rr, t="1.0", time_tag=None):
    if not rr:
        return row(label, DASH, "(not run)", "", "", "", "")
    tt = num(ttime(time_tag), 3) if time_tag and ttime(time_tag) else DASH
    return row(label, str(len(rr)), f_stat(l2(rr, t)), f_stat(l2(rr, t, "linf")), m(f"{max(r['max_abs_u'] for r in rr):.4f}"), incr_cell(rr), tt)


def cl_row(label, l2v, linf, tsec=None):
    return row(label, DASH, f_sci(l2v), f_sci(linf), DASH, DASH, num(tsec, 2) if tsec else DASH)


def tab_forward():
    H = row(BS + "textbf{Method}", m(BS + "boldsymbol{n}"), m(BS + "boldsymbol{" + BS + "|e" + BS + "|_{L^2}}"), m(BS + "boldsymbol{" + BS + "|e" + BS + "|_{L^" + BS + "infty}}"),
            m(BS + "boldsymbol{" + BS + "max|u|}"), BS + "textbf{E incr.}", BS + "textbf{Time (s)}")[:-2]
    R = []
    R.append(group_row(r"B1: $\varepsilon=0.01$, $T=1$ (errors at $t=1$)", 7))
    for (N, dt) in ((80, 5e-4), (80, 1e-4)):
        c = cheb("compat", "0.01", N, dt)
        R.append(cl_row(f"Chebyshev, $N={N}$, $" + BS + f"delta t={"10^{-4}" if dt == 1e-4 else '5' + BS + 'times10^{-4}'}$", c["l2"], c["linf"], cheb_time(0.01, N, dt)[0]))
    o = CL["compat_eps0.01"]["ode"]["1.0"]
    R.append(cl_row(r"Reaction-only ODE ($\varepsilon=0$)", o["l2"], o["linf"]))
    R.append(pinn_row("Vanilla PINN", van, time_tag="timing_vanilla_s42"))
    R.append(pinn_row("Vanilla PINN + L-BFGS", lbvan))
    R.append(pinn_row(r"Spectral-PINN as submitted, $K=32$", so32, time_tag="timing_soft_K32_s42"))
    R.append(pinn_row(r"\quad + energy penalty ($\lambda_E=50$)", so32E))
    R.append(pinn_row(r"Spectral-PINN as submitted, $K=64$", so64))
    R.append(pinn_row(r"Spectral-PINN revised, $K=32$", fh32, time_tag="timing_hard_K32_s42"))
    R.append(pinn_row(r"\quad + energy penalty ($\lambda_E=50$)", fh32E))
    R.append(pinn_row(r"\quad + L-BFGS", lbh32))
    R.append(pinn_row(r"Spectral-PINN revised, $K=64$", fh64, time_tag="timing_hard_K64_s42"))
    R.append(pinn_row(r"\quad + L-BFGS", lbh64))
    R.append(BS + "midrule")
    R.append(group_row(r"B2: $\varepsilon=0.1$, $T=1$ (errors at $t=1$)", 7))
    c = cheb("compat", "0.1", 80, 5e-4)
    R.append(cl_row(r"Chebyshev, $N=80$, $\delta t=5\times10^{-4}$", c["l2"], c["linf"], cheb_time(0.1, 80, 5e-4)[0]))
    o = CL["compat_eps0.1"]["ode"]["1.0"]
    R.append(cl_row(r"Reaction-only ODE ($\varepsilon=0$)", o["l2"], o["linf"]))
    R.append(pinn_row("Vanilla PINN", b2v))
    R.append(pinn_row("Vanilla PINN + L-BFGS", lbb2v))
    R.append(pinn_row(r"Spectral-PINN as submitted, $K=32$", b2s))
    R.append(pinn_row(r"Spectral-PINN revised, $K=32$", b2h))
    R.append(pinn_row(r"\quad + L-BFGS", lbb2h))
    R.append(BS + "midrule")
    R.append(group_row(r"B3: $\varepsilon=0.05$, $T=5$ (errors at $t=5$)", 7))
    if CL2_B3:
        for (N, dt, lab) in ((80, 5e-4, r"5\times10^{-4}"), (128, 1e-4, r"10^{-4}")):
            c = [x for x in CL2_B3["cheb"] if x["N"] == N and abs(x["dt"] - dt) < 1e-12][0]
            R.append(cl_row(f"Chebyshev, $N={N}$, $" + BS + "delta t=" + lab + "$", c["errors"]["5.0"]["l2"], c["errors"]["5.0"]["linf"]))
        o = CL2_B3["ode"]["5.0"]
        R.append(cl_row(r"Reaction-only ODE ($\varepsilon=0$)", o["l2"], o["linf"]))
    R.append(pinn_row("Vanilla PINN", b3v, "5.0"))
    R.append(pinn_row(r"Spectral-PINN revised, $K=32$", b3K32, "5.0"))
    R.append(pinn_row(r"\quad + energy penalty ($\lambda_E=50$)", b3K32E, "5.0"))
    R.append(pinn_row(r"Spectral-PINN revised, $K=64$", b3K64, "5.0"))
    R.append(pinn_row(r"\quad + energy penalty ($\lambda_E=50$)", b3K64E, "5.0"))
    cap = (r"Forward problem: errors against the certified reference on the 200-point grid, mean $\pm$ SD over $n$ seeds. "
           r"``As submitted'': soft Neumann penalty and 300 equispaced collocation points; ``revised'': exact Neumann basis~\eqref{eq:neumann-basis} and Gauss--Lobatto collocation. "
           r"$\max|u|$: largest value over all runs of a row on a $1001\times(100T+1)$ space--time grid. "
           r"E incr.: total number of energy increases on the grid of spacing $0.02$ (in parentheses: number of runs with at least one). "
           r"Time: training time of one run (Adam phase, one CPU thread, separate timing runs on the same machine for the classical solver and the PINNs; ``---'': not measured separately).")
    return wide(cap, "tab:forward", "p{5.2cm}cCCccc", H, R)


def tab_K():
    H = row(m(BS + "boldsymbol{K}"), BS + "textbf{As submitted}", m(BS + "boldsymbol{" + BS + "lambda_{" + BS + "rm BC}=K^4}"), BS + "textbf{Soft, CGL}",
            BS + "textbf{Hard, equisp.}", BS + "textbf{Revised}", m(BS + "boldsymbol{" + BS + "max|u_x(" + BS + "pm1)|}"))[:-2]

    def cell(rr):
        if not rr:
            return DASH
        return m(tex_sci(np.mean(l2(rr)), 2)) + f" ({len(rr)})"

    R = []
    for K in (32, 64, 128, 256):
        soft = runs(f"fwd_compat_K{K}_s*") or runs(f"var_K{K}_uniform_soft")
        fix = runs(f"vs_K{K}_bcfix_s*")
        cs = runs(f"var_K{K}_cgl_soft")
        uh = runs(f"var_K{K}_uniform_hard")
        hd = runs(f"fh_K{K}_s*") or (runs("kh_K256_x1") if K == 256 else [])
        wall = num(max(r["balance"]["max_abs_ux_wall"] for r in soft)) if soft and "balance" in soft[0] else DASH
        R.append(row(str(K), cell(soft), cell(fix), cell(cs), cell(uh), cell(hd), wall))
    R.append(BS + "midrule")
    R.append(group_row(r"Revised variant, seed 42, larger training budget", 7))
    for tag, K, fac in (("kh_K64_x2", 64, 2), ("kh_K128_x2", 128, 2)):
        r = run1(tag)
        base = by_seed(runs(f"fh_K{K}_s*")).get(42)
        R.append(row(str(K), r"\multicolumn{4}{l}{$" + str(fac) + r"\times$ epochs (12\,000 + 2000 warm-up)}",
                     (f_sci(r["errors"]["1.0"]["l2"]) if r else "(not run)") + (f" (default: {f_sci(base['errors']['1.0']['l2'])})" if base else ""), DASH))
    cap = (r"$L^2$ error at $t=1$ on B1 versus the mode number $K$ for the four combinations of boundary treatment (soft penalty / exact Neumann basis) and collocation grid "
           r"(equispaced / Gauss--Lobatto, CGL); mean over the number of seeds in parentheses. ``As submitted'' = soft penalty with effective weight $K^{-4}$, equispaced grid; "
           r"``$\lambda_{\rm BC}=K^4$'' = the same with the penalty rescaled to unit effective weight; ``Revised'' = exact Neumann basis, CGL grid. "
           r"The last column is the largest wall slope of the as-submitted networks (zero for the exact basis by construction).")
    return wide(cap, "tab:K", "cCCCCCc", H, R)


def tab_energy():
    H = row(BS + "textbf{Configuration}", m(BS + "boldsymbol{n}"), m(BS + "boldsymbol{" + BS + "|e" + BS + "|_{L^2}}"),
            BS + "textbf{E incr. } " + m(BS + "boldsymbol{0.02/0.002/0.0002}"), BS + "textbf{Epochs with } " + m(BS + "boldsymbol{" + BS + "mathcal L_E>0}"),
            m(BS + "boldsymbol{" + BS + "max" + BS + "lambda_E" + BS + "mathcal L_E}"))[:-2]
    R = [group_row(r"B1, revised $K=32$: penalty weight", 6)]
    pats = [(0.0, "fh_K32_s*"), (1.0, "ah_lamE1_s*"), (10.0, "ah_lamE10_s*"), (50.0, "fh_K32_E50_s*"), (500.0, "ah_lamE500_s*"), (5000.0, "ah_lamE5000_s*")]

    def erow(label, rr):
        if not rr:
            return row(label, DASH, "(not run)", "", "", "")
        inc = "/".join(incr_cell(rr, k) for k in ("dt0.02", "dt0.002", "dt0.0002") if k in rr[0]["energy"])
        nz = sum(r["penalty"]["n_epochs_nonzero"] for r in rr)
        mx = max(r["penalty"]["max_lamE_LE"] for r in rr)
        if all(r["penalty"]["lam_E"] == 0 for r in rr):
            return row(label, str(len(rr)), f_stat(l2(rr)), inc, DASH, DASH)
        return row(label, str(len(rr)), f_stat(l2(rr)), inc, str(nz), f_sci(mx, 1) if mx > 0 else "0")

    for lam, p in pats:
        R.append(erow(f"$" + BS + f"lambda_E={lam:g}$", runs(p)))
    R.append(BS + "midrule")
    R.append(group_row(r"Regimes with a poor trajectory ($K=32$ revised unless stated; 3 seeds)", 6))
    for lab, base in ((r"300 epochs (100 warm-up)", "rh_short"), (r"no warm-up", "rh_nowarm"), (r"1500 epochs, lr $10^{-4}$", "rh_pde0")):
        for lam in (0, 50):
            R.append(erow(lab + f", $" + BS + f"lambda_E={lam}$", runs(f"{base}_E{lam}_s*")))
    R.append(erow(r"$K=64$ as submitted, $\lambda_E=0$", [r for r in so64 if seed(r) in {seed(q) for q in runs("reg_K64_E50_s*")}]))
    R.append(erow(r"$K=64$ as submitted, $\lambda_E=50$", runs("reg_K64_E50_s*")))
    cap = (r"Energy penalty: accuracy at $t=1$ (B1), number of energy increases of the trained trajectory on grids of spacing $0.02$, $0.002$ and $2\times10^{-4}$ "
           r"(total over the seeds; in parentheses the number of runs with at least one increase), number of training epochs in which the penalty was nonzero (summed over the seeds), "
           r"and the largest value of the weighted penalty $\lambda_E\mathcal L_E$ during training (``---'': $\lambda_E=0$, for which $\mathcal L_E$ is not evaluated during training).")
    return wide(cap, "tab:energy", "p{4.6cm}cCCcc", H, R)


def tab_balance():
    H = row(BS + "textbf{Configuration}", m(BS + "boldsymbol{n}"), m(BS + "boldsymbol{(" + BS + "frac1T" + BS + "int" + BS + "|r" + BS + "|^2)^{1/2}}"),
            m(BS + "boldsymbol{" + BS + "frac12" + BS + "int" + BS + "|r" + BS + "|^2}"), m(BS + "boldsymbol{" + BS + "int|" + BS + "beta|}"),
            m(BS + "boldsymbol{E(0)-E(T)}"), m(BS + "boldsymbol{" + BS + "max|u_x(" + BS + "pm1)|}"), BS + "textbf{Identity}")[:-2]
    cfg = [(r"B1 revised $K=32$", fh32), (r"B1 revised $K=32$, $\lambda_E=50$", fh32E), (r"B1 revised $K=64$", fh64), (r"B1 revised $K=128$", fh128),
           (r"B1 as submitted $K=32$", so32), (r"B1 as submitted $K=64$", so64), (r"B1 soft, $\lambda_{\rm BC}=K^4$, $K=64$", vs64),
           (r"B2 revised $K=32$", b2h), (r"B3 revised $K=64$ ($T=5$)", b3K64)]
    R = []
    for lab, rr in cfg:
        rr = [r for r in rr if "balance" in r]
        if not rr:
            R.append(row(lab, DASH, "(not run)", "", "", "", "", ""))
            continue
        g = lambda k: [r["balance"][k] for r in rr]
        dE = [r["balance"]["E_first"] - r["balance"]["E_last"] for r in rr]
        R.append(row(lab, str(len(rr)), f_sci(np.mean(g("rms_residual_L2")), 1), f_sci(np.mean(g("int_half_r2")), 1), f_sci(np.mean(g("int_abs_bterm")), 1),
                     f_sci(np.mean(dE), 2), f_sci(max(g("max_abs_ux_wall")), 1), f_sci(max(g("identity_rel_err")), 0)))
    cap = (r"A-posteriori energy balance~\eqref{eq:balance-int} of the trained networks (means over the seeds): time-RMS of the residual norm, the integrated residual term, "
           r"the integrated boundary term $\beta=\varepsilon^2[u_xu_t]_{-1}^{1}$, the actual energy decrease, the largest wall slope, and the largest relative deviation "
           r"between $dE/dt$, evaluated as $\int(\varepsilon^2u_xu_{xt}+f(u)u_t)\,dx$ with $\dot{\mathbf c}$ from automatic differentiation, and the right-hand side of~\eqref{eq:balance-main} over 201 times.")
    return wide(cap, "tab:balance", "p{4.4cm}cCCCCCc", H, R)


def tab_inverse():
    out = []
    # --- table 1: joint Adam init sweep
    H = row(m(BS + "boldsymbol{" + BS + "varepsilon_{" + BS + "rm init}/" + BS + "varepsilon_{" + BS + "rm true}}"),
            r"\multicolumn{2}{c}{$\boldsymbol{\varepsilon_{\rm true}=0.05}$}", r"\multicolumn{2}{c}{$\boldsymbol{\varepsilon_{\rm true}=0.1}$}")[:-2]
    H += BS + BS + "\n" + row("", r"\textbf{as submitted}", r"\textbf{revised}", r"\textbf{as submitted}", r"\textbf{revised}")[:-2]
    R = []
    for i in (0.5, 1.0, 1.5, 2.0, 4.0):
        cells = []
        for e in (0.05, 0.1):
            for var in ("soft", "hard"):
                r = run1(f"inv_adam_{var}_e{e}_i{i}")
                cells.append(num(r["ratio_to_true"], 3) if r else DASH)
        R.append(row(num(i), *cells))
    R.append(BS + "midrule")
    cells = []
    for e in (0.05, 0.1):
        for var in ("soft", "hard"):
            rr = [run1(f"inv_adam_{var}_e{e}_i{i}") for i in (0.5, 1.0, 1.5, 2.0, 4.0)]
            rr = [r for r in rr if r]
            cells.append(m(f"{np.mean([r['ratio_to_init'] for r in rr]):.3f}" + BS + f"pm{np.std([r['ratio_to_init'] for r in rr], ddof=1):.3f}") if len(rr) > 1 else DASH)
    R.append(row(r"$\hat\varepsilon/\varepsilon_{\rm init}$", *cells))
    cap = (r"Joint Adam (procedure of the first version): estimate $\hat\varepsilon/\varepsilon_{\rm true}$ for five starting values, noise-free data, seed 42. "
           r"Last row: $\hat\varepsilon/\varepsilon_{\rm init}$, mean $\pm$ SD over the five starting values.")
    out.append(narrow(cap, "tab:inverse1", "cCCCC", H, R))
    # --- table 2: ablation of joint Adam and variable projection
    H = row(BS + "textbf{Setting}", m(BS + "boldsymbol{" + BS + "varepsilon_{" + BS + "rm true}}"), m(BS + "boldsymbol{" + BS + "hat" + BS + "varepsilon/" + BS + "varepsilon_{" + BS + "rm true}}"),
            m(BS + "boldsymbol{" + BS + "hat" + BS + "varepsilon/" + BS + "varepsilon_{" + BS + "rm init}}"), BS + "textbf{Data RMS misfit}")[:-2]
    R = [group_row(r"Joint Adam, revised network, $\varepsilon_{\rm init}=2\varepsilon_{\rm true}$", 5)]
    ab = [("default ($\\lambda_{\\rm data}=2$, 40 points, 5 times, $K=32$)", "inv_adam_hard_e0.05_i2.0"), (r"$\lambda_{\rm data}=2000$", "inv_ab_ld2000"),
          (r"$\lambda_{\rm data}=2000$, $\varepsilon_{\rm init}=\varepsilon_{\rm true}$", "inv_ab_ld2000_i1.0"), ("10 observation points", "inv_ab_nx10"),
          ("160 observation points", "inv_ab_nx160"), ("20 observation times", "inv_ab_nt20"), (r"$K=16$", "inv_ab_K16"), (r"$K=64$", "inv_ab_K64"),
          (r"$\lambda_{\rm PDE}=5$", "inv_ab_lpde5"), (r"$\lambda_{\rm PDE}=500$", "inv_ab_lpde500"),
          (r"learning-rate multiplier of $\log\varepsilon$: 1 (default 10)", "inv_ab_mult1"), (r"learning-rate multiplier of $\log\varepsilon$: 100", "inv_ab_mult100"),
          (r"16\,000 epochs (default 4000)", "inv_ab_ep16000")]
    for lab, tag in ab:
        r = run1(tag)
        R.append(row(lab, "0.05", num(r["ratio_to_true"], 3), num(r["ratio_to_init"], 3), f_sci(r["data_rms_misfit"], 1)) if r else row(lab, "0.05", "(not run)", "", ""))
    R.append(BS + "midrule")
    R.append(group_row(r"Variable projection (data-fit phase of 3000 epochs, $\lambda_{\rm data}$ as given), noise-free data", 5))
    for e in (0.05, 0.1):
        for ld in (20, 200, 2000, 20000):
            r = run1(f"inv_vp_e{e}_ld{ld}")
            R.append(row(f"$" + BS + f"lambda_{{" + BS + f"rm data}}={ld if ld < 10000 else str(ld)[:-3] + BS + "," + str(ld)[-3:]}$", f"{e:g}", num(r["ratio_to_true"], 3), DASH, f_sci(r["data_rms_misfit"], 1)) if r else row(f"ld={ld}", f"{e:g}", "(not run)", "", ""))
    for lab, tag in (("10 observation points", "inv_vp_e0.05_nx10"), ("160 observation points", "inv_vp_e0.05_nx160"), ("3 observation times", "inv_vp_e0.05_nt3"),
                     ("20 observation times", "inv_vp_e0.05_nt20"), (r"$\lambda_{\rm PDE}=5$", "inv_vp_e0.05_lpde5"), (r"$\lambda_{\rm PDE}=500$", "inv_vp_e0.05_lpde500"),
                     ("no data-fit phase", "inv_vp_e0.05_nophaseA"), (r"$K=64$", "inv_vp_e0.05_K64")):
        r = run1(tag)
        R.append(row(lab + r" ($\lambda_{\rm data}=2000$)", "0.05", num(r["ratio_to_true"], 3), DASH, f_sci(r["data_rms_misfit"], 1)) if r else row(lab, "0.05", "(not run)", "", ""))
    cap = (r"Inverse problem, noise-free data: sensitivity of the joint-Adam estimate and of the variable-projection estimate to the settings (one run each, seed 42).")
    out.append(narrow(cap, "tab:inverse2", "p{6.2cm}cccc", H, R))
    # --- table 3: noise, neural vs classical
    H = row(m(BS + "boldsymbol{" + BS + "varepsilon_{" + BS + "rm true}}"), m(BS + "boldsymbol{" + BS + "sigma}"), BS + "textbf{Var. projection (5 real.)}",
            BS + "textbf{Classical LS (same 5)}", BS + "textbf{Classical LS (25), RMSE}", BS + "textbf{Fisher s.e.}")[:-2]
    R = []
    for e in (0.01, 0.05, 0.1):
        for nz in (0.0, 0.01, 0.05):
            if e == 0.01:
                rr = [run1(f"inv_vp_e0.01_n{nz}")] if nz in (0.0, 0.01) else []
            elif nz == 0.0:
                rr = [run1(f"inv_vp_e{e}_ld2000")]
            else:
                rr = [run1(f"inv_vpnoise_e{e}_n{nz}_r{s}") for s in range(1, 6)]
            rr = [r for r in rr if r]
            if rr:
                v = [r["rel_err_pct"] for r in rr]
                neu = m(num(np.mean(v)) + (BS + "pm" + num(np.std(v, ddof=1)) if len(v) > 1 else "")) + BS + "%"
            else:
                neu = DASH
            if CL2_INV:
                c = CL2_INV[str(e)]["runs"][str(nz)]
                if nz > 0:
                    p5 = [100 * abs(c["eps_hat"][str(s)] - e) / e for s in range(1, 6)]
                    cl5 = m(num(np.mean(p5)) + BS + "pm" + num(np.std(p5, ddof=1))) + BS + "%"
                    cl25 = pct(c["rel_rmse_pct"]); fis = pct(c["fisher_rel_se_pct"])
                else:
                    cl5, cl25, fis = pct(c["rel_err_pct_mean"], 1) if c["rel_err_pct_mean"] > 1e-3 else m("<10^{-3}") + BS + "%", DASH, DASH
            else:
                cl5 = cl25 = fis = DASH
            R.append(row(f"{e:g}", f"{nz:g}", neu, cl5, cl25, fis))
        R.append(BS + "midrule") if e != 0.1 else None
    cap = (r"Inverse problem with noisy data: relative error $|\hat\varepsilon-\varepsilon_{\rm true}|/\varepsilon_{\rm true}$ of the variable-projection Spectral-PINN "
           r"($\lambda_{\rm data}=2000$; mean $\pm$ SD over five noise realizations, one run for $\varepsilon_{\rm true}=0.01$ and for noise-free data) and of the classical least-squares estimator on the same "
           r"realizations and on 25 realizations (root-mean-square error); Fisher s.e.: Cram\'er--Rao standard error $\sigma/\|\partial u_{\rm obs}/\partial\log\varepsilon\|$ of $\log\hat\varepsilon$.")
    out.append(wide(cap, "tab:inverse3", "ccCCCC", H, R))
    return "\n".join(out)


def tab_2d():
    H = row(m(BS + "boldsymbol{t}"), BS + "textbf{Fourier-Spectral-PINN}", BS + "textbf{Vanilla PINN (5 seeds)}")[:-2]
    R = []
    for t in ("0.0", "0.2", "0.5", "1.0"):
        v = [r["violations"][t] for r in ms2["results"]]
        R.append(row(t, f_sci(ex2["violations"][t]["fourier_spectral"], 1), f_stat(v)))
    cap = (r"Reflection-symmetry violation $S(t)$ (Section~\ref{sec:2d}) of the Fourier-Spectral-PINN (one run; $M=8$, $\varepsilon=0.3$) and of a vanilla PINN with the same "
           r"periodic input encoding (mean $\pm$ SD over five seeds). Data of the first version, unchanged.")
    a = narrow(cap, "tab:2d", "cCC", H, R)
    H = row(m(BS + "boldsymbol{M}"), m(BS + "boldsymbol{" + BS + "varepsilon}"), BS + "textbf{Epochs}", m(BS + "boldsymbol{S(1)}"), m(BS + "boldsymbol{" + BS + "sqrt{S(1)}}"))[:-2]
    R = [row(str(r["M"]), f"{r['eps']:g}", num(r["epochs"], 4), f_sci(r["violations"]["1.0"], 2), f_sci(np.sqrt(r["violations"]["1.0"]), 1)) for r in st2["results"]]
    cap = (r"Fourier-Spectral-PINN: reflection-symmetry violation at $t=1$ for larger mode numbers $M$ and a sharper interface ($\varepsilon=0.1$). "
           r"$\sqrt{S}$ is the root-mean-square asymmetry, to be compared with the machine epsilon $2.2\times10^{-16}$.")
    b = narrow(cap, "tab:stress-2d", "cCCCC", H, R)
    return a + "\n" + b


def tab_orig():
    H = row(BS + "textbf{Method}", m(BS + "boldsymbol{n}"), r"\textbf{Original IC} $\boldsymbol{x^2\cos(\pi x)}$", r"\textbf{Compatible IC}~\eqref{eq:IC}")[:-2]
    R = []
    for eps in ("0.01", "0.1"):
        if eps == "0.1":
            continue
        c1, c2 = cheb("orig", eps, 80, 5e-4), cheb("compat", eps, 80, 5e-4)
        R.append(row(r"Chebyshev, $N=80$, $\delta t=5\times10^{-4}$", DASH, f_sci(c1["l2"]), f_sci(c2["l2"])))
        R.append(row(r"Reaction-only ODE ($\varepsilon=0$)", DASH, f_sci(CL["orig_eps" + eps]["ode"]["1.0"]["l2"]), f_sci(CL["compat_eps" + eps]["ode"]["1.0"]["l2"])))
    for lab, po, pc in ((r"Spectral-PINN as submitted, $K=32$", "fwd_orig_K32_s*", "fwd_compat_K32_s*"), (r"Spectral-PINN as submitted, $K=64$", "fwd_orig_K64_s*", "fwd_compat_K64_s*")):
        ro, rc = runs(po), runs(pc)
        sc = set(seed(r) for r in ro)
        rc = [r for r in rc if seed(r) in sc]
        R.append(row(lab, str(len(ro)), f_stat(l2(ro)) if ro else "(not run)", f_stat(l2(rc)) if rc else DASH))
    cap = (r"Original (Neumann-incompatible) versus compatible initial condition, B1 ($\varepsilon=0.01$, $T=1$): $L^2$ error at $t=1$, mean $\pm$ SD over the same seeds "
           r"(42, 123, 7), each scored against the certified reference for its own initial condition.")
    return narrow(cap, "tab:orig", "p{5.4cm}cCC", H, R)


if __name__ == "__main__":
    tabs = {"tab_forward.tex": tab_forward, "tab_K.tex": tab_K, "tab_energy.tex": tab_energy, "tab_balance.tex": tab_balance,
            "tab_inverse.tex": tab_inverse, "tab_2d.tex": tab_2d, "tab_orig.tex": tab_orig}
    for fn, f in tabs.items():
        try:
            open(os.path.join(OUT_DIR, fn), "w", encoding="utf-8").write("% generated by make_tex_r1.py -- do not edit\n" + f())
        except Exception as ex:
            import traceback
            print("TABLE FAILED", fn, repr(ex)); traceback.print_exc()
    lines = ["% generated by make_tex_r1.py -- do not edit",
             BS + "makeatletter",
             BS + "newcommand{" + BS + "N}[1]{" + BS + "@ifundefined{num@#1}{" + BS + "textbf{??[" + BS + "detokenize{#1}]}}{" + BS + "csname num@#1" + BS + "endcsname}}",
             BS + "makeatother"]
    for k, v in sorted(NUM.items()):
        lines.append(BS + "expandafter" + BS + "def" + BS + "csname num@" + k + BS + "endcsname{" + v + "}")
    open(os.path.join(OUT_DIR, "numbers_r1.tex"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
    json.dump(NUM, open(os.path.join(ROOT, "data_r1", "numbers.json"), "w"), indent=1)
    used = set()
    for f in glob.glob(os.path.join(OUT_DIR, "*.tex")):
        if os.path.basename(f).startswith(("numbers_r1", "tab_")) or f.endswith("_v1.tex"):
            continue
        used |= set(re.findall(re.escape(BS) + r"N\{([^}]+)\}", open(f, encoding="utf-8").read()))
    missing = sorted(k for k in used if k not in NUM)
    print(len(NUM), "numbers written;", len(used), "keys used in the text;", "missing:", missing)
    for line in REPORT:
        print(" *", line)
