"""Aggregate data_r1/runs/*.json into the numbers/tables of the revised manuscript.
usage: python r1_tables.py            (prints everything; writes data_r1/summary.json)"""
import os
import re
import json
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "..", "data_r1", "runs")
OUT = os.path.join(HERE, "..", "data_r1", "summary.json")


def load(pattern):
    out = []
    for f in sorted(glob.glob(os.path.join(RUNS, pattern + ".json"))):
        out.append(json.load(open(f)))
    return out


def get(d, path):
    for p in path.split("/"):
        d = d[p]
    return d


def stat(runs, path):
    v = np.array([get(r, path) for r in runs], float)
    return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0, "n": int(len(v)),
            "min": float(v.min()), "max": float(v.max()), "values": v.tolist()}


def boot_ratio(a, b, n=20000, seed=0):
    """95% bootstrap CI of mean(a)/mean(b) with independent resampling of the seeds."""
    rng = np.random.default_rng(seed)
    a = np.asarray(a); b = np.asarray(b)
    r = [rng.choice(a, len(a)).mean() / rng.choice(b, len(b)).mean() for _ in range(n)]
    return float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))


def sci(x, d=2):
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    return f"{x / 10 ** e:.{d}f}\\times10^{{{e}}}"


def pm(s, d=2):
    e = int(np.floor(np.log10(abs(s["mean"])))) if s["mean"] != 0 else 0
    return f"({s['mean'] / 10 ** e:.{d}f}\\pm{s['sd'] / 10 ** e:.{d}f})\\times10^{{{e}}}"


def summarize(pattern, label=None):
    runs = [r for r in load(pattern) if "errors" in r]
    if not runs:
        return None
    s = {"n": len(runs), "tags": [r["job"]["tag"] for r in runs]}
    s["l2"] = stat(runs, "errors/1.0/l2")
    s["linf"] = stat(runs, "errors/1.0/linf")
    s["train_time"] = stat(runs, "train_time")
    s["max_abs_u"] = stat(runs, "max_abs_u")
    s["viol_dt0.02"] = [r["energy"]["dt0.02"]["n_increases"] for r in runs]
    if "dt0.002" in runs[0]["energy"]:
        s["viol_dt0.002"] = [r["energy"]["dt0.002"]["n_increases"] for r in runs]
        s["viol_dt0.0002"] = [r["energy"]["dt0.0002"]["n_increases"] for r in runs]
    if "penalty" in runs[0]:
        s["penalty_nonzero_epochs"] = [r["penalty"]["n_epochs_nonzero"] for r in runs]
        s["penalty_max_lamE_LE"] = [r["penalty"]["max_lamE_LE"] for r in runs]
        s["final_l_bc"] = [r["penalty"]["final_l_bc"] for r in runs]
    if "balance" in runs[0]:
        for k in ("max_abs_ux_wall", "int_half_r2", "int_abs_bterm", "bound_E_increase", "identity_rel_err", "rms_residual_L2", "E_first", "E_last"):
            s["bal_" + k] = stat([{"x": r["balance"][k]} for r in runs], "x")
    return s


if __name__ == "__main__":
    groups = {}
    for pat in sorted({re.sub(r"_s\d+$", "_s*", os.path.basename(f)[:-5]) for f in glob.glob(os.path.join(RUNS, "*.json"))}):
        s = summarize(pat)
        if s:
            groups[pat] = s
    json.dump(groups, open(OUT, "w"), indent=1)
    for k, s in groups.items():
        print(f"{k:34s} n={s['n']} L2={s['l2']['mean']:.3e}+-{s['l2']['sd']:.2e} Linf={s['linf']['mean']:.3e} t={s['train_time']['mean']:.0f}s "
              f"max|u|={s['max_abs_u']['max']:.5f} viol(0.02)={s['viol_dt0.02']}")
