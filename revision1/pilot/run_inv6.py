import subprocess, sys, os
from multiprocessing.pool import ThreadPool
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "abl6"); os.makedirs(OUT, exist_ok=True)
cfgs = []
def add(tag, e="0.05", **kw): cfgs.append((tag, e, kw))
H = dict(grid="cgl", hard="1")
for e in ("0.05", "0.1"):
    add(f"hard_adam_i2_e{e}", e=e, **H)
    add(f"hard_adam_i1_e{e}", e=e, eps_init_f=1.0, **H)
    add(f"hard_vp_ld200_e{e}", e=e, mode="varpro", phaseA=3000, lam_data=200, lam_dataA=1000, **H)
    add(f"hard_vp_ld2000_e{e}", e=e, mode="varpro", phaseA=3000, lam_data=2000, lam_dataA=1000, **H)
def run(c):
    tag, e, kw = c
    out = os.path.join(OUT, tag + ".json")
    if os.path.exists(out): return
    r = subprocess.run([sys.executable, os.path.join(HERE, "pilot_inv2.py"), "compat", e, out] + [f"{k}={v}" for k, v in kw.items()], capture_output=True, text=True)
    last = (r.stdout.strip().splitlines() or [r.stderr[-400:]])[-1]
    print(tag, last.split("eps_hat=")[-1], flush=True)
ThreadPool(int(sys.argv[1])).map(run, cfgs)
