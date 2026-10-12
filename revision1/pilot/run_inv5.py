import subprocess, sys, os
from multiprocessing.pool import ThreadPool
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "abl5"); os.makedirs(OUT, exist_ok=True)
cfgs = []
def add(tag, e="0.05", **kw): cfgs.append((tag, e, kw))
for e in ("0.05", "0.1"):
    for ld in (20, 200, 2000):
        add(f"vp_ld{ld}_e{e}", e=e, mode="varpro", phaseA=3000, lam_data=ld, lam_dataA=1000)
add("vp_ld2000_lpde5_e0.05", e="0.05", mode="varpro", phaseA=3000, lam_data=2000, lam_dataA=1000, lam_pde=5)
add("adam_ld2000_e0.05_mult1000", e="0.05", mode="adam", lam_data=2000, eps_lr_mult=1, eps_beta2=0.9)
add("adam_ld2000_init1.0", e="0.05", mode="adam", lam_data=2000, eps_init_f=1.0)
def run(c):
    tag, e, kw = c
    out = os.path.join(OUT, tag + ".json")
    if os.path.exists(out): return
    r = subprocess.run([sys.executable, os.path.join(HERE, "pilot_inv2.py"), "compat", e, out] + [f"{k}={v}" for k, v in kw.items()], capture_output=True, text=True)
    last = (r.stdout.strip().splitlines() or [r.stderr[-400:]])[-1]
    print(tag, last.split("eps_hat=")[-1], flush=True)
ThreadPool(int(sys.argv[1])).map(run, cfgs)
