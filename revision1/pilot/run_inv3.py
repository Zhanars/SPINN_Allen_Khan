import subprocess, sys, os
from multiprocessing.pool import ThreadPool
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "abl3"); os.makedirs(OUT, exist_ok=True)
cfgs = []
def add(tag, e="0.05", **kw): cfgs.append((tag, e, kw))
add("adam_lbfgs2000", lbfgs=2000)
add("lbfgs_only_3000", epochs=0, lbfgs=3000)
add("adam_epsb2_0.9", eps_beta2=0.9)
add("adam_epsb2_0.9_lbfgs", eps_beta2=0.9, lbfgs=2000)
add("lbfgs_only_i0.5", epochs=0, lbfgs=3000, eps_init_f=0.5)
add("adam_lbfgs_i0.5", lbfgs=2000, eps_init_f=0.5)
def run(c):
    tag, e, kw = c
    out = os.path.join(OUT, tag + ".json")
    if os.path.exists(out): return
    r = subprocess.run([sys.executable, os.path.join(HERE, "pilot_inv2.py"), "compat", e, out] + [f"{k}={v}" for k, v in kw.items()], capture_output=True, text=True)
    print(tag, (r.stdout.strip().splitlines() or [r.stderr[-400:]])[-1], flush=True)
ThreadPool(int(sys.argv[1])).map(run, cfgs)
