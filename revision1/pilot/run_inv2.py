import subprocess, sys, os
from multiprocessing.pool import ThreadPool
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "abl2"); os.makedirs(OUT, exist_ok=True)
cfgs = []
def add(tag, e="0.05", **kw): cfgs.append((tag, e, kw))
# variant screening on eps_true = 0.05, init 2x and 0.5x
for init in (2.0, 0.5):
    add(f"adam_i{init}", eps_init_f=init)
    add(f"phaseA_adam_i{init}", eps_init_f=init, phaseA=2000)
    add(f"varpro_i{init}", eps_init_f=init, mode="varpro")
    add(f"phaseA_varpro_i{init}", eps_init_f=init, mode="varpro", phaseA=2000)
def run(c):
    tag, e, kw = c
    out = os.path.join(OUT, tag + ".json")
    if os.path.exists(out): return
    r = subprocess.run([sys.executable, os.path.join(HERE, "pilot_inv2.py"), "compat", e, out] + [f"{k}={v}" for k, v in kw.items()], capture_output=True, text=True)
    print(tag, (r.stdout.strip().splitlines() or [r.stderr[-300:]])[-1], flush=True)
ThreadPool(int(sys.argv[1])).map(run, cfgs)
