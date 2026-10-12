import subprocess, sys, os
from multiprocessing.pool import ThreadPool
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "abl4"); os.makedirs(OUT, exist_ok=True)
cfgs = []
def add(tag, e="0.05", **kw): cfgs.append((tag, e, kw))
T20 = ",".join(str(round(0.05 * i, 3)) for i in range(1, 21)); T40 = ",".join(str(round(0.025 * i, 3)) for i in range(1, 41))
for e in ("0.05", "0.1"):
    for tn, T in (("t20", T20), ("t40", T40)):
        for nz in (0.0, 0.01):
            add(f"pAvp_{tn}_e{e}_n{nz}", e=e, mode="varpro", phaseA=3000, obs_times=T, noise=nz)
    add(f"pAvp_t5_e{e}_n0", e=e, mode="varpro", phaseA=3000)
def run(c):
    tag, e, kw = c
    out = os.path.join(OUT, tag + ".json")
    if os.path.exists(out): return
    r = subprocess.run([sys.executable, os.path.join(HERE, "pilot_inv2.py"), "compat", e, out] + [f"{k}={v}" for k, v in kw.items()], capture_output=True, text=True)
    last = (r.stdout.strip().splitlines() or [r.stderr[-400:]])[-1]
    print(tag, last.split("eps_hat=")[-1], flush=True)
ThreadPool(int(sys.argv[1])).map(run, cfgs)
