"""Run many pilot_inv.py configurations in a process pool (one thread each). usage: python run_inv_ablation.py <set> <nproc>"""
import subprocess, sys, os, json, itertools
from multiprocessing.pool import ThreadPool

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "abl")
os.makedirs(OUT, exist_ok=True)
which = sys.argv[1]
nproc = int(sys.argv[2])
E = "0.05"
cfgs = []  # (tag, ic, data, eps_true, {kw})


def add(tag, **kw):
    ic = kw.pop("ic", "compat"); data = kw.pop("data", "spec"); e = kw.pop("eps", E)
    cfgs.append((tag, ic, data, e, kw))


if which == "set1":
    add("base")
    for f in (0.5, 1.0, 1.5, 4.0):
        add(f"init{f}", eps_init_f=f)
    for v in (20, 200, 2000):
        add(f"ldata{v}", lam_data=v)
    for v in (0.5, 5, 500):
        add(f"lpde{v}", lam_pde=v)
    add("nowarm", warmup=0)
    add("ep16000", epochs=16000)
    add("lr5e-4", lr=5e-4)
    for m in (1, 100):
        add(f"epsmult{m}", eps_lr_mult=m)
    add("obs160", n_obs_x=160)
    add("obs10", n_obs_x=10)
    add("times20", obs_times=",".join(str(round(0.05 * i, 2)) for i in range(1, 21)))
    add("K16", K=16)
    add("K64", K=64)
    add("eps0.1", eps="0.1")
    add("eps0.2", eps="0.2")
    add("eps0.01", eps="0.01")


def run(c):
    tag, ic, data, e, kw = c
    out = os.path.join(OUT, f"{which}_{tag}.json")
    if os.path.exists(out):
        return
    args = [sys.executable, os.path.join(HERE, "pilot_inv.py"), ic, data, e, out] + [f"{k}={v}" for k, v in kw.items()]
    r = subprocess.run(args, capture_output=True, text=True)
    print(tag, r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-300:], flush=True)


ThreadPool(nproc).map(run, cfgs)
