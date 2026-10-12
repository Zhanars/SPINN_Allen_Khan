"""Convergence study of the trusted reference at eps=0.01 for long times (T up to 10)."""
import os, sys, json, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import refsol
xe = np.linspace(-1, 1, 4001)
times = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
out = {}
for ic in ("orig", "compat"):
    sols = {}
    for N in (256, 384, 512, 768, 1024):
        info = {}
        t0 = time.time()
        sols[N] = refsol.reference(0.01, ic, 10.0, times, xe, N=N, info=info)
        print(ic, N, f"{time.time()-t0:.1f}s", info, flush=True)
    base = sols[1024]
    out[ic] = {str(N): {str(t): float(np.max(np.abs(sols[N][t] - base[t]))) for t in times} for N in sols if N != 1024}
    for N in out[ic]:
        print(ic, "N=", N, "vs 1024:", {t: f"{v:.2e}" for t, v in out[ic][N].items()}, flush=True)
json.dump(out, open(os.path.join(HERE, "long_cert.json"), "w"), indent=1)
