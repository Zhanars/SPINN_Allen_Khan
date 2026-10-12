"""Further certification of the long-time (T=10, eps=0.01) trusted solution: N=1024 vs N=1536, tolerance, BDF."""
import os, sys, json, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import refsol
xe = np.linspace(-1, 1, 4001)
times = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
out = {}
for ic in ("orig", "compat"):
    t0 = time.time()
    base = refsol.reference(0.01, ic, 10.0, times, xe, N=1536, rtol=1e-12, atol=1e-14)
    print(ic, "base N=1536 rtol1e-12", f"{time.time()-t0:.0f}s", flush=True)
    rec = {}
    for tag, kw in {"N1024_default": dict(N=1024), "N1024_rtol1e-9": dict(N=1024, rtol=1e-9, atol=1e-11),
                    "N1024_BDF_1e-10": dict(N=1024, rtol=1e-10, atol=1e-12, method="BDF"),
                    "N1280_default": dict(N=1280)}.items():
        t0 = time.time()
        r = refsol.reference(0.01, ic, 10.0, times, xe, **kw)
        rec[tag] = {str(t): float(np.max(np.abs(r[t] - base[t]))) for t in times}
        print(ic, tag, {k: f"{v:.1e}" for k, v in rec[tag].items()}, f"{time.time()-t0:.0f}s", flush=True)
    out[ic] = rec
json.dump(out, open(os.path.join(HERE, "long_cert2.json"), "w"), indent=1)
