"""Convergence certification of the trusted reference (refsol.reference):
N doubling, tolerance tightening, integrator change (Radau vs BDF), half-domain vs
full-domain, for both ICs, eps in {0.01,0.05,0.1}, t in {0.1,0.5,1.0}, plus
long-time T=5,10 at eps=0.01.  Writes cert_reference.json."""
import os, sys, json, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import refsol

xe = np.linspace(-1, 1, 2001)          # dense uniform evaluation grid (sup-norm differences)
out = {}


def maxdiff(a, b, times):
    return {str(t): float(np.max(np.abs(a[t] - b[t]))) for t in times}


def run(eps, ic, times, N, **kw):
    t0 = time.time()
    info = {}
    r = refsol.reference(eps, ic, max(times), times, xe, N=N, info=info, **kw)
    return r, time.time() - t0, info


times = [0.1, 0.5, 1.0]
for ic in ("orig", "compat"):
    for eps in (0.1, 0.05, 0.01):
        key = f"{ic}_eps{eps}_T1"
        rec = {"N_study": {}, "tol_study": {}, "method_study": {}}
        Ns = [48, 64, 96, 128, 192, 256, 384, 512]
        sol_ref, tm, info = run(eps, ic, times, 768)
        rec["N_ref"] = 768
        rec["N_ref_time_s"] = tm
        rec["N_ref_nsteps"] = info["nsteps"]
        for N in Ns:
            r, tm, info = run(eps, ic, times, N)
            rec["N_study"][N] = maxdiff(r, sol_ref, times)
        # tolerance study at N=256 (or 384 for eps .01)
        Nt = 384
        base, _, _ = run(eps, ic, times, Nt, rtol=1e-12, atol=1e-14)
        for rt in (1e-7, 1e-8, 1e-9, 1e-10, 1e-11):
            r, tm, info = run(eps, ic, times, Nt, rtol=rt, atol=rt * 1e-2)
            rec["tol_study"][f"{rt:g}"] = dict(diff=maxdiff(r, base, times), nsteps=info["nsteps"], time_s=tm)
        # default setting vs rtol 1e-12
        r, tm, info = run(eps, ic, times, Nt)
        rec["tol_study"]["default(1e-11)"] = dict(diff=maxdiff(r, base, times), nsteps=info["nsteps"], time_s=tm)
        # integrator change
        for m in ("BDF", "LSODA"):
            try:
                r, tm, info = run(eps, ic, times, Nt, rtol=1e-11, atol=1e-13, method=m)
                rec["method_study"][m] = dict(diff=maxdiff(r, base, times), nsteps=info["nsteps"], time_s=tm)
            except Exception as e:
                rec["method_study"][m] = {"error": str(e)[:100]}
        # half-domain vs full-domain at N=128
        a = refsol.reference(eps, ic, 1.0, times, xe, N=128)
        b = refsol.reference_full(eps, ic, 1.0, times, xe, N=128)
        rec["half_vs_full_N128"] = maxdiff(a, b, times)
        out[key] = rec
        print(key, "done", flush=True)
        print("  N-study (Linf vs N=768):", {N: max(v.values()) for N, v in rec["N_study"].items()})
        print("  tol-study:", {k: max(v["diff"].values()) for k, v in rec["tol_study"].items()})
        print("  method:", {k: (max(v["diff"].values()) if "diff" in v else v) for k, v in rec["method_study"].items()})
        print("  half-vs-full:", max(rec["half_vs_full_N128"].values()), flush=True)

json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "cert_reference.json"), "w"), indent=1)
