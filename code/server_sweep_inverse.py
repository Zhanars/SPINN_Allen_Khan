"""
Large-scale server run: systematic characterization of the inverse-problem
bias(K, eps) relationship hypothesized in the manuscript (Section 4.3) --
that the eps recovery bias is caused by finite Chebyshev truncation K,
not by optimization failure. Runs eps_true x K on a process pool sized to
the machine's core count.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import sys
import json
import time
import itertools
from multiprocessing import Pool, cpu_count

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)


def run_one(args):
    import torch
    torch.set_num_threads(2)
    from inverse_problem import recover_eps
    eps_true, K, seed = args
    eps_init = eps_true * 2.0
    t0 = time.time()
    eps_hat, hist, elapsed = recover_eps(
        eps_true, K=K, noise_std=0.0, epochs=4000,
        lam_pde=50.0, lam_data=2.0, eps_init=eps_init, seed=seed, verbose=False,
    )
    rel_err = abs(eps_hat - eps_true) / eps_true * 100
    result = {"eps_true": eps_true, "K": K, "seed": seed, "eps_init": eps_init,
              "eps_hat": eps_hat, "rel_err_pct": rel_err, "elapsed": elapsed}
    print(f"[done] eps_true={eps_true} K={K} seed={seed} -> eps_hat={eps_hat:.5f} "
          f"rel_err={rel_err:.1f}% ({elapsed:.0f}s)", flush=True)
    return result


if __name__ == "__main__":
    eps_list = [0.005, 0.01, 0.02, 0.05, 0.1]
    K_list = [16, 32, 64, 96]
    seeds = [42]
    combos = [(e, k, s) for e, k, s in itertools.product(eps_list, K_list, seeds)]
    n_proc = min(len(combos), max(1, cpu_count() // 2))
    print(f"Running {len(combos)} configs on {n_proc} processes ({cpu_count()} cores total)")
    t0 = time.time()
    with Pool(processes=n_proc) as pool:
        results = pool.map(run_one, combos)
    elapsed_total = time.time() - t0
    print(f"Total wall time: {elapsed_total:.0f}s")
    outpath = os.path.join(DATA_DIR, "inverse_bias_sweep.json")
    with open(outpath, "w") as fh:
        json.dump({"results": results, "wall_time": elapsed_total}, fh, indent=2)
    print("saved", outpath)
