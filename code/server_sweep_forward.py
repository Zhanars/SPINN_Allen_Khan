"""
Large-scale server run: multi-seed forward-problem accuracy for vanilla
PINN, Spectral-PINN, and Spectral-PINN+energy-constrained, to report
mean +/- std instead of a single-seed point estimate. Also probes
Spectral-PINN accuracy at a larger K=64 for the same eps=0.01 problem.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import sys
import json
import time
from multiprocessing import Pool, cpu_count

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

EPS = 0.01
SEEDS = [42, 123, 7, 2024, 99]


def run_vanilla(seed):
    import torch
    torch.set_num_threads(2)
    from vanilla_pinn import train_vanilla_pinn, eval_on_grid
    from common import reference_solution_fd, l2_error, linf_error
    import numpy as np
    x_ref, ref_snaps = reference_solution_fd(EPS)
    net, hist, elapsed = train_vanilla_pinn(EPS, epochs=8000, seed=seed, verbose=False)
    x_test = np.linspace(-1, 1, 200)
    u_pred = eval_on_grid(net, x_test, 1.0)
    u_ref_i = np.interp(x_test, x_ref, ref_snaps[1.0])
    r = {"method": "vanilla", "seed": seed,
         "l2": l2_error(u_pred, u_ref_i, x_test[1] - x_test[0]),
         "linf": linf_error(u_pred, u_ref_i), "elapsed": elapsed}
    print("[done]", r, flush=True)
    return r


def run_spectral(args):
    seed, K, energy = args
    import torch
    torch.set_num_threads(2)
    from spectral_pinn import train_spectral_pinn, eval_on_grid
    from common import reference_solution_fd, l2_error, linf_error
    import numpy as np
    x_ref, ref_snaps = reference_solution_fd(EPS)
    model, hist, elapsed = train_spectral_pinn(EPS, K=K, epochs=6000, energy_constrained=energy, seed=seed, verbose=False)
    x_test = np.linspace(-1, 1, 200)
    u_pred = eval_on_grid(model, 1.0, x_test)
    u_ref_i = np.interp(x_test, x_ref, ref_snaps[1.0])
    tag = ("spectral_energy" if energy else "spectral") + (f"_K{K}" if K != 32 else "")
    r = {"method": tag, "seed": seed, "K": K,
         "l2": l2_error(u_pred, u_ref_i, x_test[1] - x_test[0]),
         "linf": linf_error(u_pred, u_ref_i), "elapsed": elapsed}
    print("[done]", r, flush=True)
    return r


if __name__ == "__main__":
    jobs = []
    jobs += [("vanilla", s) for s in SEEDS]
    jobs += [("spectral", (s, 32, False)) for s in SEEDS]
    jobs += [("spectral_energy", (s, 32, True)) for s in SEEDS]
    jobs += [("spectral_K64", (s, 64, False)) for s in SEEDS]

    n_proc = min(len(jobs), max(1, cpu_count() // 2))
    print(f"Running {len(jobs)} jobs on {n_proc} processes ({cpu_count()} cores total)")
    t0 = time.time()
    results = []
    with Pool(processes=n_proc) as pool:
        async_results = []
        for kind, arg in jobs:
            if kind == "vanilla":
                async_results.append(pool.apply_async(run_vanilla, (arg,)))
            else:
                async_results.append(pool.apply_async(run_spectral, (arg,)))
        for ar in async_results:
            results.append(ar.get())
    elapsed_total = time.time() - t0
    print(f"Total wall time: {elapsed_total:.0f}s")
    outpath = os.path.join(DATA_DIR, "forward_multiseed.json")
    with open(outpath, "w") as fh:
        json.dump({"results": results, "wall_time": elapsed_total}, fh, indent=2)
    print("saved", outpath)
