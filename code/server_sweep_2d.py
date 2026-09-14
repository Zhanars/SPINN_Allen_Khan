"""
Large-scale server run: multi-seed distribution of the vanilla-PINN
reflection-symmetry violation S(t) in 2D, to complement the single-seed
number in the manuscript with mean +/- std across random initializations.
The Fourier-Spectral-PINN is not re-run per seed since its violation is
architecturally exact (~1e-31) regardless of seed.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import json
import time
from multiprocessing import Pool, cpu_count

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

SEEDS = [42, 123, 7, 2024, 99]


def run_vanilla_2d(seed):
    import torch
    torch.set_num_threads(3)
    from experiment_2d_symmetry import train_vanilla_pinn_2d, symmetry_violation_vanilla
    net, vtime = train_vanilla_pinn_2d(epochs=3000, seed=seed, verbose=False)
    viol = {str(t): symmetry_violation_vanilla(net, t) for t in [0.0, 0.2, 0.5, 1.0]}
    r = {"seed": seed, "elapsed": vtime, "violations": viol}
    print("[done]", r, flush=True)
    return r


if __name__ == "__main__":
    n_proc = min(len(SEEDS), max(1, cpu_count() // 3))
    print(f"Running {len(SEEDS)} seeds on {n_proc} processes ({cpu_count()} cores total)")
    t0 = time.time()
    with Pool(processes=n_proc) as pool:
        results = pool.map(run_vanilla_2d, SEEDS)
    elapsed_total = time.time() - t0
    print(f"Total wall time: {elapsed_total:.0f}s")
    outpath = os.path.join(DATA_DIR, "symmetry_2d_multiseed.json")
    with open(outpath, "w") as fh:
        json.dump({"results": results, "wall_time": elapsed_total}, fh, indent=2)
    print("saved", outpath)
