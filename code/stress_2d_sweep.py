"""
Stress axis 4: does the Fourier-Spectral-PINN's *architectural* (exact)
symmetry guarantee survive at larger mode counts M and sharper (smaller
eps) 2D interfaces, and how does its PDE training behave there?
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import json
import time
from multiprocessing import Pool, cpu_count

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)


def run_one(args):
    import torch
    torch.set_num_threads(3)
    from experiment_2d_symmetry import train_fourier_spectral_pinn, symmetry_violation_fourier
    M, eps, n_grid, epochs = args
    t0 = time.time()
    model, ftime = train_fourier_spectral_pinn(eps=eps, M=M, epochs=epochs, warmup_epochs=500,
                                                n_grid=n_grid, verbose=False)
    viol = {str(t): symmetry_violation_fourier(model, t) for t in [0.0, 0.5, 1.0]}
    r = {"M": M, "eps": eps, "n_grid": n_grid, "epochs": epochs,
         "elapsed": ftime, "violations": viol}
    print("[done]", r, flush=True)
    return r


if __name__ == "__main__":
    jobs = [
        (8, 0.3, 48, 3000),    # baseline (already have data, re-run for consistency)
        (16, 0.3, 64, 3000),   # larger M, same budget
        (16, 0.3, 64, 6000),   # larger M, scaled budget
        (8, 0.1, 48, 3000),    # sharper interface, baseline M
        (16, 0.1, 64, 6000),   # sharper interface, larger M, scaled budget
    ]
    n_proc = min(len(jobs), max(1, cpu_count() // 3))
    print(f"Running {len(jobs)} 2D-scale configs on {n_proc} processes ({cpu_count()} cores)")
    t0 = time.time()
    with Pool(processes=n_proc) as pool:
        results = pool.map(run_one, jobs)
    elapsed_total = time.time() - t0
    print(f"Total wall time: {elapsed_total:.0f}s")
    with open(os.path.join(DATA_DIR, "stress_2d_sweep.json"), "w") as fh:
        json.dump({"results": results, "wall_time": elapsed_total}, fh, indent=2)
    print("saved stress_2d_sweep.json")
