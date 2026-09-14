"""
Stress axis 1: does the Spectral-PINN "break down" at large K, or is the
K=64-under-fixed-budget degradation (found in server_sweep_forward.py) fixed
by scaling the training budget with K? Tests K in {64,128} at both the
original (K=32) budget and a budget matched proportionally to K.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import json
import time
from multiprocessing import Pool, cpu_count

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

EPS = 0.01
BASE_K, BASE_WARMUP, BASE_EPOCHS = 32, 1000, 6000


def run_one(args):
    import torch
    torch.set_num_threads(3)
    from spectral_pinn import train_spectral_pinn, eval_on_grid
    from common import reference_solution_fd, l2_error, linf_error
    import numpy as np
    K, tag, scale = args
    warmup = int(BASE_WARMUP * scale)
    epochs = int(BASE_EPOCHS * scale)
    x_ref, ref_snaps = reference_solution_fd(EPS)
    dx_ref = x_ref[1] - x_ref[0]
    t0 = time.time()
    model, hist, elapsed = train_spectral_pinn(
        EPS, K=K, epochs=epochs, warmup_epochs=warmup, energy_constrained=False,
        seed=42, verbose=False)
    x_test = np.linspace(-1, 1, 200)
    u_pred = eval_on_grid(model, 1.0, x_test)
    u_ref_i = np.interp(x_test, x_ref, ref_snaps[1.0])
    r = {"K": K, "tag": tag, "scale": scale, "warmup": warmup, "epochs": epochs,
         "l2": l2_error(u_pred, u_ref_i, x_test[1] - x_test[0]),
         "linf": linf_error(u_pred, u_ref_i), "elapsed": elapsed}
    print("[done]", r, flush=True)
    return r


if __name__ == "__main__":
    jobs = [
        (64, "K64_1x", 1.0),
        (64, "K64_2x_matched", 2.0),
        (128, "K128_1x", 1.0),
        (128, "K128_4x_matched", 4.0),
        (256, "K256_1x", 1.0),
        (256, "K256_8x_matched", 8.0),
    ]
    n_proc = min(len(jobs), max(1, cpu_count() // 3))
    print(f"Running {len(jobs)} K-sweep configs on {n_proc} processes ({cpu_count()} cores)")
    t0 = time.time()
    with Pool(processes=n_proc) as pool:
        results = pool.map(run_one, jobs)
    elapsed_total = time.time() - t0
    print(f"Total wall time: {elapsed_total:.0f}s")
    with open(os.path.join(DATA_DIR, "stress_k_sweep.json"), "w") as fh:
        json.dump({"results": results, "wall_time": elapsed_total,
                   "baseline_K32": {"l2": 1.38e-3, "linf": 3.32e-3, "note": "from Table 1, K=32 mean over 5 seeds"}},
                  fh, indent=2)
    print("saved stress_k_sweep.json")
