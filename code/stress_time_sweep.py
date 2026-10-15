"""
Stress axis 3: long time horizons. The continuous energy-dissipation
theorem holds for all t > 0; does the *trained* energy-constrained
Spectral-PINN keep dE/dt <= 0 and reasonable accuracy out to T=5, T=10,
or does error/energy-violation accumulate over a long trajectory?
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import json
import time
from multiprocessing import Pool, cpu_count

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

EPS = 0.01


def run_one(T):
    import torch
    torch.set_num_threads(4)
    import spectral_pinn
    from spectral_pinn import train_spectral_pinn, eval_on_grid
    from chebyshev_spectral import clenshaw_curtis_weights
    from common import reference_solution_fd, l2_error, linf_error
    import numpy as np

    nt_col = int(100 * T)
    n_ckpt = int(10 * T)
    epochs = int(6000 + 1500 * (T - 1))
    save_times = list(np.round(np.arange(0.0, T + 1e-9, 1.0), 3))
    x_ref, ref_snaps = reference_solution_fd(EPS, T=T, save_times=save_times)
    dx_ref = x_ref[1] - x_ref[0]

    t0 = time.time()
    model, hist, elapsed = train_spectral_pinn(
        EPS, K=32, epochs=epochs, nt_col=nt_col, warmup_epochs=1000,
        energy_constrained=True, n_checkpoints=n_ckpt, seed=42, verbose=False, T=T)

    errs = {}
    for tv in save_times:
        if tv == 0.0:
            continue
        x_test = np.linspace(-1, 1, 200)
        u_pred = eval_on_grid(model, tv, x_test)
        u_ref_i = np.interp(x_test, x_ref, ref_snaps[tv])
        errs[str(tv)] = {"l2": l2_error(u_pred, u_ref_i, x_test[1] - x_test[0]),
                          "linf": linf_error(u_pred, u_ref_i)}

    t_grid = np.linspace(0, T, int(20 * T) + 1)
    E_traj = []
    with torch.no_grad():
        for tv in t_grid:
            c = model.coeffs(torch.tensor([[tv]]))
            E_traj.append(float(model.energy(c)[0]))
    n_viol = sum(1 for i in range(len(E_traj) - 1) if E_traj[i + 1] > E_traj[i] + 1e-10)

    r = {"T": T, "epochs": epochs, "nt_col": nt_col, "n_ckpt": n_ckpt,
         "elapsed": elapsed, "errors": errs, "n_energy_violations": n_viol,
         "n_eval_points": len(E_traj), "E_first": E_traj[0], "E_last": E_traj[-1]}
    print("[done]", r, flush=True)
    return r


if __name__ == "__main__":
    T_list = [1.0, 5.0, 10.0]
    n_proc = min(len(T_list), max(1, cpu_count() // 4))
    print(f"Running {len(T_list)} T-sweep configs on {n_proc} processes ({cpu_count()} cores)")
    t0 = time.time()
    with Pool(processes=n_proc) as pool:
        results = pool.map(run_one, T_list)
    elapsed_total = time.time() - t0
    print(f"Total wall time: {elapsed_total:.0f}s")
    with open(os.path.join(DATA_DIR, "stress_time_sweep.json"), "w") as fh:
        json.dump({"results": results, "wall_time": elapsed_total}, fh, indent=2)
    print("saved stress_time_sweep.json")
