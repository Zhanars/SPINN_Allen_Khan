"""
Stress axis 2: stiffness. As eps shrinks the Allen-Cahn interface sharpens;
does the Spectral-PINN (fixed K=32, fixed budget) and the vanilla PINN
degrade gracefully or catastrophically? The FD reference grid and Chebyshev
baseline resolution are both scaled with eps so the "ground truth" stays
trustworthy at small eps.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
# numpy's BLAS backend ignores torch.set_num_threads and, left unbounded,
# spawns up to one thread per core *per process*; with 12 parallel worker
# processes this causes catastrophic oversubscription (observed on the
# 64-core server: load average > 230). Cap BLAS threads explicitly.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "2")
import json
import time
from multiprocessing import Pool, cpu_count
from stress_common import reference_nx_for_eps, chebyshev_n_for_eps

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

EPS_LIST = [0.01, 0.005, 0.002, 0.001]


def run_vanilla(eps):
    import torch
    torch.set_num_threads(3)
    from vanilla_pinn import train_vanilla_pinn, eval_on_grid
    from common import reference_solution_fd, l2_error, linf_error
    import numpy as np
    nx_ref = reference_nx_for_eps(eps)
    x_ref, ref_snaps = reference_solution_fd(eps, Nx=nx_ref)
    net, hist, elapsed = train_vanilla_pinn(eps, epochs=8000, seed=42, verbose=False)
    x_test = np.linspace(-1, 1, 200)
    u_pred = eval_on_grid(net, x_test, 1.0)
    u_ref_i = np.interp(x_test, x_ref, ref_snaps[1.0])
    r = {"method": "vanilla", "eps": eps, "nx_ref": nx_ref,
         "l2": l2_error(u_pred, u_ref_i, x_test[1] - x_test[0]),
         "linf": linf_error(u_pred, u_ref_i), "elapsed": elapsed}
    print("[done]", r, flush=True)
    return r


def run_spectral(eps):
    import torch
    torch.set_num_threads(3)
    from spectral_pinn import train_spectral_pinn, eval_on_grid
    from common import reference_solution_fd, l2_error, linf_error
    import numpy as np
    nx_ref = reference_nx_for_eps(eps)
    x_ref, ref_snaps = reference_solution_fd(eps, Nx=nx_ref)
    model, hist, elapsed = train_spectral_pinn(eps, K=32, epochs=6000, warmup_epochs=1000,
                                                energy_constrained=False, seed=42, verbose=False)
    x_test = np.linspace(-1, 1, 200)
    u_pred = eval_on_grid(model, 1.0, x_test)
    u_ref_i = np.interp(x_test, x_ref, ref_snaps[1.0])
    r = {"method": "spectral_K32", "eps": eps, "nx_ref": nx_ref,
         "l2": l2_error(u_pred, u_ref_i, x_test[1] - x_test[0]),
         "linf": linf_error(u_pred, u_ref_i), "elapsed": elapsed}
    print("[done]", r, flush=True)
    return r


def run_chebyshev(eps):
    from chebyshev_spectral import solve_chebyshev
    from common import reference_solution_fd, l2_error, linf_error
    import numpy as np
    nx_ref = reference_nx_for_eps(eps)
    x_ref, ref_snaps = reference_solution_fd(eps, Nx=nx_ref)
    N = chebyshev_n_for_eps(eps)
    dt = min(5e-4, 0.01 * eps)
    t0 = time.time()
    x, snaps = solve_chebyshev(eps, N=N, dt=dt)
    elapsed = time.time() - t0
    u1 = snaps[1.0]
    u_ref_interp = np.interp(x[::-1], x_ref, ref_snaps[1.0])[::-1]
    r = {"method": "chebyshev", "eps": eps, "N": N, "dt": dt, "nx_ref": nx_ref,
         "l2": l2_error(u1, u_ref_interp, x_ref[1] - x_ref[0]),
         "linf": linf_error(u1, u_ref_interp), "elapsed": elapsed}
    print("[done]", r, flush=True)
    return r


def dispatch(args):
    method, eps = args
    if method == "vanilla":
        return run_vanilla(eps)
    if method == "spectral":
        return run_spectral(eps)
    if method == "chebyshev":
        return run_chebyshev(eps)


if __name__ == "__main__":
    jobs = []
    for eps in EPS_LIST:
        jobs += [("vanilla", eps), ("spectral", eps), ("chebyshev", eps)]
    n_proc = min(len(jobs), max(1, cpu_count() // 3))
    print(f"Running {len(jobs)} eps-sweep configs on {n_proc} processes ({cpu_count()} cores)")
    t0 = time.time()
    with Pool(processes=n_proc) as pool:
        results = pool.map(dispatch, jobs)
    elapsed_total = time.time() - t0
    print(f"Total wall time: {elapsed_total:.0f}s")
    with open(os.path.join(DATA_DIR, "stress_eps_sweep.json"), "w") as fh:
        json.dump({"results": results, "wall_time": elapsed_total}, fh, indent=2)
    print("saved stress_eps_sweep.json")
