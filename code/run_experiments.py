"""Driver script: runs one named experiment and saves results to ../data/*.json"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import sys
import json
import time
import numpy as np
import torch

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)


def save_json(name, obj):
    path = os.path.join(DATA_DIR, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)
    print(f"[saved] {path}")


def run_vanilla(threads=2, epochs=8000):
    torch.set_num_threads(threads)
    from vanilla_pinn import train_vanilla_pinn, eval_on_grid
    from common import reference_solution_fd, l2_error, linf_error
    eps = 0.01
    x_ref, ref_snaps = reference_solution_fd(eps)
    dx_ref = x_ref[1] - x_ref[0]
    net, hist, elapsed = train_vanilla_pinn(eps, epochs=epochs, seed=42)
    x_test = np.linspace(-1, 1, 200)
    errs = {}
    for t_val in [0.2, 0.5, 0.8, 1.0]:
        u_pred = eval_on_grid(net, x_test, t_val)
        u_ref_i = np.interp(x_test, x_ref, ref_snaps[t_val])
        errs[str(t_val)] = {"l2": l2_error(u_pred, u_ref_i, x_test[1] - x_test[0]),
                             "linf": linf_error(u_pred, u_ref_i),
                             "u_pred": u_pred.tolist()}
    save_json("vanilla_pinn.json", {"elapsed": elapsed, "errors": errs, "loss_history": hist[::50],
                                     "x_test": x_test.tolist()})


def run_spectral(threads=2, epochs=6000, energy=False, tag="spectral_pinn"):
    torch.set_num_threads(threads)
    from spectral_pinn import train_spectral_pinn, eval_on_grid
    from common import reference_solution_fd, l2_error, linf_error
    eps = 0.01
    x_ref, ref_snaps = reference_solution_fd(eps)
    dx_ref = x_ref[1] - x_ref[0]
    model, hist, elapsed = train_spectral_pinn(eps, K=32, epochs=epochs, energy_constrained=energy, seed=42)
    x_test = np.linspace(-1, 1, 200)
    errs = {}
    for t_val in [0.2, 0.5, 0.8, 1.0]:
        u_pred = eval_on_grid(model, t_val, x_test)
        u_ref_i = np.interp(x_test, x_ref, ref_snaps[t_val])
        errs[str(t_val)] = {"l2": l2_error(u_pred, u_ref_i, x_test[1] - x_test[0]),
                             "linf": linf_error(u_pred, u_ref_i),
                             "u_pred": u_pred.tolist()}
    # energy trajectory on a fine time grid, for the figure
    t_grid = np.linspace(0, 1, 51)
    E_traj = []
    from chebyshev_spectral import clenshaw_curtis_weights
    with torch.no_grad():
        for tv in t_grid:
            c = model.coeffs(torch.tensor([[tv]]))
            E_traj.append(float(model.energy(c)[0]))
    n_viol = sum(1 for i in range(len(E_traj) - 1) if E_traj[i + 1] > E_traj[i] + 1e-10)
    save_json(f"{tag}.json", {"elapsed": elapsed, "errors": errs,
                               "loss_history": hist["loss"][::50],
                               "energy_history": hist["l_energy"][::50] if energy else [],
                               "x_test": x_test.tolist(),
                               "E_traj": E_traj, "t_grid": t_grid.tolist(),
                               "n_energy_violations": n_viol})


def run_chebyshev_baseline():
    from chebyshev_spectral import solve_chebyshev
    from common import reference_solution_fd, l2_error, linf_error
    eps = 0.01
    x_ref, ref_snaps = reference_solution_fd(eps)
    dx_ref = x_ref[1] - x_ref[0]
    results = []
    for N, dt in [(40, 2e-3), (60, 1e-3), (80, 5e-4)]:
        t0 = time.time()
        x, snaps = solve_chebyshev(eps, N=N, dt=dt)
        elapsed = time.time() - t0
        u1 = snaps[1.0]
        u_ref_interp = np.interp(x[::-1], x_ref, ref_snaps[1.0])[::-1]
        results.append({"N": N, "dt": dt, "l2": l2_error(u1, u_ref_interp, dx_ref),
                         "linf": linf_error(u1, u_ref_interp), "elapsed": elapsed})
    save_json("chebyshev_baseline.json", {"results": results})


def run_inverse(threads=2, epochs=4000):
    torch.set_num_threads(threads)
    from inverse_problem import recover_eps
    eps_true = 0.01
    out = {}
    for noise in [0.0, 0.01, 0.05]:
        eps_hat, hist, elapsed = recover_eps(eps_true, noise_std=noise, epochs=epochs, verbose=False)
        out[str(noise)] = {"eps_hat": eps_hat, "rel_err_pct": abs(eps_hat - eps_true) / eps_true * 100,
                            "elapsed": elapsed, "eps_history": hist[::50]}
        print(f"noise={noise}: eps_hat={eps_hat:.5f} rel_err={out[str(noise)]['rel_err_pct']:.2f}%")
    save_json("inverse_problem.json", out)


def run_2d(threads=2, epochs=3000):
    torch.set_num_threads(threads)
    from experiment_2d_symmetry import (train_fourier_spectral_pinn, train_vanilla_pinn_2d,
                                         symmetry_violation_fourier, symmetry_violation_vanilla)
    fmodel, ftime = train_fourier_spectral_pinn(epochs=epochs, verbose=False)
    vnet, vtime = train_vanilla_pinn_2d(epochs=epochs, verbose=False)
    out = {"fourier_time": ftime, "vanilla_time": vtime, "violations": {}}
    for t_val in [0.0, 0.2, 0.5, 1.0]:
        Sf = symmetry_violation_fourier(fmodel, t_val)
        Sv = symmetry_violation_vanilla(vnet, t_val)
        out["violations"][str(t_val)] = {"fourier_spectral": Sf, "vanilla": Sv}
        print(f"t={t_val}: S_fourier={Sf:.3e}  S_vanilla={Sv:.3e}")
    save_json("experiment_2d.json", out)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which == "vanilla":
        run_vanilla()
    elif which == "spectral":
        run_spectral(energy=False, tag="spectral_pinn")
    elif which == "spectral_energy":
        run_spectral(energy=True, tag="spectral_pinn_energy")
    elif which == "chebyshev":
        run_chebyshev_baseline()
    elif which == "inverse":
        run_inverse()
    elif which == "2d":
        run_2d()
    else:
        print("Unknown experiment:", which)
