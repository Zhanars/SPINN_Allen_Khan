"""Round-1 experiment runner.
usage:  python r1_run.py <group> [nproc]     groups are defined in r1_jobs.py
Every job writes data_r1/runs/<tag>.json (+ data_r1/models/<tag>.pt for the forward-problem networks, Spectral-PINN and vanilla PINN); finished jobs are skipped."""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import sys
import json
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RUNS = os.path.join(HERE, "..", "data_r1", "runs")
MODELS = os.path.join(HERE, "..", "data_r1", "models")
os.makedirs(RUNS, exist_ok=True)
os.makedirs(MODELS, exist_ok=True)


def run_job(job):
    import torch
    torch.set_num_threads(job.get("threads", 1))
    import r1_core
    tag = job["tag"]
    out = os.path.join(RUNS, tag + ".json")
    if os.path.exists(out):
        return tag, "skipped"
    method = job["method"]
    ic = job.get("ic", "compat"); eps = job.get("eps", 0.01); T = float(job.get("T", 1.0)); seed = job.get("seed", 42)
    u0 = r1_core.IC[ic]
    t_wall = time.time()
    rec = {"job": job}
    if method == "spectral":
        from spectral_pinn import train_spectral_pinn
        lamE = float(job.get("lam_energy", 0.0))
        epochs = job.get("epochs", int(6000 + 1500 * (T - 1)))
        nt_col = job.get("nt_col", int(100 * T)); n_ck = job.get("n_ckpt", int(10 * T))
        model, hist, elapsed = train_spectral_pinn(
            eps, K=job.get("K", 32), epochs=epochs, nt_col=nt_col, warmup_epochs=job.get("warmup", 1000),
            energy_constrained=lamE > 0, lam_energy=lamE if lamE > 0 else 50.0, n_checkpoints=n_ck, seed=seed,
            verbose=False, u0=u0, T=T, lr=job.get("lr", 1e-3), lam_bc=job.get("lam_bc", 1.0),
            grid=job.get("grid", "uniform"), hard_bc=job.get("hard_bc", False), lbfgs_iters=job.get("lbfgs", 0))
        rec["train_time"] = elapsed
        rec["epochs"] = epochs
        le = np.array(hist["l_energy"])
        rec["penalty"] = {"lam_E": lamE, "n_epochs_nonzero": int((le > 0).sum()), "max_LE": float(le.max()),
                          "max_lamE_LE": float(lamE * le.max()), "final_loss": float(hist["loss"][-1]),
                          "final_l_pde": float(hist["l_pde"][-1]), "final_l_ic": float(hist["l_ic"][-1]), "final_l_bc": float(hist["l_bc"][-1])}
        rec["loss_every50"] = hist["loss"][::50]
        torch.save(model.net.state_dict(), os.path.join(MODELS, tag + ".pt"))
        rec.update(r1_core.eval_spectral(model, eps, ic, T=T, err_times=tuple(t for t in (0.2, 0.5, 0.8, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0) if t <= T + 1e-9),
                                         diag=job.get("diag", True)))
    elif method == "vanilla":
        from vanilla_pinn import train_vanilla_pinn
        epochs = job.get("epochs", 8000)
        net, hist, elapsed = train_vanilla_pinn(eps, seed=seed, epochs=epochs, verbose=False, u0=u0, T=T, nx_col=job.get("nx_col", 80), nt_col=job.get("nt_col", 80), lbfgs_iters=job.get("lbfgs", 0))
        rec["train_time"] = elapsed; rec["epochs"] = epochs
        rec.update(r1_core.eval_vanilla(net, eps, ic, T=T, err_times=tuple(t for t in (0.2, 0.5, 0.8, 1.0, 2.0, 3.0, 5.0) if t <= T + 1e-9)))
        torch.save(net.state_dict(), os.path.join(MODELS, tag + ".pt"))
    elif method == "inverse":
        import r1_inverse
        rec.update(r1_inverse.run(job))
    else:
        raise ValueError(method)
    rec["wall_time"] = time.time() - t_wall
    json.dump(rec, open(out, "w"), indent=1)
    return tag, "done %.0fs" % rec["wall_time"]


if __name__ == "__main__":
    import r1_jobs
    from multiprocessing import Pool
    group = sys.argv[1]
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    jobs = r1_jobs.JOBS[group]()
    print(f"group {group}: {len(jobs)} jobs on {nproc} processes", flush=True)
    t0 = time.time()
    with Pool(processes=nproc, maxtasksperchild=1) as pool:
        for tag, st in pool.imap_unordered(run_job, jobs):
            print(f"[{time.time() - t0:7.0f}s] {tag}: {st}", flush=True)
    print("finished", group, "in %.0fs" % (time.time() - t0))
