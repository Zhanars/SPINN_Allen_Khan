"""Check that the FD / Chebyshev replicas in compare_baselines.py reproduce the repo
functions (orig IC), and reproduce the numbers stored in data/chebyshev_baseline.json."""
import os, sys, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, r"D:\Videos\KazNPU\Symmetry_PINN_2026\code")
import common, chebyshev_spectral
sys.path.insert(0, HERE)
import compare_baselines as cb, refsol
import torch; torch.set_num_threads(2)

x_r, s_r = common.reference_solution_fd(0.01)
x_m, s_m = cb.fd_solution(0.01, refsol.u0_orig, "current", save_times=np.round(np.arange(0, 1.0 + 1e-9, 0.1), 6))
print("FD current replica max diff over all saved t:", max(np.max(np.abs(s_r[k] - s_m[k])) for k in s_m))
print("n snaps repo/mine", len(s_r), len(s_m))
x_c, s_c = chebyshev_spectral.solve_chebyshev(0.01, N=80, dt=5e-4)
x_c2, s_c2 = cb.cheb_baseline(0.01, refsol.u0_orig, 80, 5e-4)
print("Cheb replica max diff:", max(np.max(np.abs(s_c[k] - s_c2[k])) for k in s_c2), "x diff", np.max(np.abs(x_c-x_c2)))
# repo numbers
ref_i = np.interp(x_c[::-1], x_r, s_r[1.0])[::-1]
dx_ref = x_r[1] - x_r[0]
print("repo-style N=80:", common.l2_error(s_c[1.0], ref_i, dx_ref), common.linf_error(s_c[1.0], ref_i))
print("json            :", json.load(open(r"D:\Videos\KazNPU\Symmetry_PINN_2026\data\chebyshev_baseline.json"))["results"][2])
# repo u0 vs refsol u0
x = np.linspace(-1, 1, 1001)
print("u0 identical:", np.max(np.abs(common.u0_fn(x) - refsol.u0_orig(x))))
