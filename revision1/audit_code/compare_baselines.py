"""Part (b)/(c): errors of
 (i)  current FD reference  (repo common.reference_solution_fd, Nx=500, dt=1e-4,
      boundary rows diag 1+lam, off-diag -lam)
 (ii) corrected ghost-node FD (same Nx, dt; boundary rows diag 1+2lam, off-diag -2lam)
 (iv) repo solve_chebyshev(N=80, dt=5e-4)
relative to (iii) the trusted spectral reference (refsol.reference).
Metrics (all on the uniform 200-point grid of the PINN evaluation, dx=2/199, FD data
linearly interpolated exactly as in the repo, Chebyshev data barycentric):
  L2_200 = sqrt(sum e^2 dx)   (repo formula)
  L2_cont = sqrt(int_{-1}^{1} e^2 dx)  on a 20001-point fine grid (same interpolation rule)
  Linf_200 = max |e| on the 200-point grid
Writes compare_baselines.json."""
import os, sys, json, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
from scipy.linalg import solve_banded
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import refsol

REPO_CODE = r"D:\Videos\KazNPU\Symmetry_PINN_2026\code"


def f_ac(u):
    return u**3 - u


# ---- FD solvers (copy of repo common.reference_solution_fd with IC and BC option) ----
def fd_solution(eps, u0fun, bc="current", Nx=500, dt=1e-4, T=1.0, save_times=None):
    dx = 2.0 / Nx
    x = np.linspace(-1, 1, Nx + 1)
    u = u0fun(x).copy()
    Nt = int(round(T / dt))
    lam = eps**2 * dt / dx**2
    if save_times is None:
        save_times = np.round(np.arange(0.0, T + 1e-9, 0.1), 6)
    save_times = [float(s) for s in save_times]
    snaps = {}
    if 0.0 in save_times:
        snaps[0.0] = u.copy()
    ab = np.zeros((3, Nx + 1))
    ab[0, 1:] = -lam
    ab[1, :] = 1 + 2 * lam
    ab[2, :-1] = -lam
    if bc == "current":
        ab[1, 0] = 1 + lam
        ab[1, -1] = 1 + lam
    elif bc == "ghost":
        ab[0, 1] = -2 * lam       # A[0,1]
        ab[2, -2] = -2 * lam      # A[N,N-1]
    else:
        raise ValueError(bc)
    pending = sorted(s for s in save_times if s > 0)
    for n in range(Nt):
        rhs = u - dt * f_ac(u)
        u = solve_banded((1, 1), ab, rhs)
        t_cur = round((n + 1) * dt, 6)
        while pending and abs(t_cur - pending[0]) < dt / 2:
            snaps[pending.pop(0)] = u.copy()
    return x, snaps


# ---- Chebyshev IMEX baseline (copy of repo solve_chebyshev with IC parameter) ----
def cheb_baseline(eps, u0fun, N=80, dt=5e-4, T=1.0, save_times=None):
    # repo uses common.cheb_diff_matrix (Trefethen); same matrix rebuilt here
    D, x = repo_cheb_diff_matrix(N)
    D2 = D @ D
    u = u0fun(x).copy()
    Nt = int(round(T / dt))
    A = np.eye(N + 1) - dt * eps**2 * D2
    A[0, :] = D[0, :]
    A[-1, :] = D[-1, :]
    if save_times is None:
        save_times = np.round(np.arange(0.0, T + 1e-9, 0.1), 6)
    save_times = [float(s) for s in save_times]
    snaps = {}
    if 0.0 in save_times:
        snaps[0.0] = u.copy()
    pending = sorted(s for s in save_times if s > 0)
    for n in range(Nt):
        rhs = u - dt * f_ac(u)
        rhs[0] = 0.0
        rhs[-1] = 0.0
        u = np.linalg.solve(A, rhs)
        t_cur = round((n + 1) * dt, 6)
        while pending and abs(t_cur - pending[0]) < dt / 2:
            snaps[pending.pop(0)] = u.copy()
    return x, snaps


def repo_cheb_diff_matrix(N):
    j = np.arange(N + 1)
    x = np.cos(np.pi * j / N)
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0
    c = c * (-1.0)**j
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


XU = np.linspace(-1, 1, 200)
DXU = XU[1] - XU[0]
XF = np.linspace(-1, 1, 20001)


def metrics_on_200(e200):
    return dict(L2_200=float(np.sqrt(np.sum(e200**2) * DXU)), Linf_200=float(np.max(np.abs(e200))))


def l2_cont(efine):
    return float(np.sqrt(np.trapezoid(efine**2, XF)))


def fd_vs_trusted(x_fd, u_fd, trusted_fn):
    """Errors of FD data (linear interpolation, as in the repo) vs the trusted solution."""
    u200 = np.interp(XU, x_fd, u_fd)
    e200 = u200 - trusted_fn(XU)
    m = metrics_on_200(e200)
    efine = np.interp(XF, x_fd, u_fd) - trusted_fn(XF)
    m["L2_cont"] = l2_cont(efine)
    m["Linf_cont"] = float(np.max(np.abs(efine)))
    # nodal error (no interpolation): trapezoid weights on the FD nodes
    en = u_fd - trusted_fn(x_fd)
    m["Linf_nodes"] = float(np.max(np.abs(en)))
    m["L2_nodes_trap"] = float(np.sqrt(np.trapezoid(en**2, x_fd)))
    return m


def cheb_vs_trusted(x_c, u_c, trusted_fn):
    u200 = refsol.bary_interp(x_c, u_c, XU)
    e200 = u200 - trusted_fn(XU)
    m = metrics_on_200(e200)
    efine = refsol.bary_interp(x_c, u_c, XF) - trusted_fn(XF)
    m["L2_cont"] = l2_cont(efine)
    m["Linf_cont"] = float(np.max(np.abs(efine)))
    N = len(x_c) - 1
    w = refsol.clenshaw_curtis_weights(N)
    en = u_c - trusted_fn(x_c)
    m["Linf_nodes"] = float(np.max(np.abs(en)))
    m["L2_nodes_CC"] = float(np.sqrt(np.sum(w * en**2)))
    return m


def run_case(eps, ic, times=(0.1, 0.5, 1.0), T=1.0, Ntrusted=256, extra=True):
    u0fun = refsol._IC[ic]
    times = list(times)
    rec = {"eps": eps, "ic": ic}
    # trusted: nodal data at times, evaluated through bary_interp
    _, xn, nodes = refsol.reference(eps, ic, T, times, XU, N=Ntrusted, return_nodes=True)
    trusted = {t: (lambda xx, t=t: refsol.bary_interp(xn, nodes[t], xx)) for t in times}
    rec["trusted_stats"] = {}
    for t in times:
        u = nodes[t]
        rec["trusted_stats"][str(t)] = dict(max=float(u.max()), min=float(u.min()),
                                            max_abs=float(np.abs(u).max()))
    # (i) current FD, (ii) ghost FD
    t0 = time.time()
    x_fd, snap_cur = fd_solution(eps, u0fun, "current", T=T, save_times=times)
    t_cur = time.time() - t0
    x_fd, snap_gh = fd_solution(eps, u0fun, "ghost", T=T, save_times=times)
    # (iv) Chebyshev N=80 dt=5e-4
    t0 = time.time()
    x_c, snap_c = cheb_baseline(eps, u0fun, N=80, dt=5e-4, T=T, save_times=times)
    t_c = time.time() - t0
    rec["time_fd_s"] = t_cur
    rec["time_cheb_s"] = t_c
    extras = {}
    if extra:
        # diagnostics to separate time vs space error of the FD schemes
        _, snap_gh_dt = fd_solution(eps, u0fun, "ghost", T=T, dt=2.5e-5, save_times=times)
        x2, snap_gh_nx = fd_solution(eps, u0fun, "ghost", Nx=2000, T=T, dt=1e-4, save_times=times)
        x3, snap_cur_dt = fd_solution(eps, u0fun, "current", T=T, dt=2.5e-5, save_times=times)
        x4, snap_cur_nx = fd_solution(eps, u0fun, "current", Nx=2000, T=T, dt=1e-4, save_times=times)
    for t in times:
        r = {}
        r["(i)_FD_current"] = fd_vs_trusted(x_fd, snap_cur[t], trusted[t])
        r["(ii)_FD_ghost"] = fd_vs_trusted(x_fd, snap_gh[t], trusted[t])
        r["(iv)_Cheb_N80_dt5e-4"] = cheb_vs_trusted(x_c, snap_c[t], trusted[t])
        # pairwise: current vs ghost on the 200 grid
        e = np.interp(XU, x_fd, snap_cur[t]) - np.interp(XU, x_fd, snap_gh[t])
        r["(i)-(ii)"] = metrics_on_200(e)
        # interpolation-only floor: trusted sampled at FD nodes then linearly interpolated
        e = np.interp(XU, x_fd, trusted[t](x_fd)) - trusted[t](XU)
        r["linear_interp_floor_Nx500"] = metrics_on_200(e)
        if extra:
            r["diag_ghost_dt2.5e-5"] = fd_vs_trusted(x_fd, snap_gh_dt[t], trusted[t])
            r["diag_ghost_Nx2000_dt1e-4"] = fd_vs_trusted(x2, snap_gh_nx[t], trusted[t])
            r["diag_current_dt2.5e-5"] = fd_vs_trusted(x3, snap_cur_dt[t], trusted[t])
            r["diag_current_Nx2000_dt1e-4"] = fd_vs_trusted(x4, snap_cur_nx[t], trusted[t])
        # the repo's own (mismatched) evaluation of the Chebyshev baseline against the repo FD ref:
        ref_i = np.interp(x_c[::-1], x_fd, snap_cur[t])[::-1]
        dx_ref = x_fd[1] - x_fd[0]
        e = snap_c[t] - ref_i
        r["repo_style_cheb_vs_FDcurrent"] = dict(L2=float(np.sqrt(np.sum(e**2) * dx_ref)),
                                                 Linf=float(np.max(np.abs(e))))
        rec[str(t)] = r
    return rec


if __name__ == "__main__":
    out = {}
    for ic in ("orig", "compat"):
        for eps in (0.01, 0.05, 0.1):
            t0 = time.time()
            out[f"{ic}_eps{eps}"] = run_case(eps, ic)
            print(ic, eps, f"{time.time()-t0:.0f}s", flush=True)
            for t in ("0.1", "0.5", "1.0"):
                r = out[f"{ic}_eps{eps}"][t]
                print("  t=", t, {k: (round(v["L2_200"], 8), round(v["Linf_200"], 8)) for k, v in r.items()
                                  if isinstance(v, dict) and "L2_200" in v})
    json.dump(out, open(os.path.join(HERE, "compare_baselines.json"), "w"), indent=1)
