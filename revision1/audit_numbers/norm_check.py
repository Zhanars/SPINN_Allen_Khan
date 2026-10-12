"""Norm-consistency check for the Chebyshev-collocation row of Table 1 and the Chebyshev column of Table 4.

The experiment drivers call  l2_error(u1, u_ref_interp, dx_ref)  with u1 on the (N+1) Chebyshev-Gauss-Lobatto
nodes and dx_ref = x_ref[1]-x_ref[0] of the fine UNIFORM FD grid, i.e. L2 = sqrt(sum_j e_j^2 * dx_ref).
The PINN rows use 200 uniform points and dx=2/199 (sum*dx ~ integral over [-1,1] = 2).
This script re-runs the exact same solvers (imports the project code read-only, bytecode writing disabled) and
computes the L2 error with (A) the code's formula, (B) Clenshaw-Curtis weights at the nodes,
(C) the SAME protocol as the PINN rows (barycentric interpolation of the Chebyshev solution to the 200-point
uniform grid, Riemann sum with dx=2/199), (D) trapezoid on the 501-point reference grid.
It also estimates the error of the reference itself against finer FD solutions.
Run:  python norm_check.py > norm_check_output.txt
"""
import sys, os, time, json
sys.dont_write_bytecode = True
CODE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))
sys.path.insert(0, CODE)
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from common import (reference_solution_fd, l2_error, linf_error, cheb_diff_matrix, u0_fn, f_ac)
from chebyshev_spectral import solve_chebyshev, clenshaw_curtis_weights

DATA = os.path.abspath(os.path.join(CODE, "..", "data"))


def bary_interp(x_nodes, f_nodes, x_new):
    """Barycentric interpolation for Chebyshev-Gauss-Lobatto nodes x_j = cos(j pi/N)."""
    N = len(x_nodes) - 1
    w = (-1.0) ** np.arange(N + 1)
    w[0] *= 0.5
    w[-1] *= 0.5
    out = np.empty_like(x_new)
    for i, xv in enumerate(x_new):
        d = xv - x_nodes
        k = np.where(np.abs(d) < 1e-15)[0]
        if len(k):
            out[i] = f_nodes[k[0]]
        else:
            t = w / d
            out[i] = np.sum(t * f_nodes) / np.sum(t)
    return out


def trap(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def l2_all_norms(x_nodes, u_nodes, x_ref, u_ref_on_ref_grid, dx_ref_code):
    """Return dict of L2/Linf under the alternative protocols. x_nodes descending (x[0]=1)."""
    u_ref_nodes = np.interp(x_nodes[::-1], x_ref, u_ref_on_ref_grid)[::-1]
    e = u_nodes - u_ref_nodes
    N = len(x_nodes) - 1
    w = clenshaw_curtis_weights(N)
    A = float(np.sqrt(np.sum(e ** 2) * dx_ref_code))                    # code formula
    B = float(np.sqrt(np.sum(w * e ** 2)))                              # CC quadrature at nodes
    x_test = np.linspace(-1, 1, 200)
    u_cheb_test = bary_interp(x_nodes, u_nodes, x_test)
    e_test = u_cheb_test - np.interp(x_test, x_ref, u_ref_on_ref_grid)
    C = float(np.sqrt(np.sum(e_test ** 2) * (x_test[1] - x_test[0])))   # PINN-row protocol
    Cinf = float(np.max(np.abs(e_test)))
    u_cheb_ref = bary_interp(x_nodes, u_nodes, x_ref)
    e_ref = u_cheb_ref - u_ref_on_ref_grid
    Dd = float(np.sqrt(trap(e_ref ** 2, x_ref)))                        # trapezoid on reference grid
    Dinf = float(np.max(np.abs(e_ref)))
    return dict(A_code=A, B_cc=B, C_pinnprotocol=C, D_trap_refgrid=Dd, Linf_nodes=float(np.max(np.abs(e))), Linf_200=Cinf, Linf_refgrid=Dinf,
                weight_sum_code=len(x_nodes) * dx_ref_code)


def fast_cheb(eps, N, dt, T=1.0):
    """Same scheme as chebyshev_spectral.solve_chebyshev but LU factorised once (identical mathematics)."""
    D, x = cheb_diff_matrix(N)
    D2 = D @ D
    u = u0_fn(x).copy()
    Nt = int(round(T / dt))
    A = np.eye(N + 1) - dt * eps ** 2 * D2
    A[0, :] = D[0, :]
    A[-1, :] = D[-1, :]
    lu = lu_factor(A)
    for n in range(Nt):
        rhs = u - dt * f_ac(u)
        rhs[0] = 0.0
        rhs[-1] = 0.0
        u = lu_solve(lu, rhs)
    return x, u


def fmt(d):
    return "  ".join("%s=%.4e" % (k, v) for k, v in d.items())


if __name__ == "__main__":
    t00 = time.time()
    eps = 0.01
    print("=== PART 1: Table 1 Chebyshev row (eps=0.01, N=40/60/80), reference = FD Nx=500, dt=1e-4 ===")
    x_ref, snaps = reference_solution_fd(eps)
    dx_ref = x_ref[1] - x_ref[0]
    uref = snaps[1.0]
    print("dx_ref=%.4e, (N+1)*dx_ref for N=80: %.4f (true measure 2.0) => code norm scale sqrt(%.4f/2)=%.4f" % (dx_ref, 81 * dx_ref, 81 * dx_ref, np.sqrt(81 * dx_ref / 2)))
    xt = np.linspace(-1, 1, 200)
    print("PINN-row protocol: 200*dx = %.4f (true measure 2.0)" % (200 * (xt[1] - xt[0])))
    stored = {r["N"]: r for r in json.load(open(os.path.join(DATA, "chebyshev_baseline.json")))["results"]}
    sol = {}
    for N, dt in [(40, 2e-3), (60, 1e-3), (80, 5e-4)]:
        t0 = time.time()
        x, snp = solve_chebyshev(eps, N=N, dt=dt)
        el = time.time() - t0
        u1 = snp[1.0]
        sol[N] = (x, u1)
        r = l2_all_norms(x, u1, x_ref, uref, dx_ref)
        print("N=%d dt=%g  (rerun %.2fs; stored l2=%.4e linf=%.4e)" % (N, dt, el, stored[N]["l2"], stored[N]["linf"]))
        print("   " + fmt(r))
        print("   ratio proper(PINN-protocol)/code-norm = %.3f ; CC/code = %.3f ; trapezoid-refgrid/code = %.3f" % (r["C_pinnprotocol"] / r["A_code"], r["B_cc"] / r["A_code"], r["D_trap_refgrid"] / r["A_code"]))

    print("\n=== PART 2: how accurate is the reference itself? (Chebyshev N=80 solution and PINN-protocol errors vs finer FD) ===")
    refs = {}
    for Nx, dtr in [(500, 1e-4), (1000, 5e-5), (2000, 2.5e-5), (4000, 1e-5)]:
        t0 = time.time()
        xr, sn = reference_solution_fd(eps, Nx=Nx, dt_ref=dtr, save_times=[1.0])
        refs[(Nx, dtr)] = (xr, sn[1.0])
        print("FD Nx=%d dt=%g computed in %.1fs" % (Nx, dtr, time.time() - t0))
    xfine, ufine = refs[(4000, 1e-5)]
    # Richardson-ish: compare successive references on the common 200-pt grid and 4001 grid
    for key in [(500, 1e-4), (1000, 5e-5), (2000, 2.5e-5)]:
        xr, ur = refs[key]
        up = np.interp(xfine, xr, ur)
        e = up - ufine
        L2 = np.sqrt(trap(e ** 2, xfine))
        print("FD%s vs FD(4000,1e-5): L2(trapezoid)=%.4e  Linf=%.4e" % (key, L2, np.max(np.abs(e))))
    xt = np.linspace(-1, 1, 200)
    # estimate temporal error separately: FD(500,1e-4) vs FD(500,1e-5)
    xr2, sn2 = reference_solution_fd(eps, Nx=500, dt_ref=1e-5, save_times=[1.0])
    e = sn2[1.0] - uref
    print("FD(500,dt=1e-5) vs FD(500,dt=1e-4) [pure time-step effect at fixed dx]: L2=%.4e Linf=%.4e" % (np.sqrt(trap(e ** 2, x_ref)), np.max(np.abs(e))))
    xr3, sn3 = reference_solution_fd(eps, Nx=2000, dt_ref=1e-4, save_times=[1.0])
    e = np.interp(x_ref, xr3, sn3[1.0]) - uref
    print("FD(2000,dt=1e-4) vs FD(500,dt=1e-4) [pure spatial effect at fixed dt]: L2=%.4e Linf=%.4e" % (np.sqrt(trap(e ** 2, x_ref)), np.max(np.abs(e))))
    # Chebyshev converged solution (N=120, dt=2.5e-5 ; first order in time so also compare dt=1.25e-5)
    for N, dt in [(80, 5e-4), (80, 1e-4), (80, 2.5e-5), (120, 1.25e-5)]:
        t0 = time.time()
        x, u1 = fast_cheb(eps, N, dt)
        r = l2_all_norms(x, u1, xfine, ufine, xfine[1] - xfine[0])
        print("Chebyshev N=%d dt=%g vs FD(4000,1e-5): PINN-protocol L2=%.4e (Linf200=%.4e) | CC-nodes L2=%.4e | trap L2=%.4e | (%.1fs)" % (N, dt, r["C_pinnprotocol"], r["Linf_200"], r["B_cc"], r["D_trap_refgrid"], time.time() - t0))
    # error of the manuscript reference (FD500,1e-4) on the PINN protocol versus the fine reference
    ub = np.interp(xt, x_ref, uref) - np.interp(xt, xfine, ufine)
    print("Manuscript reference FD(500,1e-4) vs FD(4000,1e-5) on 200-pt grid: L2(PINN protocol)=%.4e Linf=%.4e" % (np.sqrt(np.sum(ub ** 2) * (xt[1] - xt[0])), np.max(np.abs(ub))))
    # re-evaluate PINN rows' reference sensitivity using stored u_pred of single-run seed 42 (spectral_pinn.json, vanilla_pinn.json)
    for nm in ("spectral_pinn.json", "vanilla_pinn.json"):
        d = json.load(open(os.path.join(DATA, nm)))
        up = np.array(d["errors"]["1.0"]["u_pred"])
        e1 = up - np.interp(xt, x_ref, uref)
        e2 = up - np.interp(xt, xfine, ufine)
        print("%s t=1: L2 vs manuscript reference %.4e (stored %.4e), vs FD(4000,1e-5) %.4e ; Linf %.4e vs %.4e" % (
            nm, np.sqrt(np.sum(e1 ** 2) * (xt[1] - xt[0])), d["errors"]["1.0"]["l2"], np.sqrt(np.sum(e2 ** 2) * (xt[1] - xt[0])), np.max(np.abs(e1)), np.max(np.abs(e2))))
        # same error under a proper trapezoid on the 200 pts
        print("    trapezoid-weight L2 on the same 200 pts: %.4e (Riemann sum used by code: %.4e)" % (np.sqrt(trap(e1 ** 2, xt)), np.sqrt(np.sum(e1 ** 2) * (xt[1] - xt[0]))))

    print("\n=== PART 3: Table 4 Chebyshev column (N=400, eps-scaled FD reference) ===")
    se = {r["eps"]: r for r in json.load(open(os.path.join(DATA, "stress_eps_sweep.json")))["results"] if r["method"] == "chebyshev"}
    from stress_common import reference_nx_for_eps, chebyshev_n_for_eps
    for e_, in [(0.01,), (0.005,), (0.002,), (0.001,)]:
        t0 = time.time()
        nx = reference_nx_for_eps(e_)
        xr, sn = reference_solution_fd(e_, Nx=nx)
        N = chebyshev_n_for_eps(e_)
        dt = min(5e-4, 0.01 * e_)
        x, u1 = fast_cheb(e_, N, dt)
        r = l2_all_norms(x, u1, xr, sn[1.0], xr[1] - xr[0])
        print("eps=%g nx_ref=%d N=%d dt=%g (stored: l2=%.4e linf=%.4e) [%.0fs]" % (e_, nx, N, dt, se[e_]["l2"], se[e_]["linf"], time.time() - t0))
        print("   " + fmt(r))
        print("   code/PINN-protocol = %.3f ; code/CC = %.3f" % (r["A_code"] / r["C_pinnprotocol"], r["A_code"] / r["B_cc"]))
    print("total %.0fs" % (time.time() - t00))
