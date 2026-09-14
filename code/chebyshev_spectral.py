"""
Chebyshev collocation spectral solver for the 1D Allen-Cahn equation
(semi-implicit / IMEX time-stepping), used as a classical high-accuracy
baseline that plays the role of the spectral method of Shiyapov et al.
(ref. [3] in the manuscript) for the *local* Allen-Cahn equation.

Neumann BC u_x(+-1,t)=0 are imposed by "boundary bordering": the first and
last rows of the implicit system are replaced with the Neumann condition
written in terms of the first-derivative Chebyshev matrix D.
"""
import numpy as np
from common import cheb_diff_matrix, u0_fn, f_ac, F_ac, T_FINAL


def solve_chebyshev(eps, N=80, dt=1e-3, T=T_FINAL, save_times=None):
    D, x = cheb_diff_matrix(N)          # x[0]=1 ... x[N]=-1
    D2 = D @ D
    u = u0_fn(x).copy()
    Nt = int(round(T / dt))

    A = np.eye(N + 1) - dt * eps**2 * D2
    # Boundary bordering: replace row 0 (x=1) and row N (x=-1) with Neumann rows
    A[0, :] = D[0, :]
    A[-1, :] = D[-1, :]

    if save_times is None:
        save_times = np.round(np.arange(0.0, T + 1e-9, 0.1), 6)
    snaps = {}
    if 0.0 in save_times:
        snaps[0.0] = u.copy()

    for n in range(Nt):
        rhs = u - dt * f_ac(u)
        rhs[0] = 0.0
        rhs[-1] = 0.0
        u = np.linalg.solve(A, rhs)
        t_cur = round((n + 1) * dt, 6)
        for ts in save_times:
            if abs(t_cur - ts) < dt / 2 and ts not in snaps:
                snaps[ts] = u.copy()
    return x, snaps


def clenshaw_curtis_weights(N):
    """Weights for the Chebyshev-Gauss-Lobatto nodes x_j=cos(j*pi/N),
    exact for polynomials up to degree N (Trefethen, 'clencurt.m')."""
    theta = np.pi * np.arange(N + 1) / N
    w = np.zeros(N + 1)
    ii = np.arange(1, N)
    v = np.ones(N - 1)
    if N % 2 == 0:
        w[0] = 1.0 / (N**2 - 1)
        w[N] = w[0]
        for k in range(1, N // 2):
            v -= 2 * np.cos(2 * k * theta[ii]) / (4 * k**2 - 1)
        v -= np.cos(N * theta[ii]) / (N**2 - 1)
    else:
        w[0] = 1.0 / N**2
        w[N] = w[0]
        for k in range(1, (N - 1) // 2 + 1):
            v -= 2 * np.cos(2 * k * theta[ii]) / (4 * k**2 - 1)
    w[ii] = 2 * v / N
    return w


def energy(u, x, D, eps):
    ux = D @ u
    N = len(x) - 1
    w = clenshaw_curtis_weights(N)
    return float(np.sum(w * (eps**2 / 2 * ux**2 + F_ac(u))))


if __name__ == "__main__":
    from common import reference_solution_fd, l2_error, linf_error
    eps = 0.01
    print("Reference FD solve (Nx=500, dt=1e-4)...")
    x_ref, ref_snaps = reference_solution_fd(eps)
    dx_ref = x_ref[1] - x_ref[0]

    print(f"{'N':>4} {'dt':>8} {'L2 err (t=1)':>14} {'Linf err':>12}")
    for N, dt in [(40, 2e-3), (60, 1e-3), (80, 5e-4)]:
        x, snaps = solve_chebyshev(eps, N=N, dt=dt)
        u1 = snaps[1.0]
        u_ref_interp = np.interp(x[::-1], x_ref, ref_snaps[1.0])[::-1]
        l2 = l2_error(u1, u_ref_interp, dx_ref)
        linf = linf_error(u1, u_ref_interp)
        print(f"{N:>4} {dt:>8.1e} {l2:>14.4e} {linf:>12.4e}")
