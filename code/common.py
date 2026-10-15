"""
Shared utilities for the Spectral-PINN Allen-Cahn study.
Problem: u_t = eps^2 u_xx - f(u),  x in [-1,1], t in [0,T]
Neumann BC: u_x(-1,t) = u_x(1,t) = 0
IC: u(x,0) = x^2 cos(pi x)
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import numpy as np
import torch

torch.set_default_dtype(torch.float64)

EPS_DEFAULT = 0.01
T_FINAL = 1.0


def u0_fn(x):
    return x**2 * np.cos(np.pi * x)


def f_ac(u):
    return u**3 - u


def F_ac(u):
    return (u**2 - 1)**2 / 4.0


# ----------------------------------------------------------------------
# Chebyshev machinery (coefficient-space differentiation, Trefethen-style)
# ----------------------------------------------------------------------
def cheb_diff_matrix(N):
    """Standard Chebyshev differentiation matrix + Gauss-Lobatto nodes
    x_j = cos(j*pi/N), j = 0..N  (x_0=1, x_N=-1)."""
    if N == 0:
        return np.zeros((1, 1)), np.array([1.0])
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


def cheb_coeff_derivative_matrix(K):
    """Linear operator mapping Chebyshev *coefficients* a_0..a_K of a
    degree-K series to the coefficients b_0..b_K of its derivative
    (b_K = 0). Backward recurrence, see Trefethen / Boyd."""
    # Build the linear map column by column using unit vectors, applying
    # the standard backward recurrence:
    #   b[K+1] = b[K] = 0;  b[k-1] = b[k+1] + 2*k*a[k]  for k = K,...,1;  b[0] /= 2.
    Dc = np.zeros((K + 1, K + 1))
    for col in range(K + 1):
        a = np.zeros(K + 1)
        a[col] = 1.0
        b = np.zeros(K + 1)
        for k in range(K, 0, -1):
            term = 2 * k * a[k]
            b[k - 1] = (b[k + 1] if k + 1 <= K else 0.0) + term
        b[0] = b[0] / 2.0
        Dc[:, col] = b
    return Dc


def chebyshev_vandermonde(x, K):
    """Matrix T[j,k] = T_k(x_j), via three-term recurrence."""
    x = np.asarray(x)
    Tm = np.zeros((len(x), K + 1))
    Tm[:, 0] = 1.0
    if K >= 1:
        Tm[:, 1] = x
    for k in range(1, K):
        Tm[:, k + 1] = 2 * x * Tm[:, k] - Tm[:, k - 1]
    return Tm


def cheb_boundary_deriv_vectors(K):
    """Return vectors v_plus, v_minus such that sum_k c_k * v_plus[k] = u_x(1)
    and sum_k c_k * v_minus[k] = u_x(-1), using T_k'(1)=k^2, T_k'(-1)=(-1)^{k+1} k^2."""
    k = np.arange(K + 1)
    v_plus = (k**2).astype(float)
    v_minus = ((-1.0)**(k + 1) * k**2)
    return v_plus, v_minus


def fit_chebyshev_coeffs(u0_vals, x_fine, K):
    """Least-squares projection of sampled data onto degree-K Chebyshev series."""
    Tm = chebyshev_vandermonde(x_fine, K)
    c, *_ = np.linalg.lstsq(Tm, u0_vals, rcond=None)
    return c


# ----------------------------------------------------------------------
# Fine-grid semi-implicit finite-difference reference solver
# ----------------------------------------------------------------------
def reference_solution_fd(eps, Nx=500, dt_ref=1e-4, T=T_FINAL, save_times=None):
    from scipy.linalg import solve_banded
    dx = 2.0 / Nx
    x = np.linspace(-1, 1, Nx + 1)
    u = u0_fn(x).copy()
    Nt = int(round(T / dt_ref))
    lam = eps**2 * dt_ref / dx**2
    if save_times is None:
        save_times = np.round(np.arange(0.0, T + 1e-9, 0.1), 6)
    snaps = {}
    if 0.0 in save_times:
        snaps[0.0] = u.copy()
    for n in range(Nt):
        rhs = u - dt_ref * f_ac(u)
        ab = np.zeros((3, Nx + 1))
        ab[0, 1:] = -lam
        ab[1, :] = 1 + 2 * lam
        ab[2, :-1] = -lam
        ab[1, 0] = 1 + lam
        ab[1, -1] = 1 + lam
        u = solve_banded((1, 1), ab, rhs)
        t_cur = round((n + 1) * dt_ref, 6)
        for ts in save_times:
            if abs(t_cur - ts) < dt_ref / 2 and ts not in snaps:
                snaps[ts] = u.copy()
    return x, snaps


def l2_error(u, u_ref, dx):
    return float(np.sqrt(np.sum((u - u_ref)**2) * dx))


def linf_error(u, u_ref):
    return float(np.max(np.abs(u - u_ref)))
