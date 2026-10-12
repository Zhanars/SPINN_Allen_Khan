"""
refsol.py -- trusted reference solutions for the 1D Allen-Cahn benchmark

    u_t = eps^2 u_xx - (u^3 - u),   x in [-1, 1],   u_x(+-1, t) = 0  (t > 0),
    u(x, 0) = u0(x)

used in the revision of the Spectral-PINN manuscript (Symmetry, ID
symmetry-4596052).

Method (the "trusted" solver)
-----------------------------
* Space: Chebyshev collocation on the N+1 Gauss-Lobatto nodes
  x_j = cos(j*pi/N) (same node convention as code/chebyshev_spectral.py).
* Neumann BC: "boundary bordering" in the sense of the repo code, i.e. the
  discrete condition (D u)_0 = 0, (D u)_N = 0 holds exactly; in the
  semi-discrete system it is used to ELIMINATE the two boundary unknowns
  (equivalent to the index-1 DAE with algebraic boundary rows).  For t > 0 the
  boundary condition is therefore imposed exactly, also when u0 does not
  satisfy it (u0_orig has u0'(-1)=2, u0'(1)=-2).
* Symmetry: both initial conditions offered here are even in x, and the
  PDE and the BC commute with x -> -x, so by default the solver works on the
  half interval (even Chebyshev polynomials only; N/2 unknowns).  This is
  mathematically identical to the full-interval collocation (verified in the
  self-test and in audit_code.md) and ~8x cheaper.
* Time: scipy.integrate.solve_ivp, method 'Radau' (5th-order, L-stable,
  adaptive) with analytic Jacobian, rtol=1e-11, atol=1e-13 (default).  The
  integration is restarted at every requested output time, so the returned
  values are step endpoints (no dense-output interpolation error).
* Evaluation at arbitrary points: barycentric (Salzer/Berrut) polynomial
  interpolation through the Chebyshev data -- spectrally accurate.

Certified accuracy  ("certified" = demonstrated numerically by N-doubling, tolerance
tightening, a change of integrator and an INDEPENDENT second discretisation; it is not a
rigorous a-posteriori bound).  Max-norm differences, all numbers from audit_code.md sec. 3:
  * T<=1, eps in {0.01,0.05,0.1}, ic in {orig,compat}, t in {0.1,0.5,1}, default settings:
      - vs N=768 (rtol 1e-11) ........................ <= 4e-11  (rounding floor of the N=768 matrices)
      - vs rtol=1e-12 at N=384 ........................ <= 3.4e-13
      - Radau vs BDF / LSODA (rtol 1e-11) ............. <= 1.2e-10
      - half-interval vs full-interval solver (N=128) . <= 1.8e-13
      - independent 2nd-order ghost-node finite-difference method of lines (Radau),
        Nx=4000/8000 Richardson-extrapolated, compared at the FD nodes .... <= 5.4e-9
        (worst case ic=orig, eps=0.01, t=0.1; the other 17 cases <= 1e-10)
    => reference error <= 1e-8 (spectral convergence itself is ~1e-11) in every tested case.
  * eps=0.01, T=10 (interfaces of width ~sqrt2*eps form; N=1024 default): N=1024 vs N=1536
    (rtol 1e-12) <= 1.3e-11, BDF vs Radau <= 2e-9.  N=256 is NOT enough there (1e-2 error at t>=5).
  * eps in {0.005,0.002,0.001}, T=1, ic=orig (N=512 default): N=512 vs N=1536 <= 4e-12.
  The default N is chosen by `_default_N`; for (eps, T) outside the tested set use
  `reference_certified`, which also returns the N vs 1.5N difference as an error estimate.

API
---
u0_orig(x)   = x^2 cos(pi x)                      (submitted benchmark)
u0_compat(x) = (cos(pi x) - cos(pi x)^2)/2        (= sin^2(pi x/2) cos(pi x);
               u0' = 0 at x = +-1; same values -1 at +-1 and 0 at 0, +-0.5)
reference(eps, ic, T, times, x_eval, ...) -> {t: u(x_eval, t)}   (ic: 'orig' | 'compat' | callable)
reference_certified(...) -> (that dict, N-vs-1.5N error estimate)
l2_repo(e, x), l2_trapz(e, x), linf(e)   norms (repo-style discrete L2 / trapezoid L2 / max)
energy_cc(eps, x_nodes, u_nodes)         Ginzburg-Landau energy by Clenshaw-Curtis
self_test() -> runs in ~6 s (limit 60 s); `python refsol.py` runs it.
"""
import os

# keep BLAS from oversubscribing the machine (audit instructions: <= 2 cores)
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import time
import warnings
import numpy as np
from scipy.integrate import solve_ivp

__all__ = ["u0_orig", "u0_compat", "u0_orig_prime", "u0_compat_prime",
           "cheb_nodes_diff", "bary_interp", "reference", "reference_full",
           "reference_certified", "l2_repo", "l2_trapz", "linf",
           "energy_cc", "clenshaw_curtis_weights", "self_test", "ACCURACY_NOTE"]

ACCURACY_NOTE = ("certified Linf accuracy: <=1e-8 for T<=1, eps in {0.01,0.05,0.1}, both ICs "
                 "(spectral convergence ~1e-11; independent FD-Richardson check <=5.4e-9); "
                 "~1e-10 for eps=0.01, T<=10")


# ---------------------------------------------------------------------------
# initial conditions
# ---------------------------------------------------------------------------
def u0_orig(x):
    """Submitted benchmark IC: x^2 cos(pi x).  u0'(+-1) = -+2 (violates u_x=0)."""
    x = np.asarray(x, dtype=float)
    return x**2 * np.cos(np.pi * x)


def u0_orig_prime(x):
    x = np.asarray(x, dtype=float)
    return 2 * x * np.cos(np.pi * x) - np.pi * x**2 * np.sin(np.pi * x)


def u0_compat(x):
    """Neumann-compatible IC: (cos(pi x) - cos(pi x)^2)/2 = sin^2(pi x/2) cos(pi x)."""
    c = np.cos(np.pi * np.asarray(x, dtype=float))
    return 0.5 * (c - c * c)


def u0_compat_prime(x):
    x = np.asarray(x, dtype=float)
    c, s = np.cos(np.pi * x), np.sin(np.pi * x)
    return 0.5 * np.pi * (-s + 2 * c * s)


_IC = {"orig": u0_orig, "compat": u0_compat}


def _get_ic(ic):
    if callable(ic):
        return ic, False
    if ic in _IC:
        return _IC[ic], True          # named ICs are even in x
    raise ValueError(f"unknown ic {ic!r}; use 'orig', 'compat' or a callable")


# ---------------------------------------------------------------------------
# Chebyshev machinery
# ---------------------------------------------------------------------------
_cache = {}


def cheb_nodes_diff(N):
    """Gauss-Lobatto nodes x_j = cos(j pi/N) (j=0..N, x_0=1) and the first-order
    differentiation matrix D (Weideman-Reddy 'chebdif' style: node differences
    via sin() identities, negative-sum trick on the diagonal)."""
    if N in _cache:
        return _cache[N]
    j = np.arange(N + 1)
    th = np.pi * j / N
    x = np.cos(th)
    # x_j - x_k = -2 sin((th_j+th_k)/2) sin((th_j-th_k)/2)
    S = np.sin((th[:, None] + th[None, :]) / 2)
    Dl = np.sin((th[:, None] - th[None, :]) / 2)
    dX = -2.0 * S * Dl
    np.fill_diagonal(dX, 1.0)
    c = np.ones(N + 1)
    c[0] = c[N] = 2.0
    sgn = (-1.0) ** j
    D = (c[:, None] / c[None, :]) * (sgn[:, None] * sgn[None, :]) / dX
    np.fill_diagonal(D, 0.0)
    D[np.diag_indices(N + 1)] = -D.sum(axis=1)
    _cache[N] = (x, D)
    return x, D


def bary_interp(x_nodes, f_nodes, x_eval):
    """Barycentric interpolation (Chebyshev-Lobatto weights (-1)^j, halved at the
    ends) of data f given at the nodes x_nodes = cos(j pi/N)."""
    x_nodes = np.asarray(x_nodes, float)
    f_nodes = np.asarray(f_nodes, float)
    xe = np.atleast_1d(np.asarray(x_eval, float))
    N = len(x_nodes) - 1
    w = (-1.0) ** np.arange(N + 1)
    w[0] *= 0.5
    w[-1] *= 0.5
    out = np.empty_like(xe)
    # chunk to bound memory
    chunk = max(1, int(4e6 // (N + 1)))
    for a in range(0, len(xe), chunk):
        xs = xe[a:a + chunk]
        diff = xs[:, None] - x_nodes[None, :]
        exact = (diff == 0.0)
        hit = exact.any(axis=1)
        diff[exact] = 1.0
        tmp = w[None, :] / diff
        num = tmp @ f_nodes
        den = tmp.sum(axis=1)
        res = num / den
        if hit.any():
            idx = exact.argmax(axis=1)
            res[hit] = f_nodes[idx[hit]]
        out[a:a + chunk] = res
    return out


def clenshaw_curtis_weights(N):
    """Clenshaw-Curtis weights for x_j = cos(j pi/N) (Trefethen clencurt)."""
    theta = np.pi * np.arange(N + 1) / N
    w = np.zeros(N + 1)
    ii = np.arange(1, N)
    v = np.ones(N - 1)
    if N % 2 == 0:
        w[0] = w[N] = 1.0 / (N**2 - 1)
        for k in range(1, N // 2):
            v -= 2 * np.cos(2 * k * theta[ii]) / (4 * k**2 - 1)
        v -= np.cos(N * theta[ii]) / (N**2 - 1)
    else:
        w[0] = w[N] = 1.0 / N**2
        for k in range(1, (N - 1) // 2 + 1):
            v -= 2 * np.cos(2 * k * theta[ii]) / (4 * k**2 - 1)
    w[ii] = 2 * v / N
    return w


def energy_cc(eps, x_nodes, u_nodes):
    """Ginzburg-Landau energy  int eps^2/2 u_x^2 + (u^2-1)^2/4  dx  by
    Clenshaw-Curtis quadrature on the Chebyshev nodes (spectrally accurate
    for resolved data)."""
    N = len(x_nodes) - 1
    _, D = cheb_nodes_diff(N)
    ux = D @ u_nodes
    w = clenshaw_curtis_weights(N)
    return float(np.sum(w * (0.5 * eps**2 * ux**2 + 0.25 * (u_nodes**2 - 1.0)**2)))


# ---------------------------------------------------------------------------
# semi-discrete operators with Neumann elimination
# ---------------------------------------------------------------------------
def _build_even(N):
    """Reduced linear operator on the half interval (even solutions).
    Unknowns: w = u(x_1..x_m), m=N/2 (x_m = 0).  Returns (L, E, x) with
    full half-vector v=(u_0,w) = E w  and  (u_xx)_{1..m} = L v."""
    key = ("even", N)
    if key in _cache:
        return _cache[key]
    assert N % 2 == 0
    m = N // 2
    x, D = cheb_nodes_diff(N)
    D2 = D @ D
    P = np.zeros((N + 1, m + 1))
    for j in range(m):
        P[j, j] = 1.0
        P[N - j, j] = 1.0
    P[m, m] = 1.0
    DP = D @ P
    g = DP[0, :]                            # (D u)_0 = g . v = 0
    b = -g[1:] / g[0]
    E = np.vstack([b[None, :], np.eye(m)])   # (m+1, m): v = E w
    L = (D2 @ P)[1:m + 1, :]                 # (u_xx) at nodes 1..m
    A = L @ E                                # m x m
    out = (A, E, x, P)
    _cache[key] = out
    return out


def _build_full(N):
    key = ("full", N)
    if key in _cache:
        return _cache[key]
    x, D = cheb_nodes_diff(N)
    D2 = D @ D
    g0, gN = D[0, :], D[N, :]
    # solve  [g0[0] g0[N]; gN[0] gN[N]] [u_0;u_N] = -[g0[1:N]; gN[1:N]] w
    M2 = np.array([[g0[0], g0[N]], [gN[0], gN[N]]])
    R = -np.vstack([g0[1:N], gN[1:N]])
    B = np.linalg.solve(M2, R)               # (2, N-1)
    E = np.zeros((N + 1, N - 1))
    E[0, :] = B[0]
    E[N, :] = B[1]
    E[1:N, :] = np.eye(N - 1)
    A = D2[1:N, :] @ E
    out = (A, E, x, None)
    _cache[key] = out
    return out


def _integrate(eps, A, w0, times, rtol, atol, method):
    """Integrate w' = eps^2 A w - (w^3 - w) and return list of w at `times`
    (increasing, all > 0), restarting at each output time."""
    EA = eps**2 * A
    n = len(w0)

    def rhs(t, w):
        return EA @ w - (w**3 - w)

    def jac(t, w):
        J = EA.copy()
        J[np.diag_indices(n)] += 1.0 - 3.0 * w**2
        return J

    w = w0.copy()
    t0 = 0.0
    outs = []
    nfev = 0
    nsteps_est = 0
    for t1 in times:
        if t1 <= t0:
            outs.append(w.copy())
            continue
        kw = dict(method=method, rtol=rtol, atol=atol)
        if method in ("Radau", "BDF", "LSODA"):
            kw["jac"] = jac
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solve_ivp(rhs, (t0, t1), w, **kw)
        if not sol.success:
            raise RuntimeError(f"solve_ivp failed: {sol.message}")
        w = sol.y[:, -1].copy()
        nfev += sol.nfev
        nsteps_est += len(sol.t) - 1
        outs.append(w.copy())
        t0 = t1
    return outs, dict(nfev=nfev, nsteps=nsteps_est)


def _default_N(eps, T):
    """Default resolution (convergence tests: audit_code.md section 3).
    T<=1 : the solution is smooth in the interior; only the thin wall layer created by the
           IC/BC conflict needs resolution.   T>1 : interfaces of width ~sqrt2*eps may form,
           N ~ 10.24/eps (N=1024 for eps=0.01) is needed for ~1e-10."""
    if T <= 1.0:
        if eps >= 0.05:
            return 128
        if eps >= 0.02:
            return 192
        if eps >= 0.01:
            return 256
        return 512
    return int(min(2048, max(256, 2 * np.ceil(5.12 / eps))))


def reference(eps, ic, T, times, x_eval, N=None, rtol=1e-11, atol=1e-13,
              method="Radau", symmetric=None, return_nodes=False, info=None):
    """Trusted solution of u_t = eps^2 u_xx - (u^3-u), u_x(+-1,t)=0, u(.,0)=u0.

    Parameters
    ----------
    eps : float
    ic  : 'orig' | 'compat' | callable (callable => full-interval solver)
    T   : final time (max of `times` is used; T only selects default N)
    times : iterable of output times (0 allowed => returns u0(x_eval))
    x_eval : points in [-1, 1]
    N   : Chebyshev degree (even); default by `_default_N`
    symmetric : None -> True for named ICs (even in x), False for callables

    Returns dict {t: array(len(x_eval))}  (and optionally nodes/data).
    """
    u0fun, even = _get_ic(ic)
    if symmetric is None:
        symmetric = even
    if N is None:
        N = _default_N(eps, T)
    if symmetric and N % 2:
        N += 1
    times = [float(t) for t in times]
    ts = sorted(set(t for t in times if t > 0))
    x_eval = np.asarray(x_eval, float)
    if symmetric:
        A, E, x, P = _build_even(N)
        m = N // 2
        w0 = u0fun(x[1:m + 1])
        outs, st = _integrate(eps, A, w0, ts, rtol, atol, method)
        full = [P @ (E @ w) for w in outs]
    else:
        A, E, x, _ = _build_full(N)
        w0 = u0fun(x[1:N])
        outs, st = _integrate(eps, A, w0, ts, rtol, atol, method)
        full = [E @ w for w in outs]
    if info is not None:
        info.update(st)
        info["N"] = N
    res = {}
    nodes = {}
    for t, uf in zip(ts, full):
        res[t] = bary_interp(x, uf, x_eval)
        nodes[t] = uf
    for t in times:
        if t == 0.0:
            res[0.0] = np.asarray(u0fun(x_eval), float)
            nodes[0.0] = np.asarray(u0fun(x), float)
    if return_nodes:
        return res, x, nodes
    return res


def reference_certified(eps, ic, T, times, x_eval, N=None, **kw):
    """`reference` plus an error estimate: returns (res, err_est) where err_est is the max-norm
    difference between the solutions with N and 1.5*N Chebyshev points (spectral convergence
    makes this a (very) conservative estimate of the N-error; the time-integration error is
    controlled by rtol/atol and is ~1e-11)."""
    _, even = _get_ic(ic)
    if N is None:
        N = _default_N(eps, T)
    N2 = int(2 * round(0.75 * N))
    a = reference(eps, ic, T, times, x_eval, N=N, **kw)
    b = reference(eps, ic, T, times, x_eval, N=N2, **kw)
    err = max(float(np.max(np.abs(a[t] - b[t]))) for t in a)
    return a, err


def l2_repo(e, x):
    """The repo's discrete L2 norm sqrt(sum e^2 * dx) on a uniform grid (dx = x[1]-x[0])."""
    e = np.asarray(e, float)
    return float(np.sqrt(np.sum(e**2) * (x[1] - x[0])))


def l2_trapz(e, x):
    """Trapezoid-rule approximation of the continuous L2 norm sqrt(int e^2 dx)."""
    return float(np.sqrt(np.trapezoid(np.asarray(e, float)**2, x)))


def linf(e):
    return float(np.max(np.abs(e)))


def _fd_mol_ghost(eps, u0fun, Nx, times, rtol=1e-10, atol=1e-12):
    """Independent check solver: 2nd-order ghost-node finite differences + Radau (sparse)."""
    import scipy.sparse as sp
    h = 2.0 / Nx
    x = np.linspace(-1, 1, Nx + 1)
    n = Nx + 1
    up = np.ones(n - 1)
    lo = np.ones(n - 1)
    up[0] = 2.0
    lo[-1] = 2.0
    L = sp.diags([lo, -2.0 * np.ones(n), up], [-1, 0, 1], format="csc") * (eps**2 / h**2)

    def rhs(t, u):
        return L @ u - (u**3 - u)

    def jac(t, u):
        return L + sp.diags(1.0 - 3.0 * u**2, 0, format="csc")
    u = u0fun(x)
    t0 = 0.0
    out = {}
    for t1 in times:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solve_ivp(rhs, (t0, t1), u, method="Radau", rtol=rtol, atol=atol, jac=jac)
        u = sol.y[:, -1]
        out[t1] = u.copy()
        t0 = t1
    return x, out


def reference_full(eps, ic, T, times, x_eval, **kw):
    """Same as `reference` with the full-interval (non-symmetric) solver."""
    return reference(eps, ic, T, times, x_eval, symmetric=False, **kw)


# ---------------------------------------------------------------------------
# self test (< 60 s)
# ---------------------------------------------------------------------------
def self_test(verbose=True):
    t_start = time.time()
    log = print if verbose else (lambda *a, **k: None)
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        log(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")

    # 1. IC derivative values
    check("u0_orig'(-1)=2", abs(u0_orig_prime(-1.0) - 2.0) < 1e-13)
    check("u0_orig'(+1)=-2", abs(u0_orig_prime(1.0) + 2.0) < 1e-13)
    check("u0_compat'(+-1)=0", abs(u0_compat_prime(1.0)) < 1e-12 and abs(u0_compat_prime(-1.0)) < 1e-12)
    xx = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    check("u0_compat values (-1,0,0,0,-1)", np.allclose(u0_compat(xx), [-1, 0, 0, 0, -1], atol=1e-15))
    check("u0_compat = sin^2(pi x/2) cos(pi x)",
          np.allclose(u0_compat(xx), np.sin(np.pi * xx / 2)**2 * np.cos(np.pi * xx), atol=1e-15))
    h = 1e-6
    fd = (u0_orig(-1 + h) - u0_orig(-1 - h)) / (2 * h)
    check("central difference u0_orig'(-1)~2", abs(fd - 2.0) < 1e-6, f"({fd:.8f})")

    # 2. differentiation matrix + barycentric interpolation
    x, D = cheb_nodes_diff(32)
    f = np.exp(x) * np.sin(3 * x)
    df = np.exp(x) * (np.sin(3 * x) + 3 * np.cos(3 * x))
    check("cheb D exact for smooth f (N=32)", np.max(np.abs(D @ f - df)) < 1e-9,
          f"({np.max(np.abs(D @ f - df)):.1e})")
    xe = np.linspace(-1, 1, 37)
    check("barycentric interpolation (N=32)",
          np.max(np.abs(bary_interp(x, f, xe) - np.exp(xe) * np.sin(3 * xe))) < 1e-12)
    w = clenshaw_curtis_weights(32)
    check("Clenshaw-Curtis sum(w)=2 and int e^x", abs(w.sum() - 2) < 1e-13 and
          abs(w @ np.exp(x) - (np.e - 1 / np.e)) < 1e-13)

    # 3. exact solution of the linearised problem (tiny amplitude => cubic term ~1e-12 relative):
    #    u0 = a cos(k pi x) satisfies u_x=0 at +-1 and u(t) = a exp((1 - eps^2 k^2 pi^2) t) cos(k pi x)
    xe200 = np.linspace(-1, 1, 200)
    for kk in (1, 2):
        a0 = 1e-6
        epsl = 0.1
        r = reference(epsl, lambda x, kk=kk: a0 * np.cos(kk * np.pi * x), 1.0, [0.5, 1.0], xe200,
                      N=48, rtol=1e-11, atol=1e-22)
        err = max(np.max(np.abs(r[t] - a0 * np.exp((1 - epsl**2 * kk**2 * np.pi**2) * t) * np.cos(kk * np.pi * xe200)))
                  for t in (0.5, 1.0)) / a0
        check(f"linearised exact solution (mode {kk}), relative error", err < 1e-9, f"({err:.1e})")
    # 3b. N-convergence for both ICs
    times = [0.1, 0.5, 1.0]
    for ic in ("orig", "compat"):
        for eps in (0.1, 0.01):
            a = reference(eps, ic, 1.0, times, xe200, N=96, rtol=1e-10, atol=1e-12)
            b = reference(eps, ic, 1.0, times, xe200, N=192, rtol=1e-11, atol=1e-13)
            d = max(np.max(np.abs(a[t] - b[t])) for t in times)
            tol = 5e-6 if eps < 0.05 else 1e-7
            check(f"N-convergence ic={ic} eps={eps} (N=96 vs 192)", d < tol, f"(max diff {d:.2e})")
    # 3c. independent check: 2nd-order ghost-node FD (method of lines), Richardson-extrapolated
    for ic, eps in (("orig", 0.05), ("compat", 0.01)):
        xa, ra = _fd_mol_ghost(eps, _IC[ic], 1000, [0.5])
        xb, rb = _fd_mol_ghost(eps, _IC[ic], 2000, [0.5])
        rich = (4 * rb[0.5][::2] - ra[0.5]) / 3.0
        tr = reference(eps, ic, 1.0, [0.5], xa, N=256)[0.5]
        d = np.max(np.abs(rich - tr))
        check(f"independent FD-Richardson vs spectral ic={ic} eps={eps}", d < 1e-8, f"({d:.1e})")
    # 3d. default N and certified wrapper (long time, interfaces)
    res, est = reference_certified(0.02, "orig", 10.0, [3.0, 10.0], xe200)
    check("certified wrapper eps=0.02 T=10 (N=512 vs 768)", est < 1e-8, f"(est {est:.1e})")
    r10, xn10, nd10 = reference(0.02, "orig", 10.0, [10.0], xe200, return_nodes=True)
    E10 = energy_cc(0.02, xn10, nd10[10.0])
    Eint = 2 * (2 * np.sqrt(2) / 3) * 0.02          # two pinned interfaces at x=+-0.5
    check("long-time energy = two-interface value 2*(2*sqrt2/3)*eps", abs(E10 - Eint) < 1e-4,
          f"(E={E10:.6f}, sharp-interface {Eint:.6f})")
    # 4. symmetric (half-domain) solver == full-interval solver
    a = reference(0.05, "orig", 1.0, times, xe200, N=64)
    b = reference_full(0.05, "orig", 1.0, times, xe200, N=64)
    d = max(np.max(np.abs(a[t] - b[t])) for t in times)
    check("half-domain == full-domain (N=64)", d < 1e-9, f"({d:.1e})")
    # 5. BC for t>0, MBP, parity
    res, xn, nodes = reference(0.05, "orig", 1.0, [0.5], xe200, N=64, return_nodes=True)
    _, Dm = cheb_nodes_diff(64)
    ux = Dm @ nodes[0.5]
    check("Neumann BC holds for t>0", abs(ux[0]) < 1e-9 and abs(ux[-1]) < 1e-9)
    check("max principle |u|<=1", np.max(np.abs(res[0.5])) <= 1.0 + 1e-12)
    check("even parity", np.max(np.abs(res[0.5] - res[0.5][::-1])) < 1e-12)
    # 6. energy decreases along the trusted solution
    res, xn, nodes = reference(0.05, "compat", 1.0, [0.0, 0.25, 0.5, 1.0], xe200, N=64, return_nodes=True)
    E = [energy_cc(0.05, xn, nodes[t]) for t in (0.0, 0.25, 0.5, 1.0)]
    check("energy monotonically decreasing", all(np.diff(E) < 0), f"(E={['%.5f' % e for e in E]})")
    # 7. Radau vs BDF
    a = reference(0.05, "orig", 1.0, [1.0], xe200, N=64, rtol=1e-11, atol=1e-13)
    b = reference(0.05, "orig", 1.0, [1.0], xe200, N=64, rtol=1e-10, atol=1e-12, method="BDF")
    d = np.max(np.abs(a[1.0] - b[1.0]))
    check("Radau vs BDF", d < 1e-7, f"({d:.1e})")
    el = time.time() - t_start
    log(f"self_test {'OK' if ok else 'FAILED'} in {el:.1f} s;  {ACCURACY_NOTE}")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if self_test() else 1)
