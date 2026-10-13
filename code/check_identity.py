"""Numerical sanity check of the energy identity for polynomial-in-x states:
dE/dt = -||u_t||^2 + <r,u_t> + eps^2 [u_x u_t]_{-1}^{1},  r = u_t - eps^2 u_xx + f(u).
Also checks dE/dt = 1/2||r||^2 - 1/2||u_t||^2 - 1/2||G||^2 + beta and the 1/4 bound.
All integrals are done with exact Gauss-Legendre quadrature (polynomial integrands)."""
import numpy as np
from numpy.polynomial import chebyshev as C
from numpy.polynomial import legendre as L

rng = np.random.default_rng(0)
K = 7
eps = 0.3
xg, wg = L.leggauss(40)  # exact for degree <= 79; integrands have degree <= 4K = 28


def integ(vals):
    return float(np.sum(wg * vals))


# smooth coefficient trajectory c_k(t) = a_k + b_k t + d_k t^2 + e_k sin(t)
a, b, d, e = (rng.standard_normal(K + 1) for _ in range(4))


def c(t):
    return a + b * t + d * t ** 2 + e * np.sin(t)


def cdot(t):
    return b + 2 * d * t + e * np.cos(t)


def u(coef, x):
    return C.chebval(x, coef)


def ux(coef, x):
    return C.chebval(x, C.chebder(coef, 1))


def uxx(coef, x):
    return C.chebval(x, C.chebder(coef, 2))


f = lambda s: s ** 3 - s
F = lambda s: 0.25 * (s ** 2 - 1) ** 2


def energy(t):
    co = c(t)
    return integ(0.5 * eps ** 2 * ux(co, xg) ** 2 + F(u(co, xg)))


def rhs(t):
    co, cd = c(t), cdot(t)
    ut = u(cd, xg)
    r = ut - eps ** 2 * uxx(co, xg) + f(u(co, xg))
    G = -eps ** 2 * uxx(co, xg) + f(u(co, xg))
    bnd = eps ** 2 * (ux(co, 1.0) * u(cd, 1.0) - ux(co, -1.0) * u(cd, -1.0))
    n_ut2 = integ(ut ** 2)
    n_r2 = integ(r ** 2)
    n_G2 = integ(G ** 2)
    ip = integ(r * ut)
    return n_ut2, n_r2, n_G2, ip, bnd


for t in [0.0, 0.4, 1.1]:
    h = 1e-5
    dE = (energy(t + h) - energy(t - h)) / (2 * h)
    n_ut2, n_r2, n_G2, ip, bnd = rhs(t)
    form1 = -n_ut2 + ip + bnd
    form2 = 0.5 * n_r2 - 0.5 * n_ut2 - 0.5 * n_G2 + bnd
    quarter = 0.25 * n_r2 + bnd
    half = 0.5 * n_r2 - 0.5 * n_ut2 + bnd
    print(f"t={t}: dE/dt(FD)={dE:.8f} form1={form1:.8f} form2={form2:.8f} | "
          f"bound 1/2-dissip={half:.4f} bound 1/4={quarter:.4f} | ok={abs(dE-form1)<1e-6 and abs(dE-form2)<1e-6 and dE<=half+1e-9 and dE<=quarter+1e-9}")
