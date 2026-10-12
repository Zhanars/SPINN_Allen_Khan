"""Quick spectral reference for the pilot (the certified module is built separately by the audit agent).
Chebyshev collocation, Neumann by elimination of the two wall nodes, Radau time integration."""
import sys, os
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "code"))
from common import cheb_diff_matrix


def u0_orig(x):
    return x**2 * np.cos(np.pi * x)


def u0_compat(x):
    return 0.5 * (np.cos(np.pi * x) - np.cos(np.pi * x)**2)


def bary_weights(N):
    w = (-1.0) ** np.arange(N + 1)
    w[0] *= 0.5
    w[-1] *= 0.5
    return w


def bary_eval(xn, vals, xe):
    w = bary_weights(len(xn) - 1)
    out = np.empty(len(xe))
    for i, x in enumerate(xe):
        d = x - xn
        j = np.where(np.abs(d) < 1e-14)[0]
        if len(j):
            out[i] = vals[j[0]]
        else:
            q = w / d
            out[i] = np.dot(q, vals) / np.sum(q)
    return out


def reference(eps, u0, times, N=256, rtol=1e-11, atol=1e-13):
    D, x = cheb_diff_matrix(N)
    D2 = D @ D
    # boundary values from  D[0,:]u = 0 , D[N,:]u = 0
    A = np.array([[D[0, 0], D[0, N]], [D[N, 0], D[N, N]]])
    Bint = -np.linalg.solve(A, np.vstack([D[0, 1:N], D[N, 1:N]]))  # (2, N-1)
    P = np.zeros((N + 1, N - 1))
    P[1:N, :] = np.eye(N - 1)
    P[0, :] = Bint[0]
    P[N, :] = Bint[1]
    L = eps**2 * (D2 @ P)[1:N, :]
    def rhs(t, v):
        return L @ v - (v**3 - v)
    def jac(t, v):
        return L - np.diag(3 * v**2 - 1)
    v0 = u0(x[1:N])
    sol = solve_ivp(rhs, (0, max(times)), v0, method="Radau", jac=jac, t_eval=times, rtol=rtol, atol=atol)
    assert sol.success
    out = {}
    for k, t in enumerate(times):
        out[float(t)] = P @ sol.y[:, k]
    return x, out


if __name__ == "__main__":
    xe = np.linspace(-1, 1, 200)
    for name, u0 in (("orig", u0_orig), ("compat", u0_compat)):
        for N in (128, 256, 384):
            x, s = reference(0.01, u0, [1.0], N=N)
            u = bary_eval(x, s[1.0], xe)
            print(name, N, "max|u(1)|", np.abs(u).max(), "u(-1)", u[0], "u(0)", u[100])
