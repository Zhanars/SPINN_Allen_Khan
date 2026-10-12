"""Independent verification of the trusted spectral reference by a *different* spatial discretisation:
second-order ghost-node finite differences (method of lines) integrated with Radau (sparse Jacobian,
rtol 1e-10), Nx = 1000, 2000, 4000, 8000, nodal comparison with the trusted solution (no interpolation)
and Richardson extrapolation.  Writes independent_check.json."""
import os, sys, json, time, warnings
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
import scipy.sparse as sp
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import refsol


def fd_mol(eps, u0fun, Nx, times, rtol=1e-10, atol=1e-12):
    h = 2.0 / Nx
    x = np.linspace(-1, 1, Nx + 1)
    n = Nx + 1
    main = -2.0 * np.ones(n); up = np.ones(n - 1); lo = np.ones(n - 1)
    up[0] = 2.0          # row 0:  2 (u_1 - u_0)
    lo[-1] = 2.0         # row N:  2 (u_{N-1} - u_N)
    main[0] = -2.0; main[-1] = -2.0
    L = sp.diags([lo, main, up], [-1, 0, 1], format="csc") * (eps**2 / h**2)
    I = sp.identity(n, format="csc")

    def rhs(t, u):
        return L @ u - (u**3 - u)

    def jac(t, u):
        return L + sp.diags(1.0 - 3.0 * u**2, 0, format="csc")
    u = u0fun(x)
    t0 = 0.0
    res = {}
    for t1 in times:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solve_ivp(rhs, (t0, t1), u, method="Radau", rtol=rtol, atol=atol, jac=jac)
        assert sol.success
        u = sol.y[:, -1]
        res[t1] = u.copy()
        t0 = t1
    return x, res


times = [0.1, 0.5, 1.0]
out = {}
for ic in ("orig", "compat"):
    for eps in (0.1, 0.05, 0.01):
        sols = {}
        for Nx in (1000, 2000, 4000, 8000):
            t0 = time.time()
            x, r = fd_mol(eps, refsol._IC[ic], Nx, times)
            sols[Nx] = (x, r)
        rec = {}
        # trusted on the coarsest FD grid nodes
        xc = sols[1000][0]
        tr = refsol.reference(eps, ic, 1.0, times, xc, N=256)
        for t in times:
            e = {}
            for Nx in (1000, 2000, 4000, 8000):
                xx, r = sols[Nx]
                step = Nx // 1000
                e[Nx] = float(np.max(np.abs(r[t][::step] - tr[t])))
            # Richardson from 4000 and 8000 and from 2000,4000
            u4 = sols[4000][1][t][::4]; u8 = sols[8000][1][t][::8]; u2 = sols[2000][1][t][::2]
            rich48 = (4 * u8 - u4) / 3.0
            rich24 = (4 * u4 - u2) / 3.0
            e["Richardson(4000,8000)"] = float(np.max(np.abs(rich48 - tr[t])))
            e["Richardson(2000,4000)"] = float(np.max(np.abs(rich24 - tr[t])))
            rec[str(t)] = e
        out[f"{ic}_eps{eps}"] = rec
        print(ic, eps, {t: {k: f"{v:.2e}" for k, v in d.items()} for t, d in rec.items()}, flush=True)
json.dump(out, open(os.path.join(HERE, "independent_check.json"), "w"), indent=1)
