"""Extra data-derived facts relevant to the numbers audit (read-only on code/ and data/).
 1. IC compatibility with Neumann BC (u0'(+-1)) and what the trained K=32 snapshots do at the boundary
 2. Maximum-bound-principle check |u|<=1 on all stored snapshots (single-run seed 42)
 3. Trainable-parameter counts behind the 'DoF' column of Table 1
 4. Accuracy of the coefficient-space derivative matrix D_c for K=20 (Section 2.3 claim 1e-11..1e-13)
 5. Ginzburg-Landau energy of the FD reference vs the network energies quoted in the data files
Run:  python extra_checks.py > extra_checks_output.txt
"""
import sys, os, json
sys.dont_write_bytecode = True
CODE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))
DATA = os.path.abspath(os.path.join(CODE, "..", "data"))
sys.path.insert(0, CODE)
import numpy as np
from common import (reference_solution_fd, cheb_coeff_derivative_matrix, chebyshev_vandermonde,
                    cheb_boundary_deriv_vectors, fit_chebyshev_coeffs, u0_fn, F_ac)

J = lambda n: json.load(open(os.path.join(DATA, n)))
xt = np.linspace(-1, 1, 200)

print("=== 1. IC compatibility with u_x(+-1)=0 ===")
du0 = lambda x: 2 * x * np.cos(np.pi * x) - np.pi * x ** 2 * np.sin(np.pi * x)
print("u0'(-1) = %.6f, u0'(+1) = %.6f  (analytic: +2 and -2)" % (du0(-1.0), du0(1.0)))
K = 32
vp, vm = cheb_boundary_deriv_vectors(K)
xg = np.linspace(-1, 1, 300)
c0 = fit_chebyshev_coeffs(u0_fn(xg), xg, K)
print("K=32 projection of u0: max |u0 - sum c0_k T_k| on 300-pt grid = %.2e ; u_x(+1)=%.6f u_x(-1)=%.6f" % (
    np.max(np.abs(u0_fn(xg) - chebyshev_vandermonde(xg, K) @ c0)), c0 @ vp, c0 @ vm))
bc_scale = max(1.0, K ** 2)
print("code normalises boundary vectors by K^2=%d => BC residual of the IC = %.3e, squared %.3e (relative to O(1) data)" % (bc_scale, (c0 @ vp) / bc_scale, ((c0 @ vp) / bc_scale) ** 2))
print("Trained snapshots (seed-42 single runs): Chebyshev refit (K=32, 200 points) -> boundary derivatives")
for nm in ("spectral_pinn.json", "spectral_pinn_energy.json"):
    d = J(nm)
    for t in ("0.2", "0.5", "0.8", "1.0"):
        up = np.array(d["errors"][t]["u_pred"])
        c = np.linalg.lstsq(chebyshev_vandermonde(xt, K), up, rcond=None)[0]
        res = np.max(np.abs(chebyshev_vandermonde(xt, K) @ c - up))
        print("  %-24s t=%s u_x(+1)=%+.4f u_x(-1)=%+.4f (refit residual %.1e)" % (nm, t, c @ vp, c @ vm, res))

print("\n=== 2. Maximum bound principle |u|<=1 on stored snapshots ===")
for nm in ("vanilla_pinn.json", "spectral_pinn.json", "spectral_pinn_energy.json"):
    d = J(nm)
    for t in ("0.2", "0.5", "0.8", "1.0"):
        up = np.array(d["errors"][t]["u_pred"])
        print("  %-24s t=%s max|u|=%.5f  (excess over 1: %+.2e) min u=%.5f max u=%.5f" % (nm, t, np.max(np.abs(up)), np.max(np.abs(up)) - 1, up.min(), up.max()))
x_ref, snaps = reference_solution_fd(0.01)
print("  FD reference (Nx=500): max|u| over saved times t=0..1: %.6f ; max|u0| = %.6f" % (max(np.max(np.abs(v)) for v in snaps.values()), np.max(np.abs(u0_fn(x_ref)))))
st = J("stress_time_sweep.json")["results"]
print("  T-sweep: Linf errors at t>=5 within the T=10 run are 1.0-1.05, i.e. the network is far from the reference; no snapshots saved, so |u|<=1 cannot be checked for T=5,10 (UNVERIFIABLE).")
print("  Forward multi-seed / K / eps sweeps: only error norms saved, no snapshots => MBP unverifiable there.")

print("\n=== 3. Trainable parameter counts (Table 1 'DoF' column) ===")
mlp = lambda dims: sum(a * b + b for a, b in zip(dims[:-1], dims[1:]))
print("vanilla [2,64,64,64,1]: %d (manuscript 8,577)" % mlp([2, 64, 64, 64, 1]))
for K_ in (32, 64, 128, 256):
    print("CoeffNet [1,64,64,%d]: %d trainable weights (K+1=%d coefficients)" % (K_ + 1, mlp([1, 64, 64, K_ + 1]), K_ + 1))
print("Fourier 2D M=8 [1,64,64,81]: %d ; M=16 [1,64,64,289]: %d ; vanilla 2D [5,64,64,64,1]: %d" % (mlp([1, 64, 64, 81]), mlp([1, 64, 64, 289]), mlp([5, 64, 64, 64, 1])))
print("Chebyshev N=80 unknowns: 81")

print("\n=== 4. D_c accuracy (K=20) ===")
K = 20
Dc = cheb_coeff_derivative_matrix(K)
xn = np.cos(np.pi * np.arange(K + 1) / K)
xf = np.linspace(-1, 1, 2001)
for nm, f, df in [("exp(x)", np.exp, np.exp), ("sin(2x)", lambda x: np.sin(2 * x), lambda x: 2 * np.cos(2 * x)),
                  ("sin(3x)", lambda x: np.sin(3 * x), lambda x: 3 * np.cos(3 * x)),
                  ("cos(pi x)", lambda x: np.cos(np.pi * x), lambda x: -np.pi * np.sin(np.pi * x)),
                  ("x^2 cos(pi x) [the IC]", lambda x: x ** 2 * np.cos(np.pi * x), du0),
                  ("tanh(5x)", lambda x: np.tanh(5 * x), lambda x: 5 / np.cosh(5 * x) ** 2)]:
    c = np.linalg.lstsq(chebyshev_vandermonde(xn, K), f(xn), rcond=None)[0]
    b = Dc @ c
    err = np.max(np.abs(chebyshev_vandermonde(xf, K) @ b - df(xf)))
    print("  %-24s max |D_c-derivative - exact| = %.2e" % (nm, err))

print("\n=== 5. Ginzburg-Landau energy of the FD reference (eps=0.01) ===")
eps = 0.01
def GL(x, u):
    ux = np.gradient(u, x)
    integ = eps ** 2 / 2 * ux ** 2 + F_ac(u)
    return float(np.sum(0.5 * (integ[1:] + integ[:-1]) * np.diff(x)))
print("E[u0] (FD grid) = %.6f ; stored network E(0) = %.6f (spectral_pinn.json) / %.6f (energy)" % (GL(x_ref, u0_fn(x_ref)), J("spectral_pinn.json")["E_traj"][0], J("spectral_pinn_energy.json")["E_traj"][0]))
print("E[u_ref(t=1)] = %.6f ; stored network E(1) = %.6f / %.6f" % (GL(x_ref, snaps[1.0]), J("spectral_pinn.json")["E_traj"][-1], J("spectral_pinn_energy.json")["E_traj"][-1]))
Es = [GL(x_ref, snaps[t]) for t in sorted(snaps)]
print("reference energy monotone non-increasing on t=0,0.1,..,1: %s ; values %s" % (all(np.diff(Es) <= 1e-12), ["%.4f" % e for e in Es]))
x10, s10 = reference_solution_fd(eps, T=10.0, save_times=[1.0, 5.0, 10.0])
print("reference E(t=1,5,10) = %.4f %.4f %.4f  vs network (T=5 run) E_last=%.4f (t=5), (T=10 run) E_last=%.4f (t=10)" % (
    GL(x10, s10[1.0]), GL(x10, s10[5.0]), GL(x10, s10[10.0]), st[1]["E_last"], st[2]["E_last"]))

print("\n=== 6. Growth of the Chebyshev second-derivative operator (manuscript L209: 'O(K^4) growth of the condition number') ===")
from common import cheb_diff_matrix
Ks = [32, 64, 128, 256]
n_c, n_n = [], []
for K_ in Ks:
    Dcc = cheb_coeff_derivative_matrix(K_)
    n_c.append(np.linalg.norm(Dcc @ Dcc, 2))
    Dn, _ = cheb_diff_matrix(K_)
    n_n.append(np.max(np.abs(np.linalg.eigvals((Dn @ Dn)[1:-1, 1:-1]))))
    print("  K=%3d  ||D_c^2||_2 = %.4e   spectral radius of interior nodal D^2 = %.4e" % (K_, n_c[-1], n_n[-1]))
print("  log-log slope of ||D_c^2||_2 vs K (32..256): %.2f ; slope of nodal D^2 spectral radius: %.2f" % (np.polyfit(np.log(Ks), np.log(n_c), 1)[0], np.polyfit(np.log(Ks), np.log(n_n), 1)[0]))
print("  (D_c^2 is nilpotent, so its condition number is undefined; the O(K^4) growth is a property of the norm / of the nodal operator)")
