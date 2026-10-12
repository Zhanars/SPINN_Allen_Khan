"""Where does the Chebyshev-vs-FD(500,1e-4) discrepancy sit? (supports the explanation of the norm mismatch factor)"""
import sys, os
sys.dont_write_bytecode = True
CODE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))
sys.path.insert(0, CODE)
import numpy as np
from common import reference_solution_fd
from chebyshev_spectral import solve_chebyshev, clenshaw_curtis_weights
x_ref, snaps = reference_solution_fd(0.01)
xf, sf = reference_solution_fd(0.01, Nx=4000, dt_ref=1e-5, save_times=[1.0])
for N, dt in [(40, 2e-3), (80, 5e-4)]:
    x, snp = solve_chebyshev(0.01, N=N, dt=dt)
    u = snp[1.0]
    e = u - np.interp(x[::-1], x_ref, snaps[1.0])[::-1]
    w = clenshaw_curtis_weights(N)
    j = np.argmax(np.abs(e))
    print("N=%d: max|e| = %.3e at x=%.4f ; share of CC-weighted e^2 within |x|>0.9: %.3f ; within |x|>0.5: %.3f ; share of node count with |x|>0.9: %.3f" % (
        N, np.abs(e[j]), x[j], np.sum(w * e ** 2 * (np.abs(x) > 0.9)) / np.sum(w * e ** 2), np.sum(w * e ** 2 * (np.abs(x) > 0.5)) / np.sum(w * e ** 2), np.mean(np.abs(x) > 0.9)))
    # the FD reference error vs fine FD, where is it?
    d = np.interp(x_ref, xf, sf[1.0]) - snaps[1.0]
    k = np.argmax(np.abs(d))
    print("   FD(500) minus FD(4000): max|d| = %.3e at x=%.3f" % (np.abs(d[k]), x_ref[k]))
