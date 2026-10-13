"""Quadrature check of the a-posteriori diagnostics (r1_core.balance_diagnostics).

The energy integrand and every term of the energy identity are polynomials of degree <= 4K in x and are integrated exactly by the
Clenshaw-Curtis rule of order 4K+8 (4K+9 nodes) used in r1_core.Quad.  ||r||^2 is a polynomial of degree 6K (r contains u^3), so that rule is not
exact for it.  This script recomputes int_0^T ||r||^2/2 dt for trained networks with the CC rule of order 6K+8 (6K+9 nodes, exact) and reports the
relative difference.   usage: python r1_check_quad.py  -> data_r1/check_quad.json"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import sys
import json
import glob
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import refsol
import r1_core
from common import f_ac, chebyshev_vandermonde, cheb_coeff_derivative_matrix
from spectral_pinn import SpectralPINN

torch.set_default_dtype(torch.float64)
RUNS = os.path.join(HERE, "..", "data_r1", "runs")
MODELS = os.path.join(HERE, "..", "data_r1", "models")


def int_half_r2(model, eps, T, Nq, nt=201):
    K = model.K
    x = np.cos(np.pi * np.arange(Nq + 1) / Nq)
    w = torch.tensor(refsol.clenshaw_curtis_weights(Nq))
    Tq = torch.tensor(chebyshev_vandermonde(x, K)).t()
    Dc = torch.tensor(cheb_coeff_derivative_matrix(K)); Dc2 = Dc @ Dc
    t = torch.linspace(1e-6, T, nt).reshape(-1, 1).requires_grad_(True)
    c = model.coeffs(t)
    cd = torch.cat([torch.autograd.grad(c[:, k], t, grad_outputs=torch.ones_like(c[:, k]), retain_graph=True)[0] for k in range(K + 1)], dim=1).detach()
    c = c.detach()
    u = c @ Tq; ut = cd @ Tq; uxx = (c @ Dc2.t()) @ Tq
    r = ut - eps ** 2 * uxx + f_ac(u)
    n_r2 = (r ** 2 * w).sum(-1).numpy()
    return float(np.trapezoid(0.5 * n_r2, t.detach().numpy().ravel()))


if __name__ == "__main__":
    tags = ["fh_K32_s42", "fh_K32_E50_s42", "fh_K64_s42", "fh_K128_s42", "fwd_compat_K32_s42", "fwd_compat_K64_s42", "b2_hard_K32_s42", "b3_K64_E0_s42", "vs_K64_bcfix_s42"]
    out = {}
    for tag in tags:
        fj, fm = os.path.join(RUNS, tag + ".json"), os.path.join(MODELS, tag + ".pt")
        if not (os.path.exists(fj) and os.path.exists(fm)):
            print("missing", tag); continue
        job = json.load(open(fj))["job"]
        K = job.get("K", 32); eps = job.get("eps", 0.01); T = float(job.get("T", 1.0))
        model = SpectralPINN(eps, K=K, seed=job.get("seed", 42), u0=r1_core.IC[job.get("ic", "compat")], grid=job.get("grid", "uniform"), hard_bc=job.get("hard_bc", False))
        model.net.load_state_dict(torch.load(fm))
        a = int_half_r2(model, eps, T, 4 * K + 8)
        b = int_half_r2(model, eps, T, 6 * K + 8)
        stored = json.load(open(fj))["balance"]["int_half_r2"]
        out[tag] = {"cc_4K8": a, "cc_6K8_exact": b, "rel_diff": abs(a - b) / b, "stored": stored}
        print(tag, out[tag], flush=True)
    json.dump(out, open(os.path.join(HERE, "..", "data_r1", "check_quad.json"), "w"), indent=1)
