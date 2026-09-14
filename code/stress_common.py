"""Shared helpers for the large-scale stress tests (K, eps, T, 2D axes)."""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import numpy as np


def reference_nx_for_eps(eps, min_nx=500, target_dx_over_eps=0.3):
    """Scale the FD reference grid so dx << eps (resolve the interface)."""
    nx = max(min_nx, int(np.ceil((2.0 / (target_dx_over_eps * eps)) / 2) * 2))
    return min(nx, 20000)


def chebyshev_n_for_eps(eps, min_n=80, target_dx_over_eps=0.5):
    n = max(min_n, int(np.ceil(2.0 / (target_dx_over_eps * eps))))
    return min(n, 400)
