"""Evaluate trained networks on the 200-point grid at selected times (for the profile figure of B3).
usage: python r1_profiles.py  -> data_r1/profiles_b3.json   (needs data_r1/models/*.pt)"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import sys
import json
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r1_core
from common import chebyshev_vandermonde
from spectral_pinn import SpectralPINN
from vanilla_pinn import MLP, eval_on_grid

torch.set_default_dtype(torch.float64)
RUNS = os.path.join(HERE, "..", "data_r1", "runs")
MODELS = os.path.join(HERE, "..", "data_r1", "models")
TIMES = [1.0, 3.0, 5.0]


def load_spectral(tag):
    job = json.load(open(os.path.join(RUNS, tag + ".json")))["job"]
    K = job.get("K", 32)
    model = SpectralPINN(job.get("eps", 0.01), K=K, seed=job.get("seed", 42), u0=r1_core.IC[job.get("ic", "compat")],
                         grid=job.get("grid", "uniform"), hard_bc=job.get("hard_bc", False))
    model.net.load_state_dict(torch.load(os.path.join(MODELS, tag + ".pt")))
    Tm = torch.tensor(chebyshev_vandermonde(r1_core.XE, K))
    def u(t):
        with torch.no_grad():
            return (model.coeffs(torch.tensor([[float(t)]])) @ Tm.t()).numpy().ravel()
    return u


def load_vanilla(tag):
    net = MLP([2, 64, 64, 64, 1])
    net.load_state_dict(torch.load(os.path.join(MODELS, tag + ".pt")))
    return lambda t: eval_on_grid(net, r1_core.XE, float(t))


if __name__ == "__main__":
    ref = r1_core.reference(0.05, "compat", [0.2, 0.5, 0.8, 1.0, 2.0, 3.0, 5.0])
    out = {"x": r1_core.XE.tolist(), "ref": {str(t): ref[t].tolist() for t in TIMES}}
    for tag, loader in (("b3_K32_E0_s42", load_spectral), ("b3_K64_E0_s42", load_spectral), ("b3_vanilla_s42", load_vanilla)):
        if not os.path.exists(os.path.join(MODELS, tag + ".pt")):
            print("missing", tag); continue
        f = loader(tag)
        out[tag] = {str(t): f(t).tolist() for t in TIMES}
        print(tag, {t: float(np.sqrt(np.sum((np.array(out[tag][str(t)]) - ref[t]) ** 2) * r1_core.DXE)) for t in TIMES}, flush=True)
    json.dump(out, open(os.path.join(HERE, "..", "data_r1", "profiles_b3.json"), "w"))
