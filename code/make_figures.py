"""Generate the energy-dissipation figure from saved experiment JSON."""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIG_DIR, exist_ok=True)


def energy_dissipation_figure():
    spec = json.load(open(os.path.join(DATA_DIR, "spectral_pinn.json")))
    spece = json.load(open(os.path.join(DATA_DIR, "spectral_pinn_energy.json")))

    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.plot(spec["t_grid"], spec["E_traj"], "o-", ms=3, lw=1.6, color="#3498db",
             label="Spectral-PINN (no energy penalty)")
    ax.plot(spece["t_grid"], spece["E_traj"], "s--", ms=3, lw=1.6, color="#e74c3c",
             label="Spectral-PINN + energy-constrained")
    ax.set_xlabel("t")
    ax.set_ylabel(r"Ginzburg--Landau energy $E(t)$")
    ax.set_title(r"Energy dissipation along the Spectral-PINN trajectory" "\n" r"($\varepsilon=0.01$)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    outpath = os.path.join(FIG_DIR, "energy_dissipation.png")
    plt.savefig(outpath, dpi=200)
    print("saved", outpath)


if __name__ == "__main__":
    energy_dissipation_figure()
