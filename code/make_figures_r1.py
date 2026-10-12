"""Figures of the revised manuscript from data_r1 (PNG, 300 dpi, MDPI-width)."""
import os
import sys
import json
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r1_tables as T

RUNS = T.RUNS
FIG = os.path.join(HERE, "..", "manuscript", "figures_r1")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 150, "savefig.dpi": 300})


def run(tag):
    return json.load(open(os.path.join(RUNS, tag + ".json")))


def runs(pat):
    return [json.load(open(f)) for f in sorted(glob.glob(os.path.join(RUNS, pat + ".json")))]


def fig_energy():
    """Fig. 1: energy trajectories with and without the penalty, and their difference; penalty activity."""
    a, b = run("fh_K32_s42"), run("fh_K32_E50_s42")
    t = np.linspace(0, 1, len(a["E_traj_dt0.02"]))
    Ea, Eb = np.array(a["E_traj_dt0.02"]), np.array(b["E_traj_dt0.02"])
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    ax[0].plot(t, Ea, "-", color="#3498db", lw=1.8, label=r"$\lambda_E=0$")
    ax[0].plot(t, Eb, "--", color="#e74c3c", lw=1.4, label=r"$\lambda_E=50$")
    ax[0].set_xlabel("t"); ax[0].set_ylabel("E(t)"); ax[0].legend(); ax[0].set_title("(a) energy E(t)")
    ax[1].semilogy(t, np.abs(Eb - Ea) + 1e-18, color="k")
    step = np.abs(np.diff(Ea)).min()
    ax[1].axhline(step, color="gray", ls=":", label="smallest decrement per step")
    ax[1].set_xlabel("t"); ax[1].set_ylabel(r"$|E_{50}-E_{0}|$"); ax[1].legend(fontsize=6); ax[1].set_title("(b) difference")
    # (c) penalty activity vs lambda_E over seeds: number of epochs with nonzero L_E and max lam*L_E
    lams = [0.0, 1.0, 10.0, 50.0, 500.0, 5000.0]
    pats = {0.0: "fh_K32_s*", 1.0: "ah_lamE1_s*", 10.0: "ah_lamE10_s*", 50.0: "fh_K32_E50_s*", 500.0: "ah_lamE500_s*", 5000.0: "ah_lamE5000_s*"}
    xs, ys = [], []
    for lam in lams:
        rr = runs(pats[lam])
        if not rr:
            continue
        xs.append(max(lam, 0.3)); ys.append(np.mean([r["errors"]["1.0"]["l2"] for r in rr]))
    ax[2].semilogx(xs, ys, "o-", color="#2c3e50")
    ax[2].set_xlabel(r"$\lambda_E$ (0 plotted at 0.3)"); ax[2].set_ylabel(r"mean $L^2$ error at $t=1$"); ax[2].set_title("(c) error vs weight")
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_energy.png")); plt.close()


def fig_K():
    """Fig. 2: L2 error vs K, soft (as submitted) vs hard Neumann, individual seeds."""
    fig, ax = plt.subplots(1, 1, figsize=(3.4, 2.7))
    series = {"soft BC, equispaced grid (as submitted)": ("#e74c3c", "o", {32: "fwd_compat_K32_s*", 64: "fwd_compat_K64_s*", 128: "var_K128_uniform_soft"}),
              "hard Neumann basis, Gauss--Lobatto grid": ("#27ae60", "s", {32: "fh_K32_s*", 64: "fh_K64_s*", 128: "fh_K128_s*", 256: "kh_K256_x1"})}
    for lab, (c, m, d) in series.items():
        xs, mean = [], []
        for K, pat in sorted(d.items()):
            rr = runs(pat)
            if not rr:
                continue
            v = [r["errors"]["1.0"]["l2"] for r in rr]
            ax.plot([K] * len(v), v, m, color=c, alpha=0.35, ms=3)
            xs.append(K); mean.append(np.mean(v))
        if xs:
            ax.plot(xs, mean, m + "-", color=c, label=lab, ms=5)
    cl = json.load(open(os.path.join(HERE, "..", "data_r1", "classical.json")))["compat_eps0.01"]
    ax.axhline(cl["ode"]["1.0"]["l2"], color="gray", ls="--", lw=1, label=r"reaction-only ODE ($\varepsilon=0$)")
    ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.set_xlabel("K"); ax.set_ylabel(r"$L^2$ error at $t=1$ ($\varepsilon=0.01$)")
    ax.legend(fontsize=5.5, loc="upper left")
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_K.png")); plt.close()


def fig_balance():
    """Fig. 3: energy balance of the trained trajectory: E(0)-E(t) versus the a-posteriori bound."""
    r = run("fh_K32_E50_s42")
    s = r["balance"]["t_series"]
    t = np.array(s["t"]); n_r = np.array(s["n_r"]); n_ut2 = np.array(s["n_ut2"]); dEdt = np.array(s["dEdt"]); rhs = np.array(s["rhs"])
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.5))
    ax[0].plot(t, dEdt, "-", color="#2c3e50", label=r"$dE/dt$ (autodiff)")
    ax[0].plot(t, rhs, "--", color="#e67e22", label=r"$-\|u_t\|^2+\langle r,u_t\rangle+\beta$")
    ax[0].plot(t, -n_ut2, ":", color="#7f8c8d", label=r"$-\|u_t\|^2$")
    ax[0].set_xlabel("t"); ax[0].legend(fontsize=6); ax[0].set_title("(a) identity, trained network")
    for tag, c in (("fh_K32_s42", "#3498db"), ("fh_K64_s42", "#27ae60"), ("fh_K32_E50_s42", "#e74c3c")):
        try:
            q = run(tag)["balance"]["t_series"]
        except Exception:
            continue
        ax[1].semilogy(q["t"], np.array(q["n_r"]), color=c, label=tag)
    ax[1].set_xlabel("t"); ax[1].set_ylabel(r"$\|r(\cdot,t)\|_{L^2}$"); ax[1].legend(fontsize=6); ax[1].set_title("(b) PDE residual of the trained networks")
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_balance.png")); plt.close()


def fig_inverse():
    """Fig. 4: (a) joint Adam: eps_hat/eps_init is constant; (b) variable projection: error vs data weight."""
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for var, c, m in (("soft", "#e74c3c", "o"), ("hard", "#27ae60", "s")):
        for e, ls in ((0.05, "-"), (0.1, "--")):
            xs, ys = [], []
            for init in (0.5, 1.0, 1.5, 2.0, 4.0):
                f = os.path.join(RUNS, f"inv_adam_{var}_e{e}_i{init}.json")
                if os.path.exists(f):
                    d = json.load(open(f)); xs.append(init); ys.append(d["eps_hat"] / d["eps_true"])
            if xs:
                ax[0].plot(xs, ys, m + ls, color=c, label=f"{var}, eps={e}", ms=4)
    ax[0].plot([0.5, 4], [0.5, 4], "k:", lw=1, label=r"$\hat\varepsilon=\varepsilon_{\rm init}$")
    ax[0].plot([0.5, 4], [0.717 * 0.5, 0.717 * 4], "k-", lw=0.6, alpha=0.5, label=r"$0.717\,\varepsilon_{\rm init}$")
    ax[0].set_xlabel(r"$\varepsilon_{\rm init}/\varepsilon_{\rm true}$"); ax[0].set_ylabel(r"$\hat\varepsilon/\varepsilon_{\rm true}$"); ax[0].legend(fontsize=5.5)
    ax[0].set_title("(a) joint Adam: the estimate follows the start")
    for e, c in ((0.05, "#2980b9"), (0.1, "#8e44ad")):
        xs, ys = [], []
        for ld in (20, 200, 2000, 20000):
            f = os.path.join(RUNS, f"inv_vp_e{e}_ld{ld}.json")
            if os.path.exists(f):
                d = json.load(open(f)); xs.append(ld); ys.append(100 * abs(d["eps_hat"] - d["eps_true"]) / d["eps_true"])
        if xs:
            ax[1].loglog(xs, ys, "o-", color=c, label=rf"$\varepsilon={e}$", ms=4)
    ax[1].set_xlabel(r"$\lambda_{\rm data}$"); ax[1].set_ylabel("relative error of the estimate (%)"); ax[1].legend(fontsize=6)
    ax[1].set_title("(b) variable projection, noise-free data")
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_inverse.png")); plt.close()


def fig_b3():
    """Fig. B3: profiles at t = 5 (eps = 0.05): reference, vanilla PINN, revised Spectral-PINN K=32 and K=64 (needs data_r1/models)."""
    import torch
    from vanilla_pinn import MLP
    import spectral_pinn as sp
    import r1_core
    tags = {"K=64": "b3_K64_E0_s42", "K=32": "b3_K32_E0_s42", "vanilla": "b3_vanilla_s42"}
    ref = r1_core.reference(0.05, "compat", [1.0, 2.0, 3.0, 4.0, 5.0])
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    ax[0].plot(r1_core.XE, ref[5.0], "k-", lw=2, label="reference")
    for lab, tag in tags.items():
        f = os.path.join(HERE, "..", "data_r1", "models", tag + ".pt")
        if not os.path.exists(f):
            continue
        if lab == "vanilla":
            net = MLP([2, 64, 64, 64, 1]); net.load_state_dict(torch.load(f)); 
            from vanilla_pinn import eval_on_grid
            u = eval_on_grid(net, r1_core.XE, 5.0)
        else:
            K = int(lab.split("=")[1])
            m = sp.SpectralPINN(0.05, K=K, u0=r1_core.IC["compat"], grid="cgl", hard_bc=True)
            m.net.load_state_dict(torch.load(f))
            u = sp.eval_on_grid(m, 5.0, r1_core.XE)
        ax[0].plot(r1_core.XE, u, "--", lw=1.2, label=lab)
        ax[1].semilogy(r1_core.XE, np.abs(u - ref[5.0]) + 1e-12, lw=1.0, label=lab)
    ax[0].set_xlabel("x"); ax[0].set_ylabel("u(x, 5)"); ax[0].legend(fontsize=6); ax[0].set_title("(a) profile at t = 5")
    ax[1].set_xlabel("x"); ax[1].set_ylabel("|u - u_ref|"); ax[1].legend(fontsize=6); ax[1].set_title("(b) pointwise error")
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_b3.png")); plt.close()


if __name__ == "__main__":
    for f in (fig_energy, fig_K, fig_balance, fig_inverse, fig_b3):
        try:
            f(); print("ok", f.__name__)
        except Exception as ex:
            print("skip", f.__name__, repr(ex)[:120])
