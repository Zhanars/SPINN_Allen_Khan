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
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
    ax[0].plot(t, Ea, "-", color="#3498db", lw=1.8, label=r"$\lambda_E=0$")
    ax[0].plot(t, Eb, "--", color="#e74c3c", lw=1.4, label=r"$\lambda_E=50$")
    ax[0].set_xlabel("t"); ax[0].set_ylabel("E(t)"); ax[0].legend(); ax[0].set_title("(a) energy E(t)")
    ax[1].semilogy(t, np.abs(Eb - Ea) + 1e-18, color="k", label=r"$|E_{50}-E_{0}|$")
    step = np.abs(np.diff(Ea)).min()
    ax[1].axhline(step, color="gray", ls=":", label="smallest decrement per step")
    ax[1].set_xlabel("t"); ax[1].legend(fontsize=6, loc="center right"); ax[1].set_title("(b) difference of the curves")
    # (c) penalty activity vs lambda_E over seeds: number of epochs with nonzero L_E and max lam*L_E
    lams = [0.0, 1.0, 10.0, 50.0, 500.0, 5000.0]
    pats = {0.0: "fh_K32_s*", 1.0: "ah_lamE1_s*", 10.0: "ah_lamE10_s*", 50.0: "fh_K32_E50_s*", 500.0: "ah_lamE500_s*", 5000.0: "ah_lamE5000_s*"}
    xs, ys, sds = [], [], []
    for lam in lams:
        rr = runs(pats[lam])
        if not rr:
            continue
        v = [r["errors"]["1.0"]["l2"] for r in rr]
        xs.append(max(lam, 0.3)); ys.append(np.mean(v)); sds.append(np.std(v, ddof=1) if len(v) > 1 else 0.0)
        ax[2].plot([max(lam, 0.3)] * len(v), v, "o", color="#7f8c8d", ms=2.5, alpha=0.6)
    ax[2].errorbar(xs, ys, yerr=sds, fmt="s-", color="#2c3e50", ms=4, capsize=2, label="mean $\\pm$ SD (5 seeds)")
    ax[2].set_xscale("log"); ax[2].set_yscale("log"); ax[2].set_ylim(3e-4, 3e-3)
    ax[2].set_xlabel(r"$\lambda_E$ (0 plotted at 0.3)"); ax[2].set_ylabel(r"$L^2$ error at $t=1$"); ax[2].set_title("(c) accuracy vs weight"); ax[2].legend(fontsize=6)
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_energy.png")); plt.close()


def fig_K():
    """Fig. 2: L2 error vs K, soft (as submitted) vs hard Neumann, individual seeds."""
    fig, ax = plt.subplots(1, 1, figsize=(3.8, 3.6))
    series = {"soft BC, equispaced grid (as submitted)": ("#e74c3c", "o", {32: "fwd_compat_K32_s*", 64: "fwd_compat_K64_s*", 128: "var_K128_uniform_soft"}),
              r"soft BC rescaled, $\lambda_{\rm BC}=K^4$": ("#8e44ad", "^", {32: "vs_K32_bcfix_s*", 64: "vs_K64_bcfix_s*", 128: "vs_K128_bcfix_s*"}),
              "hard Neumann basis, Gauss–Lobatto grid (revised)": ("#27ae60", "s", {32: "fh_K32_s*", 64: "fh_K64_s*", 128: "fh_K128_s*", 256: "kh_K256_x1"}),
              "revised, 2x training budget": ("#16a085", "D", {64: "kh_K64_x2", 128: "kh_K128_x2"})}
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
    ax.legend(fontsize=5, loc="upper center", bbox_to_anchor=(0.45, -0.16), ncol=2, frameon=False)
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_K.png")); plt.close()


def fig_balance():
    """Fig. 3: energy balance of the trained trajectory: E(0)-E(t) versus the a-posteriori bound."""
    r = run("fh_K32_E50_s42")
    s = r["balance"]["t_series"]
    t = np.array(s["t"]); n_r = np.array(s["n_r"]); n_ut2 = np.array(s["n_ut2"]); dEdt = np.array(s["dEdt"]); rhs = np.array(s["rhs"])
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.5))
    ax[0].plot(t, dEdt, "-", color="#2c3e50", label=r"$dE/dt=\int(\varepsilon^2u_xu_{xt}+f(u)u_t)\,dx$")
    ax[0].plot(t, rhs, "--", color="#e67e22", label=r"$-\|u_t\|^2+\langle r,u_t\rangle+\beta$")
    ax[0].plot(t, -n_ut2, ":", color="#7f8c8d", label=r"$-\|u_t\|^2$")
    ax[0].set_xlabel("t"); ax[0].legend(fontsize=6); ax[0].set_title("(a) identity, trained network")
    for tag, c, ls, lab in (("fh_K32_s42", "#3498db", "-", r"revised $K=32$, $\lambda_E=0$"), ("fh_K64_s42", "#27ae60", "-", r"revised $K=64$, $\lambda_E=0$"),
                            ("fh_K32_E50_s42", "#e74c3c", "--", r"revised $K=32$, $\lambda_E=50$")):
        try:
            q = run(tag)["balance"]["t_series"]
        except Exception:
            continue
        ax[1].semilogy(q["t"], np.array(q["n_r"]), ls, color=c, label=lab)
    ax[1].set_xlabel("t"); ax[1].set_ylabel(r"$\|r(\cdot,t)\|_{L^2}$"); ax[1].legend(fontsize=6); ax[1].set_title("(b) PDE residual of the trained networks")
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_balance.png")); plt.close()


def fig_inverse():
    """Fig. 4: (a) joint Adam: eps_hat/eps_init is nearly constant; (b) variable projection: error vs data weight (noise-free);
    (c) noisy data: variable projection versus classical least squares on the same noise realisations."""
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    allr = []
    for var, c, m in (("soft", "#e74c3c", "o"), ("hard", "#27ae60", "s")):
        for e, ls in ((0.05, "-"), (0.1, "--")):
            xs, ys = [], []
            for init in (0.5, 1.0, 1.5, 2.0, 4.0):
                f = os.path.join(RUNS, f"inv_adam_{var}_e{e}_i{init}.json")
                if os.path.exists(f):
                    d = json.load(open(f)); xs.append(init); ys.append(d["eps_hat"] / d["eps_true"]); allr.append(d["ratio_to_init"])
            if xs:
                ax[0].plot(xs, ys, m + ls, color=c, label=("as submitted" if var == "soft" else "revised") + r", $\varepsilon_{\rm true}$=" + f"{e}", ms=3.5, lw=1)
    ax[0].plot([0.5, 4], [0.5, 4], "k:", lw=1, label=r"$\hat\varepsilon=\varepsilon_{\rm init}$")
    if allr:
        q = float(np.mean(allr))
        ax[0].plot([0.5, 4], [q * 0.5, q * 4], "k-", lw=0.6, alpha=0.5, label=r"$\hat\varepsilon=%.2f\,\varepsilon_{\rm init}$" % q)
    ax[0].axhline(1.0, color="gray", lw=0.6, ls="--")
    ax[0].set_xlabel(r"$\varepsilon_{\rm init}/\varepsilon_{\rm true}$"); ax[0].set_ylabel(r"$\hat\varepsilon/\varepsilon_{\rm true}$"); ax[0].legend(fontsize=4.6)
    ax[0].set_title("(a) joint Adam")
    for e, c in ((0.05, "#2980b9"), (0.1, "#8e44ad")):
        xs, ys = [], []
        for ld in (20, 200, 2000, 20000):
            f = os.path.join(RUNS, f"inv_vp_e{e}_ld{ld}.json")
            if os.path.exists(f):
                d = json.load(open(f)); xs.append(ld); ys.append(100 * abs(d["eps_hat"] - d["eps_true"]) / d["eps_true"])
        if xs:
            ax[1].loglog(xs, ys, "o-", color=c, label=r"$\varepsilon_{\rm true}=%g$" % e, ms=3.5)
    ax[1].set_xlabel(r"$\lambda_{\rm data}$"); ax[1].set_ylabel("relative error (%)"); ax[1].legend(fontsize=6)
    ax[1].set_title("(b) variable projection, no noise")
    clf = os.path.join(HERE, "..", "data_r1", "classical2_inverse.json")
    cl = json.load(open(clf)) if os.path.exists(clf) else None
    xpos, labels = 0, []
    for e in (0.05, 0.1):
        for nz in (0.01, 0.05):
            nn = [os.path.join(RUNS, f"inv_vpnoise_e{e}_n{nz}_r{s}.json") for s in range(1, 6)]
            v = [100 * abs(json.load(open(f))["eps_hat"] - e) / e for f in nn if os.path.exists(f)]
            if v:
                ax[2].plot([xpos - 0.12] * len(v), v, "o", color="#27ae60", ms=3, alpha=0.8, label="Spectral-PINN, var. projection" if xpos == 0 else None)
            if cl:
                w = [100 * abs(cl[str(e)]["runs"][str(nz)]["eps_hat"][str(s)] - e) / e for s in range(1, 6)]
                ax[2].plot([xpos + 0.12] * len(w), w, "s", color="#2c3e50", ms=3, alpha=0.8, label="classical least squares" if xpos == 0 else None)
            labels.append(f"{e:g} / {nz:g}"); xpos += 1
    ax[2].set_yscale("log"); ax[2].set_xticks(range(len(labels))); ax[2].set_xticklabels(labels, fontsize=6)
    ax[2].set_xlabel(r"$\varepsilon_{\rm true}$ / $\sigma$"); ax[2].set_ylabel("relative error (%)"); ax[2].set_ylim(1e-2, 3e2); ax[2].legend(fontsize=5, loc="lower center")
    ax[2].set_title("(c) noisy data (paired)")
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_inverse.png")); plt.close()


def fig_b3():
    """Fig. 5: B3 (eps=0.05, T=5) profiles of the reference and of the trained networks at t = 1, 3, 5 (data_r1/profiles_b3.json)."""
    d = json.load(open(os.path.join(HERE, "..", "data_r1", "profiles_b3.json")))
    x = np.array(d["x"])
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4), sharey=True)
    lab = {"b3_K32_E0_s42": ("Spectral-PINN $K=32$", "#e67e22", "--"), "b3_K64_E0_s42": ("Spectral-PINN $K=64$", "#27ae60", "-."),
           "b3_vanilla_s42": ("vanilla PINN", "#8e44ad", ":")}
    for i, t in enumerate(("1.0", "3.0", "5.0")):
        ax[i].plot(x, d["ref"][t], "k-", lw=2.2, alpha=0.35, label="reference")
        for tag, (l, c, ls) in lab.items():
            if tag in d:
                ax[i].plot(x, d[tag][t], ls, color=c, lw=1.2, label=l)
        ax[i].set_title("t = %g" % float(t)); ax[i].set_xlabel("x")
    ax[0].set_ylabel("u(x,t)"); ax[0].legend(fontsize=5.5, loc="lower center")
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_b3.png")); plt.close()


if __name__ == "__main__":
    for f in (fig_energy, fig_K, fig_balance, fig_inverse, fig_b3):
        try:
            f(); print("ok", f.__name__)
        except Exception as ex:
            print("skip", f.__name__, repr(ex)[:120])


