"""Job lists for the round-1 revision."""
SEEDS = [42, 123, 7, 2024, 99]


def J(method, tag, **kw):
    d = {"method": method, "tag": tag}
    d.update(kw)
    return d


def smoke():
    return [J("spectral", "smoke_K16", K=16, epochs=200, warmup=100, lam_energy=50.0, seed=1)]


def forward():
    """Table 1: vanilla PINN, Spectral-PINN K=32, +energy penalty (lam_E=50), K=64; 5 seeds, compat IC (main) and orig IC."""
    jobs = []
    for s in SEEDS:                                   # long jobs first
        jobs.append(J("vanilla", f"fwd_compat_vanilla_s{s}", ic="compat", seed=s))
    for s in SEEDS:
        jobs.append(J("spectral", f"fwd_compat_K64_s{s}", ic="compat", K=64, seed=s))
    for s in SEEDS:
        jobs.append(J("spectral", f"fwd_compat_K32_s{s}", ic="compat", K=32, seed=s))
        jobs.append(J("spectral", f"fwd_compat_K32_E50_s{s}", ic="compat", K=32, seed=s, lam_energy=50.0))
    for s in SEEDS:
        jobs.append(J("spectral", f"fwd_orig_K32_s{s}", ic="orig", K=32, seed=s))
        jobs.append(J("spectral", f"fwd_orig_K32_E50_s{s}", ic="orig", K=32, seed=s, lam_energy=50.0))
        jobs.append(J("spectral", f"fwd_orig_K64_s{s}", ic="orig", K=64, seed=s))
    for s in SEEDS:
        jobs.append(J("vanilla", f"fwd_orig_vanilla_s{s}", ic="orig", seed=s))
    return jobs


def ablation_lamE():
    """lam_E ablation on the compat benchmark, K=32 (lam_E=0 and 50 come from the forward group), 5 seeds."""
    jobs = []
    for lam in (1.0, 10.0, 500.0, 5000.0):
        for s in SEEDS:
            jobs.append(J("spectral", f"abl_lamE{lam:g}_s{s}", ic="compat", K=32, seed=s, lam_energy=lam))
    return jobs


def penalty_regimes():
    """Regimes in which the unpenalised network is inaccurate: does the penalty change anything there?"""
    jobs = []
    for s in (42, 123, 7):
        jobs.append(J("spectral", f"reg_K64_E50_s{s}", ic="compat", K=64, seed=s, lam_energy=50.0))
        jobs.append(J("spectral", f"reg_short_E0_s{s}", ic="compat", K=32, seed=s, epochs=300, warmup=100))
        jobs.append(J("spectral", f"reg_short_E50_s{s}", ic="compat", K=32, seed=s, epochs=300, warmup=100, lam_energy=50.0))
        jobs.append(J("spectral", f"reg_nowarm_E0_s{s}", ic="compat", K=32, seed=s, warmup=0))
        jobs.append(J("spectral", f"reg_nowarm_E50_s{s}", ic="compat", K=32, seed=s, warmup=0, lam_energy=50.0))
    return jobs


def stress_K():
    jobs = [J("spectral", "k_K64_x2", ic="compat", K=64, seed=42, epochs=12000, warmup=2000),
            J("spectral", "k_K128_x1", ic="compat", K=128, seed=42),
            J("spectral", "k_K128_x4", ic="compat", K=128, seed=42, epochs=24000, warmup=4000),
            J("spectral", "k_K256_x1", ic="compat", K=256, seed=42)]
    return jobs


def stress_T():
    jobs = []
    for T in (5.0, 10.0):
        for lam in (0.0, 50.0):
            jobs.append(J("spectral", f"T{T:g}_E{lam:g}", ic="compat", K=32, seed=42, T=T, lam_energy=lam))
    return jobs


def variants():
    """Which implementation choice drives the K-degradation?  (grid: equispaced vs CGL; BC: soft penalty as submitted vs hard Neumann basis)"""
    jobs = []
    for K in (32, 64, 128):
        for grid in ("uniform", "cgl"):
            for hard in (False, True):
                if grid == "uniform" and not hard:
                    if K == 128:
                        jobs.append(J("spectral", f"var_K{K}_uniform_soft", ic="compat", K=K, seed=42))
                    continue                              # K=32/64 uniform-soft already in the forward group (seed 42)
                jobs.append(J("spectral", f"var_K{K}_{grid}_{'hard' if hard else 'soft'}", ic="compat", K=K, seed=42, grid=grid, hard_bc=hard))
    return jobs


HARD = dict(grid="cgl", hard_bc=True)          # revised Spectral-PINN: exact Neumann basis + Gauss-Lobatto collocation


def forward_hard():
    """Revised Spectral-PINN (hard Neumann basis, CGL collocation) on the compat benchmark, eps=0.01, T=1."""
    jobs = []
    for s in SEEDS:
        jobs.append(J("spectral", f"fh_K64_s{s}", ic="compat", K=64, seed=s, **HARD))
    for s in SEEDS:
        jobs.append(J("spectral", f"fh_K32_s{s}", ic="compat", K=32, seed=s, **HARD))
        jobs.append(J("spectral", f"fh_K32_E50_s{s}", ic="compat", K=32, seed=s, lam_energy=50.0, **HARD))
    for s in (42, 123, 7):
        jobs.append(J("spectral", f"fh_K128_s{s}", ic="compat", K=128, seed=s, **HARD))
    return jobs


def b2():
    """Diffusion-active benchmark eps=0.1, T=1 (smooth, no sharp interface): vanilla, as-submitted and revised Spectral-PINN."""
    jobs = []
    for s in (42, 123, 7):
        jobs.append(J("vanilla", f"b2_vanilla_s{s}", ic="compat", eps=0.1, seed=s))
    for s in SEEDS:
        jobs.append(J("spectral", f"b2_hard_K32_s{s}", ic="compat", eps=0.1, K=32, seed=s, **HARD))
    for s in (42, 123, 7):
        jobs.append(J("spectral", f"b2_soft_K32_s{s}", ic="compat", eps=0.1, K=32, seed=s))
    return jobs


def abl_hard():
    """lam_E ablation with the revised Spectral-PINN (lam_E = 0 and 50 come from forward_hard), K=32, eps=0.01."""
    jobs = []
    for lam in (1.0, 10.0, 500.0, 5000.0):
        for s in SEEDS:
            jobs.append(J("spectral", f"ah_lamE{lam:g}_s{s}", ic="compat", K=32, seed=s, lam_energy=lam, **HARD))
    return jobs


def regimes_hard():
    """Regimes with a poor trajectory: does the penalty then change the energy behaviour?"""
    jobs = []
    for s in (42, 123, 7):
        for lam in (0.0, 50.0):
            tag = f"E{lam:g}"
            kw = dict(lam_energy=lam) if lam > 0 else {}
            jobs.append(J("spectral", f"rh_short_{tag}_s{s}", ic="compat", K=32, seed=s, epochs=300, warmup=100, **HARD, **kw))
            jobs.append(J("spectral", f"rh_nowarm_{tag}_s{s}", ic="compat", K=32, seed=s, warmup=0, **HARD, **kw))
            jobs.append(J("spectral", f"rh_pde0_{tag}_s{s}", ic="compat", K=32, seed=s, epochs=1500, warmup=1000, lr=1e-4, **HARD, **kw))
    return jobs


def regimes_soft64():
    """The as-submitted K=64 network (soft penalty with effective weight K^-4: wall slopes 3.6-6.7, L2(t=1) ~ 1e-1) with the energy penalty;
    the same configuration without the penalty is fwd_compat_K64_s* (group forward)."""
    return [J("spectral", f"reg_K64_E50_s{s}", ic="compat", K=64, seed=s, lam_energy=50.0) for s in (42, 123, 7)]


def stress_K_hard():
    return [J("spectral", "kh_K64_x2", ic="compat", K=64, seed=42, epochs=12000, warmup=2000, **HARD),
            J("spectral", "kh_K128_x2", ic="compat", K=128, seed=42, epochs=12000, warmup=2000, **HARD),
            J("spectral", "kh_K256_x1", ic="compat", K=256, seed=42, **HARD)]


def queue1():
    return forward_hard() + b2() + abl_hard() + regimes_hard() + stress_K_hard()


def inverse():
    jobs = []
    def add(tag, **kw):
        d = {"method": "inverse", "tag": "inv_" + tag, "ic": "compat"}
        d.update(kw)
        jobs.append(d)
    soft = dict(grid="uniform", hard_bc=False)
    for e in (0.05, 0.1):
        # 1. joint Adam, as submitted: eps_init sweep, soft/uniform (as submitted) and hard/CGL (revised network)
        for init in (0.5, 1.0, 1.5, 2.0, 4.0):
            add(f"adam_soft_e{e}_i{init}", eps=e, mode="adam", eps_init_f=init, **soft)
            add(f"adam_hard_e{e}_i{init}", eps=e, mode="adam", eps_init_f=init)
        # 2. variable projection, noise free, data weight sweep
        for ld in (20, 200, 2000, 20000):
            add(f"vp_e{e}_ld{ld}", eps=e, mode="varpro", phaseA=3000, lam_data=ld, lam_dataA=1000)
        # 3. noise
        for nz in (0.01, 0.05):
            for ns in (1, 2, 3, 4, 5):
                add(f"vpnoise_e{e}_n{nz}_r{ns}", eps=e, mode="varpro", phaseA=3000, lam_data=2000, lam_dataA=1000, noise=nz, noise_seed=ns)
    # 4. the stiff, almost eps-independent case eps = 0.01
    for nz in (0.0, 0.01):
        add(f"vp_e0.01_n{nz}", eps=0.01, mode="varpro", phaseA=3000, lam_data=2000, lam_dataA=1000, noise=nz)
    # 5. observation design and loss weights (eps = 0.05, noise free, variable projection)
    for nx in (10, 160):
        add(f"vp_e0.05_nx{nx}", eps=0.05, mode="varpro", phaseA=3000, lam_data=2000, lam_dataA=1000, n_obs_x=nx)
    for nt, times in (("3", [0.4, 0.8, 1.0]), ("20", [round(0.05 * i, 2) for i in range(1, 21)])):
        add(f"vp_e0.05_nt{nt}", eps=0.05, mode="varpro", phaseA=3000, lam_data=2000, lam_dataA=1000, obs_times=times)
    for lp in (5.0, 500.0):
        add(f"vp_e0.05_lpde{lp:g}", eps=0.05, mode="varpro", phaseA=3000, lam_data=2000, lam_dataA=1000, lam_pde=lp)
    add("vp_e0.05_nophaseA", eps=0.05, mode="varpro", phaseA=0, lam_data=2000)
    add("vp_e0.05_K64", eps=0.05, mode="varpro", phaseA=3000, lam_data=2000, lam_dataA=1000, K=64)
    return jobs


def inverse_abl():
    """Joint Adam (the procedure of the first version) with the revised network, eps_true=0.05, eps_init=2 eps_true:
    which settings change the estimate?  (data weight, observation design, K, learning-rate multiplier of log eps, budget, PDE weight)"""
    jobs = []
    def add(tag, **kw):
        d = {"method": "inverse", "tag": "inv_ab_" + tag, "ic": "compat", "eps": 0.05, "mode": "adam"}
        d.update(kw)
        jobs.append(d)
    add("ld2000", lam_data=2000.0)
    add("ld2000_i1.0", lam_data=2000.0, eps_init_f=1.0)
    add("nx10", n_obs_x=10)
    add("nx160", n_obs_x=160)
    add("nt20", obs_times=[round(0.05 * i, 2) for i in range(1, 21)])
    add("K16", K=16)
    add("K64", K=64)
    add("mult1", eps_lr_mult=1.0)
    add("mult100", eps_lr_mult=100.0)
    add("ep16000", epochs=16000)
    add("lpde5", lam_pde=5.0)
    add("lpde500", lam_pde=500.0)
    return jobs


def b3():
    """Interface regime: eps=0.05, T=5 (K=64 representation floor 2e-4, K=32 floor 1e-2; the reaction-only ODE is 0.46 away)."""
    kw = dict(ic="compat", eps=0.05, T=5.0, epochs=8000, nt_col=250, n_ckpt=50)
    jobs = []
    for s in (42, 123):
        for lam in (0.0, 50.0):
            e = dict(lam_energy=lam) if lam > 0 else {}
            jobs.append(J("spectral", f"b3_K64_E{lam:g}_s{s}", K=64, seed=s, **kw, **HARD, **e))
    for lam in (0.0, 50.0):
        e = dict(lam_energy=lam) if lam > 0 else {}
        jobs.append(J("spectral", f"b3_K32_E{lam:g}_s42", K=32, seed=42, **kw, **HARD, **e))
    jobs.append(J("vanilla", "b3_vanilla_s42", ic="compat", eps=0.05, T=5.0, seed=42, epochs=8000, nx_col=60, nt_col=120))
    return jobs


def variants2():
    """Soft Neumann penalty with the weight rescaled so that its effective weight on u_x(+-1)^2 is 1 (lam_BC = K^4)."""
    jobs = []
    for K in (32, 64):
        for s in (42, 123, 7):
            jobs.append(J("spectral", f"vs_K{K}_bcfix_s{s}", ic="compat", K=K, seed=s, lam_bc=float(K) ** 4))
    jobs.append(J("spectral", "vs_K128_bcfix_s42", ic="compat", K=128, seed=42, lam_bc=128.0 ** 4))
    return jobs


def queue2():
    return variants2() + lbfgs() + b3() + [J("spectral", "var_K128_uniform_soft", ic="compat", K=128, seed=42)]


def lbfgs():
    """Stronger PINN baselines: Adam followed by 2000 full-batch L-BFGS iterations (B1 and B2)."""
    jobs = []
    for s in (42, 123, 7):
        jobs.append(J("vanilla", f"lb_vanilla_s{s}", ic="compat", seed=s, lbfgs=2000))
        jobs.append(J("spectral", f"lb_hard_K32_s{s}", ic="compat", K=32, seed=s, lbfgs=2000, **HARD))
    for s in (42, 123):
        jobs.append(J("spectral", f"lb_hard_K64_s{s}", ic="compat", K=64, seed=s, lbfgs=2000, **HARD))
    jobs.append(J("vanilla", "lb_b2_vanilla_s42", ic="compat", eps=0.1, seed=42, lbfgs=2000))
    jobs.append(J("spectral", "lb_b2_hard_K32_s42", ic="compat", eps=0.1, K=32, seed=42, lbfgs=2000, **HARD))
    return jobs


def _cost(job):
    """rough single-thread cost in minutes on an idle machine (used only to order the jobs, longest first)"""
    m = job["method"]
    if m == "vanilla":
        return 25.0 * job.get("epochs", 8000) / 8000 * (job.get("nx_col", 80) * job.get("nt_col", 80)) / 6400 * (1 + job.get("lbfgs", 0) / 2000)
    if m == "inverse":
        return 3.0 * (1 + job.get("phaseA", 0) / 4000) * (job.get("K", 32) / 32) ** 1.5
    T = job.get("T", 1.0); K = job.get("K", 32)
    ep = job.get("epochs", int(6000 + 1500 * (T - 1))) + job.get("warmup", 1000)
    nt = job.get("nt_col", int(100 * T))
    return 2.2 * (ep / 7000) * (nt / 100) * ((K + 1) / 33) ** 1.4 * (1 + job.get("lbfgs", 0) / 3000)


def all_jobs():
    """Everything that is still needed for the revision, without duplicates, longest jobs first.
    (The orig-IC jobs of the group 'forward' are trimmed to 3 seeds / no vanilla: they only serve Appendix D.)"""
    seen, jobs = set(), []
    for g in (forward_trimmed, forward_hard, b2, abl_hard, regimes_hard, stress_K_hard, variants2, lbfgs, b3, inverse, var3264, var128, inverse_abl, regimes_soft64):
        for jb in g():
            if jb["tag"] not in seen:
                seen.add(jb["tag"]); jobs.append(jb)
    return sorted(jobs, key=lambda jb: -_cost(jb))


def forward_trimmed():
    return [jb for jb in forward() if not (jb["tag"].startswith("fwd_orig_vanilla") or jb["tag"].endswith(("_s2024", "_s99")) and "fwd_orig" in jb["tag"] or ("fwd_orig" in jb["tag"] and "_E50_" in jb["tag"]))]


def var3264():
    """The K=32/64 implementation variants of Table 2 (the K=128 cells other than uniform-soft are not used)."""
    return [jb for jb in variants() if not jb["tag"].startswith("var_K128_")]


def var128():
    return [J("spectral", "var_K128_uniform_soft", ic="compat", K=128, seed=42)]


def timing():
    """Clean timing on an idle machine, one job at a time (run with nproc=1, after everything else has finished)."""
    return [J("vanilla", "timing_vanilla_s42", ic="compat", seed=42),
            J("spectral", "timing_soft_K32_s42", ic="compat", K=32, seed=42),
            J("spectral", "timing_hard_K32_s42", ic="compat", K=32, seed=42, **HARD),
            J("spectral", "timing_hard_K64_s42", ic="compat", K=64, seed=42, **HARD)]


JOBS = {"all": all_jobs, "timing": timing, "inverse_abl": inverse_abl, "regimes_soft64": regimes_soft64, "lbfgs": lbfgs, "queue2": queue2, "variants2": variants2, "b3": b3, "inverse": inverse, "forward_hard": forward_hard, "b2": b2, "abl_hard": abl_hard, "regimes_hard": regimes_hard, "stress_K_hard": stress_K_hard, "queue1": queue1, "variants": variants, "var3264": var3264, "smoke": smoke, "forward": forward, "ablation_lamE": ablation_lamE, "penalty_regimes": penalty_regimes,
        "stress_K": stress_K, "stress_T": stress_T}
