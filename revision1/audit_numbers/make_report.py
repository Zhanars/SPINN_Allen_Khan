# -*- coding: utf-8 -*-
"""Generates audit_numbers.md from the row list below (so the status counts in the summary cannot drift from the table).
All numbers quoted in the rows come from compute_all_output.txt, norm_check_output.txt, extra_checks_output.txt."""
import os, collections

HERE = os.path.dirname(os.path.abspath(__file__))

# (id, location, claim as written (exact LaTeX), recomputed from data/code, status, proposed corrected wording)
R = []
def add(i, loc, claim, val, st, fix):
    assert st in ("OK", "WRONG", "UNVERIFIABLE"), st
    R.append((i, loc, claim, val, st, fix))

# ============================================================ ABSTRACT
add("ABS-1", "Abstract, L53",
    "`cuts the $L^2$ error by a factor of $\\sim$49 relative to a vanilla PINN`",
    "Ratio of 5-seed means 6.7377e-2 / 1.3785e-3 = **48.88** (energy-constrained variant 48.91). Bootstrap 95% CI [45.0, 52.6]; delta-method 48.9 +- 2.2 (1 SE). Same-seed paired ratios 46.5, 50.2, 54.7, 50.9, 42.5 (mean 49.0, median 50.2).",
    "OK",
    "\"by a factor of about 49 (ratio of the 5-seed mean L2 errors, 6.74e-2 vs 1.38e-3; 95% bootstrap interval 45-53)\". State the statistic in the abstract or in Sec. 3.1.")
add("ABS-2", "Abstract, L53",
    "`while training $\\sim$5.4$\\times$ faster`",
    "Ratio of mean elapsed times 919.33 s / 169.42 s = **5.43** (plain Spectral-PINN only). Energy-constrained variant 191.13 s -> **4.81**. Per-seed ratios 5.83, 5.90, 5.10, 4.79, 5.70. Independent eps-sweep batch (seed 42): 605.1 / 112.0 = 5.40. Single-run files (seed 42): 787.4 / 178.9 = 4.40. `elapsed` excludes the 1000-epoch warm-up (t0 set after it); vanilla runs 8000 epochs, Spectral-PINN 6000.",
    "OK",
    "\"about 5.4x less training time per run (4.8x for the energy-constrained variant; 8000 vs 6000 epochs; warm-up not timed)\". The same sentence in the abstract pairs 5.4x with the energy-constrained result, for which the factor is 4.8x.")
add("ABS-3", "Abstract, L53",
    "`with zero violations under an added energy-constrained loss`",
    "spectral_pinn_energy.json: n_energy_violations = 0, len(E_traj) = 51 (t = 0, 0.02, ..., 1), tolerance 1e-10, seed 42 only. No energy-violation count is stored for the 5-seed runs (forward_multiseed.json has only l2/linf/elapsed).",
    "OK",
    "Count is correct but covers one seed; say \"0 of 51 test times (seed 42)\" or rerun the energy check inside the multi-seed script.")
add("ABS-3b", "Abstract, L53 (implied attribution)",
    "same sentence: zero violations are presented as a result of the energy-constrained loss",
    "Control WITHOUT penalty (spectral_pinn.json): n_energy_violations = 0 of 51 as well; every consecutive step decreases by at least 1.285e-3 (tolerance 1e-10). Penalty term barely active: L_E nonzero at 1 of 120 logged epochs, max 7.56e-8 (lam_E * L_E <= 3.8e-6 vs total loss 4.1e-5). Energy curves with/without penalty differ by <= 4.5e-7 (relative 1.3e-6).",
    "WRONG",
    "\"The Spectral-PINN trajectory is energy-non-increasing at all 51 test times both with and without the penalty (seed 42); the penalty was inactive in this regime.\" Do not claim the penalty produces the zero count.")
add("ABS-4", "Abstract, L53",
    "`A stress test up to $K=256$ modes, however, shows accuracy \\emph{degrading}, not improving, even under proportionally larger training budgets`",
    "Matched-budget L2 (stress_k_sweep.json): K=64 2.25e-2, K=128 2.19e-1, K=256 1.47 (K=32: 1.38e-3). Same-budget (1x): 6.18e-2, 3.36e-1, 1.51.",
    "OK", "-")
add("ABS-5", "Abstract, L53",
    "`the bias stays flat across a $6\\times$ range of $K$`",
    "Column means 33.18, 40.05, 37.23, 40.91 % for K = 16, 32, 64, 96 (K=16 lowered by one -12.6 % run; without it 38.32). 16 -> 96 is exactly 6.0x. One seed (42) per cell.",
    "OK", "Add \"(single seed per configuration)\".")
add("ABS-6", "Abstract L53; Sec. 3.2 L184; Limitations L306; Conclusions L311",
    "`$\\sim$38\\% recovery bias` / `robust $\\sim38\\%$ overestimation bias`",
    "Mean absolute relative error over the 20 configurations 37.84 %; signed mean 36.58 %; 19 of 20 are overestimates, 1 (eps=0.01, K=16) is an UNDER-estimate (eps_hat = 0.008739).",
    "OK", "\"about 38 % (mean absolute relative error, 20 single-seed runs; 19 over- and 1 under-estimate)\".")
add("ABS-7", "Abstract, L53",
    "`the Fourier-Spectral-PINN's symmetry violation remains at machine precision across all tested scales, whereas a vanilla PINN's violation grows measurably with time`",
    "Fourier S = 8.71e-32 ... 3.74e-31 over every configuration/time stored (experiment_2d.json, stress_2d_sweep.json). Vanilla (5 seeds) mean S = 4.12e-5 (t=0) -> 2.22e-4 (t=1). The vanilla PINN was run only at the baseline eps = 0.3, not at the stress configurations.",
    "OK", "-")
add("ABS-8", "Abstract L53; Sec. 3.3 L243, L304",
    "`the energy-constrained loss keeps enforcing zero energy violations even once solution accuracy has fully collapsed`",
    "stress_time_sweep.json: n_energy_violations = 0 of 201 at T=10 with L2(t=10) = 0.990 (count OK). No run WITHOUT the penalty at T=5 or T=10 and no L_E log for the T-sweep, so \"enforcing\" cannot be attributed to the penalty. (Unpenalised K=32 net already gave 0/51 at T=1.)",
    "UNVERIFIABLE",
    "\"the network's own energy is non-increasing at all 201 test times although the L2 error is ~1 (no unpenalised control was run at T>1)\", or run the control.")

# ============================================================ METHODS
add("MET-1", "Sec. 2.1, Eq. (AC), L86",
    "`u_x(\\pm1,t)=0, \\qquad u(x,0)=u_0(x)=x^2\\cos(\\pi x)`",
    "u0'(x) = 2x cos(pi x) - pi x^2 sin(pi x): u0'(-1) = +2, u0'(+1) = -2 (numerically +2.000000 / -2.000000; the K=32 Chebyshev projection of u0 reproduces u_x(+-1) = -+2.0000, fit error 2.3e-15). IC is NOT compatible with the Neumann condition (Reviewer 1, major ii).",
    "WRONG",
    "Either change the IC to a compatible one (u0'(+-1)=0) and rerun, or state explicitly that the IC violates the BC at t=0 (absolute value of u_x equal to 2) and discuss the resulting boundary layer. Do not write that the problem has u_x(+-1,t)=0 for all t >= 0 with this IC.")
add("MET-2", "Sec. 2.1, L88", "`we use $\\varepsilon=0.01$ and $T=1$`", "common.py: EPS_DEFAULT = 0.01, T_FINAL = 1.0; all forward scripts use EPS = 0.01.", "OK", "-")
add("MET-3", "Sec. 2.3, L111",
    "`obtained errors at the level of $10^{-11}$--$10^{-13}$ for $K=20$`",
    "No data file stores this check (test function not named). Independent reproduction with D_c (K=20): exp(x) 1.2e-12, sin(2x) 6.6e-13, sin(3x) 4.2e-13, cos(pi x) 1.1e-12, x^2 cos(pi x) 1.3e-12; non-smooth-at-this-K tanh(5x) 1.0e-1.",
    "OK", "Name the test function(s) used (the claim is only true for functions resolved at K=20).")
add("MET-4", "Sec. 2.5, L133",
    "`sampled at $n_x=40$ spatial points and five times $t\\in\\{0.2,0.4,0.6,0.8,1.0\\}$ ... $\\sigma\\in\\{0,0.01,0.05\\}$ ... learning rate $10\\times$ larger`",
    "inverse_problem.py: n_obs_x=40, obs_times=(0.2,...,1.0), log_eps group lr = 10*lr; run_experiments.run_inverse: noise in [0.0, 0.01, 0.05].",
    "OK", "-")
add("MET-5", "Sec. 2.5, L133",
    "`$(\\lambda_{\\rm PDE},\\lambda_{\\rm IC},\\lambda_{\\rm BC},\\lambda_{\\rm data})=(50,10,1,2)$`",
    "server_sweep_inverse.py (Table 2 sweep): lam_pde=50, lam_data=2, defaults lam_ic=10, lam_bc=1 -> matches. run_experiments.run_inverse (the only script producing the noise runs of L182/L203) uses recover_eps defaults lam_pde=10, lam_data=5, eps_init=0.03 -> does NOT match; sweep uses eps_init = 2*eps_true.",
    "OK", "True for Table 2. State the weights/eps_init of the noise runs separately (see INV-1).")
add("MET-6", "Sec. 2.5, L133",
    "`A comparatively large $\\lambda_{\\rm PDE}$ was necessary: with equal weighting the optimizer could satisfy $\\mathcal L_{\\rm data}$ almost exactly while leaving $\\varepsilon$ far from $\\varepsilon_{\\rm true}$`",
    "No log/ablation in data/ for equal weighting.",
    "UNVERIFIABLE", "Either supply the ablation (the weights ablation the reviewers ask for) or soften to \"in preliminary runs (not reported)\".")
add("MET-7", "Sec. 2.6, L139",
    "`single-core for individual runs reported in Section~\\ref{sec:cheb-baseline}--\\ref{sec:2d}`",
    "run_experiments.py: run_vanilla/run_spectral/run_inverse/run_2d all call torch.set_num_threads(threads) with default threads=2 (not 1).",
    "WRONG", "\"two CPU threads\" (or document how the single-core run was actually made).")
add("MET-8", "Sec. 2.6, L139",
    "`the multi-seed and multi-configuration sweeps of Section~3 were parallelized, two threads per run, across a 64-core / 62~GB server`",
    "Threads per worker: forward multiseed 2, inverse sweep 2, 2D multiseed 3 (server_sweep_2d), stress_k 3, stress_eps 3 (+BLAS cap 2), stress_time 4, stress_2d 3. Hardware (64 cores / 62 GB) is not recorded in any JSON.",
    "WRONG", "\"two threads per run for the forward and inverse sweeps, three to four for the stress tests\".")
add("MET-9", "Sec. 2.6, L139",
    "`pretrained for 500--1000 epochs` ... `$K=32$ Chebyshev modes` ... `$M=8$ modes per direction`",
    "warmup_epochs: 1000 (1D, spectral_pinn.py), 500 (2D); K=32 default; M_MODES=8.", "OK", "-")
add("MET-10", "Sec. 2.6, L139; Table 1 caption L157",
    "`Reference solutions use a semi-implicit finite-difference scheme with $N_x=500$, $\\delta t=10^{-4}$`",
    "common.reference_solution_fd defaults (Nx=500, dt_ref=1e-4) -> matches. Reference error (not reported): FD(500,1e-4) vs FD(4000,1e-5): L2 5.47e-5 (trapezoid), 7.74e-5 (200-pt protocol), Linf 4.85e-4; spatial part dominates (time-step part 7.6e-6).",
    "OK", "Add the reference error to the text (it is as large as the Chebyshev row, see T1-1).")
add("MET-11", "Sec. 2.6, L139",
    "`hidden widths 64--64, or 64--64--64 for the vanilla space-time PINN`",
    "CoeffNet hidden=(64,64); vanilla arch=(2,64,64,64,1); 2D vanilla [5,64,64,64,1]; Fourier hidden (64,64).", "OK", "-")
add("MET-12", "Sec. 2.6 (2D extension) L137; Appendix B L464",
    "S(t) definition, `$u_0(x,y)=0.4\\cos(2x)\\cos(2y)$`, `$\\varepsilon=0.3$`, `($n=49$)`, `$h=2\\pi/(n-1)$`",
    "experiment_2d_symmetry.py: EPS2D=0.3, u0_2d = 0.4 cos(2x) cos(2y), symmetry_violation_* n=49 endpoint-inclusive, S = mean((u-u[::-1,:])^2) + mean((u-u[:,::-1])^2).",
    "OK", "-")

# ============================================================ TABLE 1 and Sec. 3.1
add("T1-1", "Table 1, Chebyshev row, L162",
    "`Chebyshev spectral collocation ($N=80$) & 81 & $8.51\\times10^{-5}$ & $5.56\\times10^{-4}$ & 0.16`",
    "Entries reproduce chebyshev_baseline.json (N=80: l2 8.5110e-5, linf 5.5627e-4, 0.158 s). BUT the L2 is `l2_error(u1, u_ref, dx_ref)` = sqrt(sum_j e_j^2 * dx_ref) over the 81 Chebyshev nodes with dx_ref = 0.004 (weight sum 0.324 instead of 2); the PINN rows use 200 uniform points with dx = 2/199 (weight sum 2.010). Same solution, consistent norms (re-run with the project code): 200-pt PINN protocol **9.31e-5** (x1.09), Clenshaw-Curtis at nodes 6.80e-5 (x0.80), trapezoid on the 501-pt grid 7.06e-5 (x0.83). For N=40: 7.58e-5 stored vs 1.80e-4 (x2.4); N=60: 7.78e-5 vs 1.14e-4 (x1.5). Linf is unaffected (same on nodes and on the 200-pt grid). Details in Section 2 below.",
    "WRONG",
    "L2 = 9.3e-5 evaluated with the same 200-point protocol as the PINN rows (or give all rows on one common grid). Add that the FD reference itself deviates from a finer FD solution by L2 7.7e-5 / Linf 4.9e-4, i.e. this row measures the discrepancy of two discretisations; against FD(4000,1e-5) the N=80, dt=5e-4 solution has L2 4.3e-5, Linf 7.1e-5.")
add("T1-1b", "Table 1, Chebyshev Linf and time, L162",
    "`$5.56\\times10^{-4}$` and `0.16`",
    "Linf 5.5627e-4 on nodes = on 200-pt grid. FD reference's own Linf error (4.85e-4) is nearly the whole entry. Time 0.158 s -> Spectral-PINN/Chebyshev wall-clock ratio 169.4/0.158 = 1070; accuracy ratio 1.3785e-3 / 8.5e-5 = 16 (code norm) or 14.8 (consistent norm 9.31e-5) against the same reference.",
    "OK", "-")
add("T1-2", "Table 1, vanilla row, L163",
    "`8{,}577 & $(6.74\\pm0.55)\\times10^{-2}$ & $(8.68\\pm0.43)\\times10^{-2}$ & $919\\pm7$`",
    "L2 6.7377e-2 +- 5.458e-3; Linf 8.6794e-2 +- 4.320e-3; time 919.33 +- 6.98 (SD with n-1, n=5). Parameter count of [2,64,64,64,1] = 8577.",
    "OK", "State that SD uses n-1.")
add("T1-3", "Table 1, Spectral-PINN K=32, L164",
    "`$(1.38\\pm0.08)\\times10^{-3}$ & $(3.32\\pm0.48)\\times10^{-3}$ & $169\\pm16$`",
    "L2 1.3785e-3 +- 7.68e-5; Linf 3.3193e-3 +- 4.76e-4; time 169.42 +- 16.43.", "OK", "-")
add("T1-4", "Table 1, energy-constrained row, L165",
    "`$(1.38\\pm0.07)\\times10^{-3}$ & $(3.31\\pm0.48)\\times10^{-3}$ & $191\\pm27$`",
    "L2 1.3777e-3 +- 7.36e-5; Linf 3.3130e-3 +- 4.77e-4; time 191.13 +- 26.74.", "OK", "-")
add("T1-5", "Table 1, K=64 row, L166",
    "`$(6.00\\pm1.21)\\times10^{-2}$ & $(1.10\\pm0.22)\\times10^{-1}$ & $267\\pm2$`",
    "L2 6.0031e-2 +- 1.214e-2; Linf 1.1038e-1 +- 2.187e-2; time 266.55 +- 1.50. Config: warm-up 1000 + 6000 epochs (train_spectral_pinn default warmup_epochs=1000, server_sweep_forward epochs=6000), as stated in L171.", "OK", "-")
add("T1-6", "Table 1, DoF column, L162-166",
    "`81`, `8{,}577`, `$K{+}1=33$ coeff.`, `$K{+}1=65$ coeff.`",
    "Column mixes definitions. Trainable parameters: vanilla 8,577; Spectral-PINN K=32: 6,433 (1-64-64-33); K=64: 8,513 (1-64-64-65). Output coefficients: 33 / 65. Chebyshev: 81 nodal unknowns.",
    "WRONG", "Header \"Trainable parameters\": 81 (unknowns) / 8,577 / 6,433 / 6,433 / 8,513, with the coefficient counts (33, 65) in a separate column.")
add("T1-7", "Sec. 3.1, L154 and captions",
    "`mean $\\pm$ SD over 5 random seeds (all other hyperparameters fixed)`",
    "5 seeds [42, 123, 7, 2024, 99]; SD convention n-1 reproduces Tables 1 and 6 (not the Table 2 bottom row, see INV-6).", "OK", "State the SD convention once.")
add("T1-8", "Sec. 3.1 L171; Discussion L296",
    "`makes the Spectral-PINN roughly $40\\times$ \\emph{less} accurate than at $K=32$` / `degrades ... by roughly $40\\times$`",
    "6.0031e-2 / 1.3785e-3 = **43.55** (ratio of means). Same-seed ratios 42.6, 39.6, 59.4, 46.9, 30.7.", "OK",
    "\"about 44x (30-59x for individual seeds)\".")
add("T1-9", "Sec. 3.1, L171",
    "`consistently across all 5 seeds (SD comparable to the mean, but no seed approaches the $K=32$ accuracy)`",
    "K=64 SD/mean = 1.214e-2 / 6.003e-2 = **0.20** (vanilla 0.08, K=32 0.06). 'No seed approaches' is true: min K=64 L2 4.19e-2 vs max K=32 L2 1.45e-3 (29x).",
    "WRONG", "\"(SD about one fifth of the mean; the best K=64 seed, 4.2e-2, is still 29x worse than the worst K=32 seed)\".")
add("T1-10", "Sec. 3.1 L171; Discussion L296",
    "`compare $N=40,60,80$ in the underlying computation, all accurate to $\\sim10^{-4}$), ... more modes are unambiguously beneficial` / `where increasing $N$ from 40 to 80 is unambiguously beneficial (Section~\\ref{sec:cheb-baseline})`",
    "Stored (chebyshev_baseline.json): L2 7.58e-5, 7.78e-5, 8.51e-5 for N=40, 60, 80 -> error RISES with N; Linf 5.50e-4, 5.53e-4, 5.56e-4 (flat). N and dt changed together (2e-3, 1e-3, 5e-4) and the reference error (Linf 4.9e-4) equals the Linf entries. Under a consistent norm the L2 values are 1.80e-4, 1.14e-4, 9.31e-5 (decreasing), but that mixes N and dt effects. The N=40/60/80 numbers appear in no table, and Section 2.2 contains no such comparison.",
    "WRONG",
    "Delete \"unambiguously beneficial\". Either add a table (N, dt, L2/Linf on a common grid, reference error) or write \"errors against the FD reference are 7.6e-5-8.5e-5 for N=40-80 and are limited by the reference and the first-order time step\".")
add("T1-11", "Sec. 3.1, L171",
    "`scaling the training budget \\emph{in proportion to} $K$ up to $K=128$ does not restore the $K=32$ accuracy either`",
    "K=64, 2x: 2.25e-2 (16.4x worse than K=32 mean); K=128, 4x: 2.19e-1 (159x). Budget helps (x2.7 and x1.5 vs 1x) but does not restore accuracy.", "OK", "-")
add("T1-12", "Sec. 3.1, L178",
    "`Evaluated at 51 uniformly spaced times over $[0,1]$, the energy-constrained Spectral-PINN exhibited zero energy violations (an increase ... exceeding $10^{-10}$)`",
    "0 violations, 51 points (spacing 0.02), network's own E_theta (trapezoid on 300-pt grid, exact u_x), largest step change -1.285e-3 (all steps decrease), seed 42 only.", "OK", "Add \"(single seed)\".")
add("T1-13", "Sec. 3.1, L178",
    "`achieving essentially the same forward-problem accuracy as the unconstrained Spectral-PINN`",
    "Mean L2 1.3777e-3 vs 1.3785e-3 (ratio 0.9994); per seed identical to 1e-6.", "OK", "-")
add("T1-14", "Sec. 3.1, L178",
    "`the energy penalty ... therefore acts as a certificate rather than a competing objective in this regime`",
    "The unpenalised net satisfies the same property (0/51). Penalty inactive (see ABS-3b). Figure 1 curves differ by <= 4.5e-7 (relative 1.3e-6): E(0) = 0.400457, E(1) = 0.334773 / 0.334774.",
    "WRONG", "\"The penalty was inactive here: the unpenalised network already decreases its energy monotonically, which is why the two curves in Figure 1 coincide to 5e-7.\"")
add("T1-15", "Figure 1 / caption L175",
    "`Discrete Ginzburg--Landau energy along the Spectral-PINN trajectory, with and without the energy-constrained penalty`",
    "Data (E_traj, 51 pts) reproduce; the two curves coincide (max abs difference 4.54e-7). Also E[u0] (FD grid) = 0.400436 vs network E(0) 0.400457; reference E(1) = 0.334556 vs network 0.334773.",
    "OK", "Add the explanation requested by Reviewer 2 (curves coincide because the penalty is inactive).")

# ============================================================ INVERSE
add("INV-1", "Sec. 3.2, L182",
    "`$\\hat\\varepsilon=0.01402,0.01402,0.01400$ for $\\sigma=0,0.01,0.05$, i.e., $40.2\\%,40.2\\%,40.0\\%$ relative error`",
    "No data: run_experiments.py `inverse` would write data/inverse_problem.json, which does not exist. That script uses lam_pde=10, lam_data=5, eps_init=0.03 (not the (50,10,1,2) stated in Sec. 2.5). Nearest data: sweep K=32, eps=0.01 -> eps_hat 0.014044 (40.44 %); K=96 -> 0.014022 (40.22 %).",
    "UNVERIFIABLE", "Re-run and store, or drop these numbers; state the loss weights actually used.")
add("INV-2", "Sec. 3.2, L182",
    "`A single additional run at $\\varepsilon_{\\rm true}=0.05$ then gave a smaller, $14.7\\%$ error`",
    "No data. The sweep at eps_true=0.05, K=32 gives 43.36 % (different weights/eps_init).",
    "UNVERIFIABLE", "Remove or reproduce; as written it contradicts Table 2 at the nominal configuration.")
add("INV-3", "Sec. 3.2, L184",
    "`the full $5\\times4$ grid ... (20 independent noise-free runs, parallelized across a 64-core server)`",
    "inverse_bias_sweep.json: 20 entries, all seed 42, eps_init = 2*eps_true, 4000 epochs, noise 0.",
    "OK", "\"20 single-seed runs (one per configuration)\"; they are not replicates.")
add("INV-4", "Sec. 3.2 L184; Table 2 caption L187",
    "`$37.8\\pm7.0\\%$ (mean $\\pm$ SD)` / `overall mean $37.8\\pm7.0\\%$`",
    "mean 37.843, SD(n-1) = 6.977 (SD(n) = 6.800). This is the spread over different configurations, not a seed-to-seed SD.", "OK", "Say \"SD across the 20 configurations\".")
add("INV-5", "Table 2, 20 cells, L192-196",
    "e.g. `0.005 & 30.0 & 35.9 & 37.0 & 44.6` ... `0.1 & 43.6 & 43.3 & 37.3 & 40.0`",
    "All 20 cells equal rel_err_pct rounded to 0.1 (30.04, 35.90, 37.01, 44.59; 12.61, 40.44, 39.39, 40.22; 35.34, 37.22, 35.43, 39.42; 44.27, 43.36, 37.00, 40.33; 43.63, 43.35, 37.33, 39.99).",
    "OK", "The table lists absolute errors; the 12.6 entry is signed -12.6 (underestimate).")
add("INV-6", "Table 2, bottom row, L198",
    "`Mean $\\pm$ SD & $33.2\\pm11.9$ & $40.0\\pm3.1$ & $37.2\\pm1.5$ & $40.9\\pm2.2$`",
    "Means 33.18, 40.05, 37.23, 40.91. SD(n-1) = 12.94, 3.44, 1.41, 2.09; SD(n) = 11.58, 3.07, 1.27, 1.87. The printed SDs (11.9, 3.1, 1.5, 2.2) match neither convention (only 3.1 matches SD(n)); 40.05 rounds to 40.1, not 40.0. The overall 7.0 uses n-1.",
    "WRONG", "`33.2 +- 12.9`, `40.1 +- 3.4`, `37.2 +- 1.4`, `40.9 +- 2.1` (SD with n-1, as in the other tables).")
add("INV-7", "Sec. 3.2, L184",
    "`\\emph{no} systematic decrease as $K$ increases from 16 to 96 at any fixed $\\varepsilon_{\\rm true}$`",
    "Net change K=16 -> 96 (percentage points) for eps = 0.005, 0.01, 0.02, 0.05, 0.1: +14.6, +27.6, +4.1, -3.9, -3.6. At eps=0.05 and 0.1 the error falls 44.3 -> 37.0 and 43.6 -> 37.3 up to K=64, then rises. No monotone decrease anywhere. One seed per cell, so differences of a few points are not resolvable.",
    "OK", "Add \"within single-seed resolution\".")
add("INV-8", "Sec. 3.2, L184",
    "`and no consistent trend across $\\varepsilon_{\\rm true}$ at any fixed $K$`",
    "corr(log eps, error) = +0.73 (K=16), +0.84 (K=32), -0.18 (K=64), -0.66 (K=96). At K=16 and 32 the error is 43-44 % for eps >= 0.05 and 30-40 % for eps <= 0.02.",
    "WRONG", "\"only a weak, non-monotone dependence on eps_true (somewhat larger bias for eps >= 0.05 at K = 16, 32), within what a single seed per cell can resolve\".")
add("INV-9", "Sec. 3.2, L184",
    "`the one low outlier, $12.6\\%$ at $\\varepsilon_{\\rm true}=0.01,K=16$, is not reproduced at any neighboring configuration`",
    "12.61 %; neighbours 30.04 (eps=0.005), 40.44 (K=32), 35.34 (eps=0.02). eps_hat = 0.008739 < 0.01: this run UNDER-estimates (signed -12.6 %), so the 'overestimation' narrative has one sign reversal.",
    "OK", "Report the sign.")
add("INV-10", "Sec. 3.2 L184; Discussion L306; Conclusions L311",
    "`the inverse solver has a robust $\\sim38\\%$ overestimation bias ... a property of the estimator as a whole`",
    "Numeric part reproduced (INV-4). 'Robust' / 'property of the estimator' is not testable from the stored data: one seed per cell, no eps trajectories stored in the sweep file (only final eps_hat), all runs start at eps_init = 2*eps_true and 19 of 20 end at 1.30-1.45 x eps_true (between truth and init), so bias cannot be separated from incomplete convergence in 4000 epochs. No ablation of weights / warm-up / observation design exists (Reviewer 2.5).",
    "UNVERIFIABLE", "Replace by \"in this configuration the recovered eps was 30-45 % too large in 19 of 20 single-seed runs started at 2 eps_true; convergence and the cause were not investigated\" unless the ablation is added.")
add("INV-11", "Sec. 3.2, L203",
    "`recovery under $1\\%$ and $5\\%$ observation noise ($40.2\\%$ and $40.0\\%$ relative error) is indistinguishable from the noise-free case ($40.2\\%$)`",
    "No data (see INV-1). The sweep's own noise-free K=32, eps=0.01 value is 40.44 %, not 40.2 %. Noise is sigma = 0.01 / 0.05 absolute (u in [-1,1]), not '1 % / 5 %' relative.",
    "UNVERIFIABLE", "Rerun with the sweep configuration and store the JSON; write sigma = 0.01 and 0.05.")
add("INV-12", "Discussion L302; Conclusions L311 (abstract L53 says '6x')",
    "`by varying $K$ over more than a $6\\times$ range` / `across more than a $6\\times$ range of $K$`",
    "K in {16, 32, 64, 96}: 96/16 = 6.0 exactly.", "WRONG", "\"a 6x range\" (as in the abstract).")

# ============================================================ STRESS K
add("STK-1", "Table 3, all rows, L217-222; caption L212",
    "`64 & 1$\\times$ & 1000 + 6000 & $6.18\\times10^{-2}$ & 178` ... `256 & 8$\\times$ matched & 8000 + 48000 & $1.47\\times10^{0}$ & 4659`; caption baseline `$L^2=1.4\\times10^{-3}$`",
    "stress_k_sweep.json: L2 0.061756, 0.022548, 0.336276, 0.219080, 1.510083, 1.472173; epochs/warm-up as printed; times 178.4, 364.6, 351.2, 1387.1, 777.6, 4659.4 s. Baseline 1.38e-3 (Table 1 mean).", "OK", "-")
add("STK-2", "Sec. 3.3, L209",
    "`$8\\times$ budget ($8000$ warm-up $+48{,}000$ epochs, $4659\\,$s of training) yields $L^2=1.47$ ... ($L^2=1.51$ after $778\\,$s)`",
    "As in STK-1. L2 ratio 1x/8x = 1.026; Linf 5.46 vs 5.09.", "OK", "-")
add("STK-3", "Sec. 3.3, L209",
    "`statistically indistinguishable from the $1\\times$-budget result`",
    "One run per budget (seed 42); no variance available, so no statistical statement is possible. A 2.6 % difference (1.47 vs 1.51) is smaller than the 5.6 % seed-to-seed relative SD seen at K=32.",
    "UNVERIFIABLE", "\"essentially unchanged (1.47 vs 1.51, single runs)\".")
add("STK-4", "Sec. 3.3, L209",
    "`an extra $\\sim$78 minutes of training produced no meaningful improvement at all`",
    "8x run: 4659.4 s = 77.7 min TOTAL; extra over the 1x run (777.6 s = 13.0 min) = 3881.9 s = **64.7 min**.",
    "WRONG", "\"an extra ~65 minutes of training (78 minutes in total)\".")
add("STK-5", "Sec. 3.3, L209",
    "`plausibly related to the $O(K^4)$ growth of the condition number of Chebyshev second-derivative operators`",
    "Not a data claim. Numerically the 2-norm of D_c^2 = 1.49e5, 2.27e6, 3.53e7, 5.58e8 for K = 32, 64, 128, 256 (log-log slope 3.96); nodal D^2 spectral radius slope 4.00. D_c^2 is nilpotent, so the 'condition number' of D_c^2 is undefined; the K^4 growth is of the norm.",
    "OK", "Say \"norm\" instead of \"condition number\" for the coefficient-space operator.")

# ============================================================ STRESS EPS
add("STE-1", "Table 4, vanilla and Spectral-PINN columns, L235-238",
    "`$6.75\\times10^{-2}$ & $1.45\\times10^{-3}$` ... `$6.81\\times10^{-2}$ & $9.77\\times10^{-4}$`",
    "stress_eps_sweep.json: vanilla 0.067545, 0.067911, 0.068071, 0.068145; Spectral K=32 1.4522e-3, 9.888e-4, 9.799e-4, 9.768e-4 -> all equal the printed values.", "OK", "-")
add("STE-2", "Table 4, Chebyshev column, L235-238",
    "`$1.17\\times10^{-4}$`, `$3.39\\times10^{-5}$`, `$7.05\\times10^{-6}$`, `$2.66\\times10^{-6}$`",
    "Values equal stress_eps_sweep.json, but l2_error(u1, u_ref, x_ref[1]-x_ref[0]) is applied to N+1 = 401 nodes, so the weight sum is 401*dx_ref = 1.201, 0.601, 0.241, 0.120 (should be 2). Re-run: stored vs 200-pt PINN protocol 1.17e-4 vs 6.68e-5; 3.39e-5 vs 2.96e-5; 7.05e-6 vs 1.32e-5; 2.66e-6 vs 9.44e-6. Clenshaw-Curtis at nodes: 4.69e-5, 1.63e-5, 7.82e-6, 7.74e-6. The stored column falls 44x over the eps range; with a consistent norm it falls 7x.",
    "WRONG", "Replace the column by the 200-point-protocol values 6.7e-5, 3.0e-5, 1.3e-5, 9.4e-6 (or give all columns on the same refining grid).")
add("STE-3", "Sec. 3.3, L227",
    "`scaling the fine-grid reference and Chebyshev-baseline resolutions with $1/\\varepsilon$ so the ground truth itself remains resolved`",
    "FD reference nx_ref = 668, 1334, 3334, 6668 (scales). Chebyshev N = 400 for ALL four eps (stress_common.chebyshev_n_for_eps caps at 400); only dt is scaled (1e-4 ... 1e-5). For eps <= 0.005 the cap leaves N below the nominal 2/(0.5 eps) = 400, 800, 2000, 4000.",
    "WRONG", "\"scaling the FD reference grid with 1/eps (N_x = 668 ... 6668); the Chebyshev baseline uses a fixed N = 400 with dt scaled with eps\".")
add("STE-4", "Sec. 3.3, L227",
    "`the classical Chebyshev baseline's own error \\emph{decreases} with $\\varepsilon$ ... only because it is evaluated at its own (comparably coarse, $N=400$) nodes rather than this fixed grid`",
    "Main cause of the 44x decrease is the node-count x dx_ref weighting (effective measure 1.20 -> 0.12, STE-2), not the evaluation grid per se.",
    "WRONG", "Explain the decrease by the inconsistent quadrature weight, or remove the Chebyshev column.")
add("STE-5", "Sec. 3.3, L227",
    "`essentially \\emph{flat} across this $10\\times$ range of $\\varepsilon$ ($6.8\\times10^{-2}$ and $\\sim1\\times10^{-3}$, respectively ...)` ; `(grid spacing $\\approx0.01$)` ; `$N=400$`",
    "0.01 -> 0.001 is 10x. Vanilla 6.75e-2 ... 6.81e-2. Spectral 1.45e-3, 9.9e-4, 9.8e-4, 9.8e-4. Evaluation grid np.linspace(-1,1,200): spacing 2/199 = 0.01005. Chebyshev N = 400.", "OK", "-")

# ============================================================ STRESS T
add("STT-1", "Table 5, L251-253; epochs text L243",
    "`1 & 6000 & $1.5\\times10^{-3}$ & $1.5\\times10^{-3}$ & 0 / 21` ... `10 & 19{,}500 & $2.1\\times10^{-1}$ & $9.9\\times10^{-1}$ & 0 / 201`",
    "stress_time_sweep.json: L2(t=1) = 1.4517e-3, 0.11237, 0.20610; L2(t=T) = 1.4517e-3, 0.90369, 0.99004; violations 0 of 21, 101, 201; epochs 6000 = 6000+1500(T-1) -> 6000, 12000, 19500.", "OK", "-")
add("STT-2", "Sec. 3.3, L243",
    "`grows from $1.5\\times10^{-3}$ ... to $0.11$ ... to $0.21$ ... by $t\\gtrsim5$ within the $T=10$ run the error saturates near $1.0$`",
    "0.0015 -> 0.112 -> 0.206. T=10 run: L2 = 0.987 (t=5), 1.022, 1.031, 1.019, 1.001, 0.990 (t=6..10).", "OK", "-")
add("STT-3", "Sec. 3.3, L243",
    "`zero energy violations were observed at any of the $21$, $101$, or $201$ evaluation points for $T=1,5,10$`",
    "n_eval_points = 21, 101, 201 (spacing 0.05); count = increases of the network's own E_theta between consecutive points above 1e-10. Penalty checkpoints are every 0.1 (n_ckpt = 10, 50, 100), i.e. half of the test points were never seen by the penalty.", "OK", "State the evaluation spacing and the 1e-10 tolerance in the text.")
add("STT-4", "Discussion, L304",
    "`(zero violations at 201 checkpoints)`",
    "201 is the number of evaluation points (T/0.05 + 1). The penalty uses 100 intervals = 101 checkpoint times at T=10.",
    "WRONG", "\"at all 201 evaluation points (the penalty itself acts on 101 checkpoints)\".")
add("STT-5", "Discussion, L304",
    "`by $t\\approx5$, an $L^2$ error of essentially $O(1)$`", "0.987 at t=5 (T=10 run); 0.904 at t=5 (T=5 run).", "OK", "-")

# ============================================================ 2D
add("S2D-1", "Sec. 3.3, L258",
    "`the Fourier-Spectral-PINN's violation stays at machine precision ($\\sim10^{-32}$--$10^{-31}$) at every $t$`",
    "experiment_2d.json: 8.71e-32 (t=0), 9.41e-32, 1.03e-31, 1.26e-31 (t=1).", "OK", "-")
add("S2D-2", "Sec. 3.3, L258",
    "`the vanilla PINN's grows monotonically, by roughly an order of magnitude over $t\\in[0,1]$, in every one of the 5 seeds`",
    "Monotone at t = 0, 0.2, 0.5, 1 in all 5 seeds: TRUE. Growth S(1)/S(0) per seed 4.36, 6.79, 6.14, 6.20, 4.47; ratio of means 5.40 (0.73 decades).",
    "WRONG", "\"by a factor of about 5 (4.4-6.8 across seeds)\".")
add("S2D-3", "Table 6, Fourier column, L266-269",
    "`$8.7\\times10^{-32}$`, `$9.4\\times10^{-32}$`, `$1.0\\times10^{-31}$`, `$1.3\\times10^{-31}$`",
    "8.711e-32, 9.406e-32, 1.029e-31, 1.263e-31 (single run, seed 42).", "OK", "-")
add("S2D-4", "Table 6, vanilla column, L266-269",
    "`$(4.12\\pm1.13)\\times10^{-5}$`, `$(4.94\\pm1.16)\\times10^{-5}$`, `$(9.46\\pm1.81)\\times10^{-5}$`, `$(2.23\\pm0.45)\\times10^{-4}$`",
    "symmetry_2d_multiseed.json (n-1): 4.122e-5 +- 1.125e-5; 4.938e-5 +- 1.160e-5; 9.463e-5 +- 1.812e-5; **2.225e-4** +- 4.489e-5. Last mean rounds to 2.22, not 2.23.",
    "WRONG", "`(2.22 +- 0.45)e-4` (minor rounding slip; the same quantity is written (2.2 +- 0.4)e-4 in L300, which is consistent).")
add("S2D-5", "Discussion, L300",
    "`from $(4.1\\pm1.1)\\times10^{-5}$ at $t=0$ to $(2.2\\pm0.4)\\times10^{-4}$ at $t=1$` ; `$S(t)\\sim10^{-31}$--$10^{-32}$, Table~\\ref{tab:2d}`",
    "4.12 +- 1.13e-5; 2.22 +- 0.45e-4 (0.449 -> 0.4). Fourier S in [8.7e-32, 1.3e-31].", "OK", "-")
add("S2D-6", "Sec. 3.3, L274 (reviewer minor iii; PDF line 390)",
    "`$S(t)$ remains at $10^{-31}$--$10^{-31.4}$ across every configuration` ... `completely scale-invariant`",
    "Table 7 / stress_2d_sweep.json at t=1: S = 1.257e-31 ... 3.743e-31 = 10^-30.90 ... 10^-30.43. Over all stored times (0, 0.5, 1): 8.81e-32 ... 3.74e-31 = 10^-31.06 ... 10^-30.43. The quoted upper end 10^-31.4 = 4.0e-32 is 9.4x below the actual maximum. Also S(1) at eps=0.1 is 2.9x that at eps=0.3 (3.74e-31 vs 1.27e-31), so 'completely scale-invariant' is imprecise.",
    "WRONG", "\"S(1) stays between 1.3e-31 and 3.7e-31 (about 10^-30.9 to 10^-30.4) in all five configurations, i.e. at round-off level (about 3x larger for eps = 0.1 than for eps = 0.3)\".")
add("S2D-7", "Table 7, L282-286",
    "`8 & 0.3 & 3000 & $1.27\\times10^{-31}$` ... `16 & 0.1 & 6000 & $3.74\\times10^{-31}$`",
    "S(1) = 1.2734e-31, 1.2718e-31, 1.2569e-31, 3.7426e-31, 3.7421e-31; epochs 3000, 3000, 6000, 3000, 6000; M, eps as printed.", "OK", "-")
add("S2D-8", "Table 7 row 1 vs Table 6",
    "`8 & 0.3 & 3000 & $1.27\\times10^{-31}$ & (baseline, cf.\\ Table~\\ref{tab:2d})`",
    "Table 6 has 1.263e-31 (printed 1.3e-31) for the same configuration: a separate re-run (stress_2d_sweep.json) differing by 0.8 %.", "OK", "Note that the baseline was re-run.")
add("S2D-9", "Sec. 3.3, L274",
    "`with the training budget also doubled for the larger-$M$ runs`",
    "M=16 runs: 3000 (same budget), 6000, 6000 -> doubled in 2 of 3.", "WRONG", "\"doubled in two of the three M=16 runs\".")
add("S2D-10", "Discussion, L304",
    "`holding at $10^{-31}$ regardless of $M$ or $\\varepsilon$`", "1.3e-31 ... 3.7e-31 at t=1 (order of magnitude 10^-31 correct).", "OK", "-")

add("APP-1", "Sec. 2.7 L143; Appendix A L352-377",
    r"`$2a(a-b)=a^2-b^2+(a-b)^2$`; Gronwall step `$a^2\le\frac{1+\tau(L+1)}{1-\tau(L+1)}b^2$`; `$\tau=\frac{1}{2(L+1)}$`",
    "Hand check: a^2-b^2+(a-b)^2 = 2a^2-2ab = 2a(a-b); a^2 <= b^2 + c(a^2+b^2) <=> (1-c)a^2 <= (1+c)b^2, divide by 1-c>0; tau(L+1) = 1/2 < 1. Lean toolchain file = leanprover/lean4:v4.33.1 (as stated). Compilation / axiom claims not re-run here.",
    "OK", "-")

# ============================================================ CONCLUSIONS
add("CON-1", "Conclusions, L311",
    "`improves forward-problem accuracy over a vanilla space-time PINN by a factor of about 47`",
    "48.9 (ratio of 5-seed means; abstract says 49). The value 47 is reproduced only by single seed-42 pairs: 46.5 (forward_multiseed.json) / 46.6 (vanilla_pinn.json vs spectral_pinn.json).",
    "WRONG", "\"about 49\" (and cite the same statistic as in the abstract).")
add("CON-2", "Conclusions, L311", "`while training faster`", "5.4x (4.8x for energy-constrained), see ABS-2.", "OK", "-")

R_ORDER = [r[0] for r in R]

cnt = collections.Counter(r[4] for r in R)
wrong = [r for r in R if r[4] == "WRONG"]
unv = [r for r in R if r[4] == "UNVERIFIABLE"]

out = []
w = out.append
w("# Numbers audit: manuscript.tex vs data/*.json and code/*.py")
w("")
w("Manuscript: `manuscript/manuscript.tex` (symmetry-4596052, round-1 version). Every quantitative claim was recomputed with the scripts in this folder from `data/*.json` (read-only); the way each error is computed was checked in `code/` and, where in doubt, the code was re-executed (read-only import, no file in `code/`, `data/`, `manuscript/` was modified).")
w("Line numbers L<n> refer to `manuscript.tex`. The reviewer's \"line 390 / Table 7\" is the PDF line number; it corresponds to L274 and `tab:stress-2d`.")
w("")
w("## 0. Summary")
w("")
w("| status | count |")
w("|---|---|")
w("| OK | %d |" % cnt["OK"])
w("| WRONG | %d |" % cnt["WRONG"])
w("| UNVERIFIABLE (no data saved / no control run / not testable) | %d |" % cnt["UNVERIFIABLE"])
w("| total claims audited | %d |" % len(R))
w("")
w("WRONG items: " + ", ".join(r[0] for r in wrong) + ".")
w("")
w("UNVERIFIABLE items: " + ", ".join(r[0] for r in unv) + ".")
w("")
w("### The findings that matter most for the revision")
w("")
w("1. **49 vs 47 (Reviewer 2.1).** 49 = ratio of the 5-seed mean L2 errors (48.88; 95 % bootstrap interval [45.0, 52.6]). 47 = the single seed-42 pair (46.5-46.6). Per-seed paired ratios range 42.5-54.7. Use 49 everywhere and name the statistic (CON-1, ABS-1).")
w("2. **Chebyshev row of Table 1 uses a different norm than the PINN rows (T1-1).** `l2_error(u1, u_ref, dx_ref)` is applied to 81 Chebyshev nodes with the 501-point grid spacing: weight sum 0.324 instead of 2. Effect on the headline N=80 entry is moderate and of mixed sign (8.51e-5 stored; 9.31e-5 with the PINN protocol, 6.8e-5 with Clenshaw-Curtis weights), but for N=40 it is x2.4 and for the eps-sweep column of Table 4 it is up to x0.28 (STE-2). The N=40/60/80 'more modes are beneficial' statement is contradicted by the stored numbers (T1-10). The FD reference itself is only accurate to L2 7.7e-5 / Linf 4.9e-4, which is the same size as the Chebyshev 'error'.")
w("3. **The range of S(t) (Reviewer minor iii) is wrong (S2D-6).** The sentence says 10^-31 ... 10^-31.4; Table 7 gives 1.26e-31 ... 3.74e-31 at t=1 (10^-30.9 ... 10^-30.4).")
w("4. **'Zero energy violations under the penalty' is not evidence for the penalty (ABS-3b, T1-14, ABS-8).** The unpenalised network has 0/51 as well; the penalty term was non-zero in 1 of 120 logged epochs (max 7.6e-8). No multi-seed and no unpenalised long-horizon control exists. lambda_E = 50 and M = 10 checkpoints (code) appear nowhere in the manuscript (Section 4 below).")
w("5. **The noise-robustness numbers and the 14.7 % single run in Sec. 3.2 have no data file (INV-1, INV-2, INV-11).** `data/inverse_problem.json` is missing, and the script that makes those runs uses different loss weights and eps_init than the text states.")
w("6. **Table 2 bottom-row SDs are inconsistent (INV-6)**, and all 20 inverse runs start from eps_init = 2 eps_true with one seed and no stored trajectories, so the '38 % robust bias' cannot be separated from incomplete convergence (INV-10). One of the 20 cells is an under-estimate (INV-9).")
w("7. **The initial condition violates the Neumann condition** (u0'(+-1) = -+2; MET-1).")
w("")
w("## 1. Claim-by-claim table")
w("")
w("Status: **OK** = reproduced from data (rounding included); **WRONG** = contradicted by data/code or internally inconsistent; **UNVERIFIABLE** = no data file, no control run, or statistically untestable as worded. Statements that need no number (qualitative narrative, literature) are not listed.")
w("")
w("| id | location | claim as written (LaTeX) | recomputed from data/code | status | proposed corrected wording |")
w("|---|---|---|---|---|---|")
for (i, loc, claim, val, st, fix) in R:
    stx = {"OK": "OK", "WRONG": "**WRONG**", "UNVERIFIABLE": "**UNVERIFIABLE**"}[st]
    w("| %s | %s | %s | %s | %s | %s |" % (i, loc, claim, val, stx, fix))
w("")

# ---------------------------------------------------------------- norm section
w("## 2. How each 'error' is computed (code audit of the norms)")
w("")
w("| table / rows | script | grid on which e = u - u_ref is formed | L2 formula | weight sum | L-inf |")
w("|---|---|---|---|---|---|")
w("| Table 1 vanilla, Spectral K=32, energy, K=64 (`server_sweep_forward.py`); Tables 3, 5 (`stress_k_sweep.py`, `stress_time_sweep.py`); Table 4 vanilla/Spectral columns (`stress_eps_sweep.py`) | PINN prediction on `x_test = linspace(-1,1,200)`; reference = `np.interp` of the FD solution (Nx = 500, or eps-scaled) | 200 uniform points incl. both ends | `sqrt(sum(e^2) * dx)`, dx = 2/199 (Riemann sum) | 200 dx = 2.0101 (true 2) | max over the 200 points |")
w("| Table 1 Chebyshev row (`run_experiments.run_chebyshev_baseline`, `chebyshev_spectral.py __main__`) | Chebyshev-Gauss-Lobatto nodes x_j = cos(j pi/N), N = 40/60/80; reference interpolated to the nodes | N+1 nodes, strongly non-uniform (clustered at +-1) | `sqrt(sum(e^2) * dx_ref)` with dx_ref = x_ref[1]-x_ref[0] = 0.004 of the **501-point uniform FD grid** | (N+1) dx_ref = 0.164, 0.244, **0.324** (true 2) | max over the nodes |")
w("| Table 4 Chebyshev column (`stress_eps_sweep.run_chebyshev`) | N = 400 nodes, eps-scaled FD reference | 401 nodes | same formula, dx_ref = 2/nx_ref | 1.201, 0.601, 0.241, 0.120 (true 2) | max over the nodes |")
w("")
w("So the Chebyshev numbers are **not** computed with the norm of the PINN rows. Re-running exactly the project's solvers (`norm_check.py`, output in `norm_check_output.txt`):")
w("")
w("**Table 1 Chebyshev row, eps = 0.01, reference FD(500, 1e-4):**")
w("")
w("| N | dt | stored (code formula) | PINN protocol (interp. to 200 pts, dx = 2/199) | Clenshaw-Curtis at nodes | trapezoid on 501-pt grid | PINN-protocol / stored | CC / stored |")
w("|---|---|---|---|---|---|---|---|")
w("| 40 | 2e-3 | 7.578e-5 | 1.798e-4 | 1.663e-4 | 1.707e-4 | 2.37 | 2.19 |")
w("| 60 | 1e-3 | 7.780e-5 | 1.140e-4 | 9.361e-5 | 9.746e-5 | 1.46 | 1.20 |")
w("| 80 | 5e-4 | 8.511e-5 | 9.310e-5 | 6.796e-5 | 7.057e-5 | 1.09 | 0.80 |")
w("")
w("Linf is the same in all protocols (5.50e-4, 5.53e-4, 5.56e-4). The mismatch factor is not a constant (x2.4 / x1.5 / x1.1) because it depends on where the error sits relative to the node clustering: the maximum error is at the boundary node x = +1 (5.5e-4 for N=40 and 80; `err_location_output.txt`), exactly where the FD(500) reference deviates most from FD(4000) (4.85e-4 at x = -1, the boundary at which the IC violates the Neumann condition), and the share of the Clenshaw-Curtis-weighted squared error in |x| > 0.9 is 13.5 % for N=40 but 73.7 % for N=80 (30 % of the nodes lie there). The headline N = 80 entry is accidentally within +9 % / -20 % of the consistent values.")
w("")
w("**Error of the reference itself (the 'fine-grid' FD with Nx = 500, dt = 1e-4):**")
w("")
w("| comparison | L2 | Linf |")
w("|---|---|---|")
w("| FD(500,1e-4) vs FD(4000,1e-5), trapezoid on the 4001-pt grid | 5.47e-5 | 4.85e-4 |")
w("| FD(500,1e-4) vs FD(4000,1e-5), 200-pt PINN protocol | 7.74e-5 | 4.85e-4 |")
w("| FD(1000,5e-5) vs FD(4000,1e-5) | 2.41e-5 | 2.16e-4 |")
w("| FD(2000,2.5e-5) vs FD(4000,1e-5) | 8.27e-6 | 7.32e-5 |")
w("| pure time-step effect, FD(500,1e-5) vs FD(500,1e-4) | 7.59e-6 | 1.17e-5 |")
w("| pure spatial effect, FD(2000,1e-4) vs FD(500,1e-4) | 4.58e-5 | 4.12e-4 |")
w("| Chebyshev N=80, dt=5e-4 vs FD(4000,1e-5), 200-pt protocol | 4.27e-5 | 7.11e-5 |")
w("| Chebyshev N=80, dt=1e-4 vs FD(4000,1e-5), 200-pt protocol | 1.38e-5 | 7.37e-5 |")
w("")
w("Consequences: (i) the reference error (L2 5.5e-5-7.7e-5, Linf 4.9e-4) is as large as the Chebyshev 'error' in Table 1, so that row measures the discrepancy of two discretisations (spatial error of the FD reference dominates); against the finer reference the N=80, dt=5e-4 Chebyshev solution has L2 4.3e-5 and Linf 7.1e-5, still 30x more accurate in L2 than the Spectral-PINN (1.38e-3). The qualitative message (classical spectral collocation is both faster and more accurate) is unchanged; (ii) the PINN rows are not affected: seed-42 Spectral-PINN L2 is 1.4505e-3 against the manuscript reference and 1.4569e-3 against FD(4000,1e-5); vanilla 6.7613e-2 vs 6.7631e-2; the Riemann sum vs trapezoid on the 200 points changes the L2 by 1e-4 relative.")
w("")
w("**Table 4 Chebyshev column (N = 400, eps-scaled reference; re-run reproduces the stored values to 0.3 %):**")
w("")
w("| eps | nx_ref | weight sum 401 dx_ref | stored L2 | PINN protocol (200 pts) | CC at nodes | trapezoid on ref grid | stored / PINN-protocol |")
w("|---|---|---|---|---|---|---|---|")
w("| 0.01 | 668 | 1.201 | 1.174e-4 | 6.678e-5 | 4.687e-5 | 4.779e-5 | 1.76 |")
w("| 0.005 | 1334 | 0.601 | 3.385e-5 | 2.956e-5 | 1.635e-5 | 1.679e-5 | 1.15 |")
w("| 0.002 | 3334 | 0.241 | 7.053e-6 | 1.315e-5 | 7.821e-6 | 7.880e-6 | 0.54 |")
w("| 0.001 | 6668 | 0.120 | 2.659e-6 | 9.441e-6 | 7.744e-6 | 7.744e-6 | 0.28 |")
w("")
w("Under the stored (inconsistent) formula the column drops 44x over the eps range; under a consistent formula it drops 7x (6.7e-5 -> 9.4e-6). The manuscript's explanation (L227) is therefore incomplete: the effective measure shrinks from 1.2 to 0.12 because 401 nodes are weighted by the spacing of a reference grid that is refined with 1/eps.")
w("")
w("Other norm remarks: all PINN rows use the same 200 uniform points, so they are mutually consistent; the L2 there is a Riemann sum with end-point weight dx (trapezoid would change it by about 1e-4 relative). The 200-point grid (spacing 0.01) under-resolves interfaces of width O(eps) for eps <= 0.01, as the manuscript itself says for Table 4; the same applies, less severely, to Tables 1, 3, 5 at eps = 0.01.")
w("")

# ---------------------------------------------------------------- no-data section
w("## 3. Claims with NO supporting data file")
w("")
w("| claim (location) | why unsupported | what exists |")
w("|---|---|---|")
w("| Noise-level results eps_hat = 0.01402 / 0.01402 / 0.01400, 40.2 / 40.2 / 40.0 % (L182, L203) | `data/inverse_problem.json` (written by `run_experiments.py inverse`) is not in `data/`; the script's defaults (lam_pde=10, lam_data=5, eps_init=0.03) differ from the stated (50,10,1,2) | sweep values at eps = 0.01: 12.6 / 40.4 / 39.4 / 40.2 % (K = 16/32/64/96) |")
w("| Single run at eps_true = 0.05 gives 14.7 % (L182, L302) | no data | sweep K=32, eps=0.05: 43.4 % |")
w("| 'With equal weighting the optimizer could satisfy L_data almost exactly while leaving eps far from eps_true' (L133) | no ablation file | - |")
w("| lambda_E ablation / effect of lambda_E (requested by Reviewer 2.2) | none was run; lambda_E = 50 only exists as a default in `spectral_pinn.py` | - |")
w("| Zero energy violations in the 5-seed runs | `forward_multiseed.json` stores only l2/linf/elapsed | seed 42 only: `spectral_pinn_energy.json`, `stress_time_sweep.json` |")
w("| Zero violations 'enforced by the penalty' at T = 5, 10 (L243, L304) | no unpenalised run and no L_E log at T > 1 | - |")
w("| Energy-penalty activity (certificate vs competing objective, L178) | only `energy_history` of the single seed-42 run | L_E nonzero at 1 of 120 logged epochs |")
w("| 'Chebyshev N = 40, 60, 80 ... unambiguously beneficial' (L171, L296) | numbers exist only in `chebyshev_baseline.json`, in no table, and they contradict the claim | see T1-10 |")
w("| D_c derivative accuracy 1e-11 ... 1e-13 for K = 20 (L111) | no stored result, test function unnamed | reproduced independently (4e-13 ... 1.3e-12) |")
w("| Hardware/threads: 64 cores / 62 GB; 'single-core' individual runs (L139) | no log in any JSON; code uses 2-4 threads | see MET-7, MET-8 |")
w("| Accuracy of the 2D Fourier-Spectral-PINN and of the 2D vanilla PINN against a reference solution | never computed or stored; only S(t) and training times (fourier 175.8 s, vanilla 1461.2 s in `experiment_2d.json`, not quoted) | - |")
w("| Maximum bound principle (abs(u) <= 1; Reviewer 1 minor ii) for the sweeps | no snapshots stored for multi-seed, K, eps, T runs | seed-42 single-run snapshots only (Section 5) |")
w("| Convergence of the inverse runs (trajectories of eps) | sweep JSON stores only the final eps_hat | - |")
w("| Statement in Data Availability (L324): 'minimal dataset sufficient to reproduce Tables 1-6 and Figure 1' | the inverse noise numbers and eps trajectories are missing; the reference solutions are not stored (regenerated by code) | - |")
w("")

# ---------------------------------------------------------------- missing numbers
w("## 4. Numbers that are in the code but NOT in the manuscript (needed by Reviewer 2.2 and for reproducibility)")
w("")
w("- Energy penalty: **lambda_E = 50** (`lam_energy=50.0`), **M = 10 checkpoints** on [0, 1] (spacing 0.1; `n_checkpoints=10`), and 10 T checkpoints in the T-sweep; one-sided squared penalty summed over checkpoints; the violation test uses tolerance 1e-10 at 51 (T=1) or 20 T + 1 points. The text (L125-129) names lambda_E and M without values.")
w("- Other loss weights: Spectral-PINN lam_ic = 10, lam_bc = 1 with the boundary vectors divided by K^2 (so the effective weight on u_x(+-1)^2 is 1/K^4 = 9.5e-7 at K = 32); vanilla lam_ic = lam_bc = 10; 2D: lam_ic = 10, no BC term.")
w("- Optimisation: Adam lr 1e-3 (1D), 2e-3 (inverse, 2D Fourier), cosine annealing; vanilla 8000 epochs, Spectral-PINN 6000 + 1000 warm-up (warm-up not included in the reported times); 100 time-collocation points and a 300-point auxiliary x-grid (Spectral-PINN), 80 x 80 collocation points (vanilla), seeds 42, 123, 7, 2024, 99; SD uses n-1.")
w("- Inverse sweep: eps_init = 2 eps_true, 4000 epochs, seed 42, obs noise seed 1.")
w("")

# ---------------------------------------------------------------- extra facts
w("## 5. Further data-derived facts relevant to the reviewers (from `extra_checks.py`)")
w("")
w("1. **Neumann BC is not satisfied by the trained Spectral-PINN (seed-42 snapshots, refit with K = 32):** |u_x(+-1)| = 0.376/0.381 (t=0.2), 0.241/0.222 (0.5), 0.132/0.106 (0.8), 0.261/0.231 (t=1), versus 0 required; the IC itself has |u_x| = 2. Cause visible in the code: the boundary vectors are divided by K^2 = 1024, so an error of 0.3 in u_x costs only ~1e-7 in the loss.")
w("2. **Maximum bound principle (|u| <= 1), seed-42 single runs at t = 0.2, 0.5, 0.8, 1.0:** Spectral-PINN max|u| = 0.9967, 0.9942, 0.9949, 0.9971 (energy variant identical to 1e-5); vanilla PINN 0.9500, 0.9839, 0.9967, **1.00079** (slightly exceeds 1 at t = 1); FD reference <= 1.000000 at all saved times. Nothing is stored for other seeds, K, eps or T.")
w("3. **Long-horizon energies:** the reference Ginzburg-Landau energy is 0.3346 (t=1), 0.0513 (t=5), 0.0187 (t=10). The penalised network ends at E = 0.2491 (T=5 run, t=5) and 0.2631 (T=10 run, t=10): monotone, but 5x and 14x above the true energy. This supports the 'dissociation' message with a number the manuscript does not give (E[u0] = 0.4004 matches the network's E(0) = 0.4005 and the t=1 values 0.3346 vs 0.3348 agree to 0.07 %).")
w("4. **Wall-clock times are not comparable across tables** (machine load/threads differ): identical configuration K = 64, seed 42, 1000 + 6000 epochs took 268.1 s in Table 1 and 178.4 s in Table 3; K = 32 energy-constrained seed 42 took 162.3 s (Table 1 batch), 189.7 s (single run), 90.6 s (T = 1 row of Table 5); plain K = 32 took 155.7 s (Table 1 batch, seed 42), 178.9 s (single run), 112.0 s (eps sweep). Only within-batch ratios are meaningful (5.43 and 5.40 in the two batches that ran vanilla and Spectral-PINN together).")
w("5. The L2/Linf of the PINN rows depend on the 200-point evaluation grid; the K = 64 seed-42 value is 6.1756e-2 in both Table 1 data and Table 3, as it must be.")
w("")
w("## 6. Out of scope for this audit")
w("")
w("Lean claims (compiles without sorry, only the three standard axioms, Lean 4.33.1): the `lean-toolchain` file says `leanprover/lean4:v4.33.1` (matches); the compile/axiom claims were not re-run here (separate Lean audit). Literature claims and the novelty discussion are not number claims.")
w("")
w("## 7. Files")
w("")
w("- `compute_all.py` / `compute_all_output.txt`: recomputation of Tables 1-7, ratios, statistics, energy counts from `data/*.json`.")
w("- `norm_check.py` / `norm_check_output.txt`: re-execution of the Chebyshev and FD solvers (project code imported read-only) with alternative norms; reference-error study.")
w("- `err_location.py` / `err_location_output.txt`: where the Chebyshev-vs-reference discrepancy sits.")
w("- `extra_checks.py` / `extra_checks_output.txt`: IC/BC compatibility, maximum bound principle, parameter counts, D_c accuracy, energies, K^4 growth.")
w("- `make_report.py`: generates this file (so that the counts above always match the table).")

with open(os.path.join(HERE, "audit_numbers.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
print(dict(cnt), len(R))
print("WRONG:", [r[0] for r in wrong])
print("UNVERIFIABLE:", [r[0] for r in unv])
