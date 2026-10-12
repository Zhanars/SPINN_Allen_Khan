# HANDOFF: Symmetry revision round 1 (symmetry-4596052), state of 2026-10-12

Work stopped here on purpose: the experiments must be run on the server **ai.atu** (64 cores), which is not reachable from the
machine where the code was written. Everything needed is in this repository. Continue on the PC that can reach ai.atu.

## 1. Decision and what the reviewers said

MDPI *Symmetry* rejected the manuscript "to allow revision" (resubmit with the same ID, highlighted changes, point-by-point response).
Reviewer PDFs: `manuscript/round 1/`. Submitted version (for reference): `manuscript/v1_submitted/`.

| # | Reviewer point | Status |
|---|---|---|
| R1-i | Novelty: Xia-Boettcher-Chou 2023 (coefficient networks, Chebyshev) and Kuetuek-Yuecel 2025 (energy penalty) | Accepted. Verified, cited, claims toned down (Intro, drafted). Lit review: `revision1/lit/lit_review.md` |
| R1-ii | Initial condition x^2 cos(pi x) violates Neumann (u_x(-1,0)=2) | **Confirmed.** New IC u0=(cos(pi x)-cos^2(pi x))/2. Pilot: no change in conclusions of the submitted code |
| R1-iii | Lean theorem not connected to the trained network | Answered: a-posteriori energy balance dE/dt=-|u_t|^2+<r,u_t>+eps^2[u_x u_t], Lean-verified (`revision1/lean/PerturbedFlow.lean`, 41 theorems, standard axioms only, compiles in ~100 s) + residual measured for every trained net (`balance` in the run JSONs) |
| R1-minor | abstract too long; MBP discussion; S(t) range wrong | Abstract: TODO (<=200 words). MBP: `max_abs_u` recorded in every run. S(t): text fixed in 30_results.tex draft (true range 1.3e-31..3.7e-31) |
| R2-1 | 49 vs 47 | 49 = ratio of 5-seed means; 47 = seed 42 only. Bootstrap CI code in `code/r1_tables.py` |
| R2-2 | lambda_E not reported | lambda_E = 50 (10 checkpoints per unit time). Ablation lambda_E in {0,1,10,50,500,5000}: group `abl_hard` |
| R2-3 | Figure 1 curves identical | Penalty is inactive (audit: nonzero in 1 of 120 logged epochs; curves differ by 4.5e-7). Group `regimes_hard` looks for regimes where it matters |
| R2-4 | niche vs classical Chebyshev solver | Honest answer drafted in `manuscript/r1_src/40_discussion.tex` (classical solver wins the forward problem; benchmarks B2, B3 added) |
| R2-5 | 38 % inverse bias unexplained | **Retracted.** It is an optimisation artefact (see 2) |

## 2. Findings that change the paper (all from our own audit/runs)

1. **"Accuracy degrades with K" was an implementation artefact.** The Neumann penalty was divided by K^2 and then squared (effective weight lambda_BC/K^4); at K=64 the wall slope of the trained network is ~4 and L2 6e-2 (seed 42, pilot). With the exact Neumann basis (`hard_bc=True`, phi_k = T_k - k^2/(k+2)^2 T_{k+2}) K=64 gives L2 1.7e-3 vs 9.2e-4 at K=32. Gauss-Lobatto collocation (`grid="cgl"`) is a secondary improvement.
2. **The ~38 % inverse "bias" is the start value.** Joint Adam gives eps_hat ~ 0.68-0.72 x eps_init for every init and eps_true (data weight, K, snapshots irrelevant). Classical least squares recovers eps exactly (noise-free), Fisher-limited otherwise; eps=0.01 is not identifiable for sigma>=0.01. A working neural estimator: variable projection (closed-form eps^2) + data-fit phase + large data weight: pilot (hard basis, lambda_data=2000): eps_hat/eps_true = 0.97 (eps=0.05), 0.96 (eps=0.1).
3. **Energy penalty is inactive**; unpenalised nets (also the vanilla PINN) already have 0 energy increases. MBP: networks overshoot |u|<=1 by up to ~1e-2 (not guaranteed).
4. **Benchmark B1 (eps=0.01,T=1) is reaction dominated**: the closed-form diffusion-free ODE already has L2 1.9e-3 (compat IC), as good as the PINN. Added B2 (eps=0.1,T=1: ODE error 0.147) and B3 (eps=0.05,T=5, interfaces; K=32 representation floor 1.2e-2, K=64 floor 1.9e-4).
5. Reference solver of v1 had first-order wall rows; Chebyshev-baseline row used a wrong norm; 22 of 77 numerical claims of v1 were wrong/unverifiable (`revision1/audit_numbers/audit_numbers.md`).
6. Certified reference: `code/refsol.py` (Chebyshev + Radau, certified <=1e-8; report `revision1/audit_code/audit_code.md`).

## 3. Code map (all new code in `code/`)

* `r1_core.py` evaluation (errors on the 200-pt grid, MBP, energy violations at 3 resolutions, energy balance diagnostics), `refsol.py` reference.
* `r1_run.py <group> [nproc]` runs a job group; each job writes `data_r1/runs/<tag>.json` (+ `data_r1/models/<tag>.pt`); finished tags are skipped.
* `r1_jobs.py` job groups. **`all`** = everything still needed (179 jobs, 35 done, ~960 core-minutes), longest first. Other groups: `forward`, `forward_hard`, `b2`, `abl_hard`, `regimes_hard`, `stress_K_hard`, `variants2`, `lbfgs`, `b3`, `inverse`, `timing`.
* `r1_inverse.py` inverse problem (modes `adam`, `varpro`).
* `r1_tables.py` aggregation (mean+-SD, bootstrap ratios), `make_tex_r1.py` (**only the `ode.*`/`cheb.*` number macros are implemented; the table/number generators for the other `\N{...}` keys used in `manuscript/r1_src/30_results.tex` still have to be written**), `make_figures_r1.py` (4 figures, first draft).
* `spectral_pinn.py` options: `grid="uniform"|"cgl"`, `hard_bc`, `lbfgs_iters`, `T`, `u0`. `vanilla_pinn.py`: `T`, `u0`, `lbfgs_iters`. Defaults reproduce the submitted code.
* Old pilot/diagnostic scripts: `revision1/pilot/`. Audits: `revision1/audit_*`, Lean: `revision1/lean/`, literature: `revision1/lit/`.

## 4. What to run on ai.atu

```bash
git clone https://github.com/Zhanars/SPINN_Allen_Khan.git && cd SPINN_Allen_Khan/code
pip install torch numpy scipy matplotlib          # CPU torch is enough; Python 3.10+
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1   # IMPORTANT (thread oversubscription bug seen before)
nohup python r1_run.py all 56 > ../data_r1/all.log 2>&1 &     # ~1000 core-minutes: expect 30-60 min wall (longest single job ~50 min)
# when finished and the machine is idle (clean timings for the tables):
python r1_run.py timing 1
```

Before trusting the inverse results run first `inv_vp_e0.05_ld2000` and check `ratio_to_true` is about 0.9-1.0: the least-squares eps^2 inside `r1_inverse.py` was just changed to use the same quadrature weights as the PDE loss (the pilot used unweighted sums consistently; the first harness version mixed them and gave 0.1-0.2). If it is far off, revert to unweighted sums for both loss and LS.

Not yet run at all: `variants2` (soft penalty with lambda_BC=K^4 control), `lbfgs` (stronger PINN baselines), `b3`, `abl_hard`, `regimes_hard`, `stress_K_hard`, `inverse` (only the 10 joint-Adam init-sweep runs exist), `b2`. Already finished (in `data_r1/runs`, 42 JSONs): `fh_K32`, `fh_K32_E50`, `fh_K64` (5 seeds each), `fh_K128` (s42, s7), `fwd_compat_vanilla` s42/123/7/99, `fwd_compat_K64` s42/123/7/2024, `var_K32/K64_*` pilot variants (single seed), `inv_adam_*` (10 runs, eps_true=0.05 only), `smoke`. Preliminary numbers from these (L2 at t=1, B1): vanilla ~7e-2, revised K=32 ~9.8e-4, revised K=64 ~1.5e-3, as-submitted K=64 ~9.6e-2. Their `train_time` was measured under CPU oversubscription: use the `timing` group for the time columns. Delete a JSON to force a rerun.

`data_r1/classical.json` (Chebyshev baseline + ODE) is final; `python r1_classical.py` regenerates it.

## 4b. Update (same day, after the first commit): what was written afterwards without new experiments

* `code/make_tex_r1.py` is now complete: it writes `manuscript/r1_src/numbers_r1.tex` and all table fragments (`tab_forward/K/energy/balance/inverse/2d/orig.tex`) from `data_r1`; missing data never stop the build (they show as `--` or `[n/a: key]`). Run it after the experiments, then `python manuscript/build_r1.py both` (prints pages, LaTeX errors, undefined references, number of TODO markers and missing numbers; currently 31 pages, 0 errors, 3 TODO + 16 missing numbers because the data are not there yet). `code/make_figures_r1.py` makes 5 figures (energy, K, balance, inverse, B3 profiles; the last needs the saved models).
* Written: Introduction, Methods (all subsections), Results/Discussion drafts (check against data!), Conclusions, abstract (`abstract_r1.tex`), Appendices (Lean A, perturbed flow B, symmetry C, original IC D, certified reference E), back matter (GenAI statement extended: please confirm the wording), `submission/response_to_reviewers_round1.md` (draft with TODOs). Lean: `PerturbedFlow.lean` is now part of the Lean project (`lake build` OK, 8710 jobs).
* `code/r1_classical.py` also computes the classical solver on B3 (stored in `data_r1/classical.json`); `data_r1/ident.json` = classical least-squares / Cramer-Rao study (from `revision1/pilot/classical_inverse.py`).
* The `[TODO: ...]` markers in the manuscript text mark the sentences that can only be written after the data exist (regimes of the energy penalty; B2/B3 sentence in the Conclusions); `[n/a: key]` markers are numbers that will be filled automatically.

## 5. What remains to write (after the runs)

1. Implement the generators for the remaining `\N{...}` keys and the tables `tab_forward.tex, tab_K.tex, tab_energy.tex, tab_balance.tex, tab_inverse.tex, tab_2d.tex` (referenced from `30_results.tex`; the 2D numbers come from the unchanged `data/symmetry_2d_multiseed.json`, `data/stress_2d_sweep.json`, `data/experiment_2d.json`). Re-read every sentence of `30_results.tex` and `40_discussion.tex` against the actual numbers: they were drafted before the data existed and make claims (e.g. "the vanilla PINN overshoots more", "penalty does nothing in all regimes", "variable projection error falls with lambda_data") that **must be checked and rewritten if the data say otherwise**.
2. `50_conclusions.tex` and `abstract_r1.tex` (<=200 words), `73_app_orig.tex` (original-IC replication; placeholder `%%APP_ORIG%%` in `70_appendix.tex`), update `60_backmatter.tex` (Data availability, GenAI statement), Lean appendix text (`70_appendix.tex` contains `revision1/lean/perturbed_flow.tex`; adapt its Remark about the K^-4 penalty to the revised method, copy `PerturbedFlow.lean` into the Lean project and compile).
3. `python manuscript/build_r1.py both` builds `manuscript_r1_highlighted.tex/pdf` (changes in blue via `\rev{}`) and the clean version, reorders the bibliography by first citation (all 31 new references are verified in `revision1/lit/verify_*.md`; add pages 26548--26560 for Krishnapriyan2021).
4. Response to reviewers (point by point, honest: list the three retracted conclusions and the other corrections), update `submission/summary.md`, re-send the PDF/zip, commit/push.
5. Resubmit at https://susy.mdpi.com/user/manuscripts/upload?pre_hash_key=760fc201bdbd8cda5fc7fe42d548ad3d&resubmit=1 (author action; keep manuscript ID).

Machine notes: GPU (RTX 4060) exists on the writing PC but everything is float64 with tiny networks, so it was not used. The authors' Windows shell mangles backslashes in bash heredocs: write LaTeX/Lean/Python files with an editor, not `cat <<EOF`.
