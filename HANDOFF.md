# HANDOFF: Symmetry revision round 1 (symmetry-4596052), state of 2026-10-15

**Status: the revision is complete up to the author actions in Section 5.** All experiments have been run (most on ai.atu, the 42 runs listed in `data_r1/runs_workstation.txt` on the i5-9400F workstation), every result of these runs and of the classical baselines, every data table and every figure
is generated from `data_r1/` (2D tables: `data/`), both manuscript versions, the response to the reviewers and the upload packages are built.
(The previous handoff, written before the experiments, is in git history: commit 03649f2.)

## 1. Reviewer points and where they are answered

| # | Reviewer point | Answer (revised manuscript) |
|---|---|---|
| R1-i | Novelty (Xia et al. 2023; Kuetuek-Yuecel 2025) | Cited; claims toned down; contributions restated (Intro, Sec. 4.1) |
| R1-ii | IC violates Neumann | New compatible IC; Appendix D repeats the B1 runs of the as-submitted Spectral-PINN (K=32/64, 3 seeds) and of the classical solver with the original IC (conclusions unchanged; the revised network cannot represent the original IC) |
| R1-iii | Lean not connected to the trained net | A-posteriori energy balance (Sec. 2.5, App. B, `PerturbedFlow.lean`, 39 theorems) + measured for all 120 forward-problem Spectral-PINNs (Sec. 3.4, Table 4) |
| R1-minor | abstract; MBP; S(t) range | Abstract < 200 words; max\|u\| for every network trained on the 1D forward problem (Sec. 3.4); S(t) range corrected (Sec. 3.6) |
| R2-1 | 49 vs 47 | Ratio of seed means + bootstrap CI: 81 [66, 99] on B1 (14 after L-BFGS; on B2 2.9 with Adam, reversed after L-BFGS) |
| R2-2 | lambda_E | lambda_E = 50 stated; ablation 0..5000 (Table 3): no effect |
| R2-3 | Figure 1 curves identical | Penalty practically inactive (positive in 5 of 30 000 epochs; curves differ by 2.2e-6); Fig. 3(b,c) |
| R2-4 | niche vs classical solver | Classical solver 15x (B1) / 84x (B2) more accurate; B3: PINNs fail; Sec. 4.1 states it plainly |
| R2-5 | 38 % bias | Retracted: estimate = 0.70 x eps_init (start value + optimizer step size); variable projection; classical LS; identifiability |

Response letter: `manuscript/response_r1.tex/.pdf` (references resolve against `manuscript_r1_clean.aux` via `xr-hyper`).

## 2. Findings of the runs that changed the drafted text (all incorporated)

1. Rescaling the soft Neumann penalty to unit weight (lambda_BC = K^4) does **not** fix it: L2 0.35 (K=32), 0.40 (K=64), 1.27 (K=128) - ill-conditioned loss. Only the exact basis works.
2. Revised network at fixed budget still degrades with K (8.3e-4, 1.5e-3, 6.7e-3, 3.5e-2 for K=32..256); 2x budget: K=64 2.2e-4, K=128 1.5e-3.
3. B3 (eps=0.05, T=5): all PINNs fail (L2 ~0.9 at t=5, already 5e-2 at t=0.2) - wrong trajectory, not a representation floor.
4. B2: soft penalty harmful already at K=32 (5.1e-2, max|u| 1.30); after L-BFGS the vanilla PINN (1.31e-3) beats the revised Spectral-PINN (1.97e-3), one seed each.
5. Energy penalty: no substantial effect except in training without warm-up, where it makes results worse (4.1e-2 vs 6.5e-3); it is active on B3 (65-799 epochs per run) but does not change the failure there, and for the as-submitted K=64 network (13-28 epochs per run) it lowers the error by 4.5 % (9.49e-2 vs 9.94e-2, same seeds), less than the seed-to-seed scatter. All 120 Spectral-PINNs trained for the forward problem (and all vanilla PINNs on B1/B2) have monotone energy, even with L2 ~ 1; only the vanilla PINN on B3 has increases (169 of 250 intervals).
6. MBP: revised nets max|u| <= 1.0015 on B1 (1.0082 including the penalty regimes, 1.0068 on B3), vanilla up to 1.0084, inaccurate nets up to 2.05.
7. Inverse: joint Adam displacement is set by the lr of log eps (multipliers 1/10/100 -> |log| 0.10/0.38/2.5), ~99 % of it in the first 200 epochs. Variable projection: 0.62 %/0.33 % noise-free, but 11 %/3 % low with sigma=0.01, collapses with sigma=0.05, fails with K=64 or 10 obs. points; classical LS on the same noise: 1.8 %/0.6 % (Cramer-Rao limited).
8. ||r||^2 (degree 6K) is integrated by Clenshaw-Curtis on 4K+9 nodes (order 4K+8) to 1e-12 relative for the seed-42 network of each configuration of Table 4 (checked against 6K+9 nodes, `data_r1/check_quad.json`).

## 2b. Parallel commit ae7911f

While the experiments ran, a second session on the writing PC pushed ae7911f (generators, drafts, response draft) without the new data. It was merged with the data-based versions kept for every overlapping file. Kept from it: the B3 entry of `data_r1/classical.json` (`code/r1_classical.py`), `data_r1/ident.json` (pilot Cramer-Rao study), the `\TODO` macro in the preamble, and `submission/response_to_reviewers_round1.md` (marked SUPERSEDED). The generator uses `data_r1/classical2_*.json` for B3, the classical inverse estimator and the timings; these two extra files are not used.

## 3. Code and data map (additions of this session)

* `code/r1_jobs.py`: new groups `inverse_abl` (joint-Adam ablation, 12 jobs) and `regimes_soft64` (as-submitted K=64 + penalty); both are part of `all`, as is `var3264` (the K=32/64 implementation variants of Table 2). `all` has 200 jobs and reproduces every run of the article except the timings (group `timing`) and the excluded smoke run `smoke_K16`. Run only the groups `all` and `timing`: the other groups (e.g. `variants`, `forward`, `ablation_lamE`, `penalty_regimes`, `stress_K`, `stress_T`) also contain jobs that are not used, and `make_tex_r1.py` would pick up their records (Table 2, energy and balance counts, Appendix D).
* `code/r1_classical2.py b3|inverse|floor|timing` -> `data_r1/classical2_*.json` (classical solver on B3, classical LS inverse with the *same* noise realisations as the neural runs + 20 extra, polynomial representation floors, server timing of the classical solver).
* `code/r1_profiles.py` -> `data_r1/profiles_b3.json` (Fig. 1, B3 profiles); `code/r1_check_quad.py` -> `data_r1/check_quad.json`; `code/finalize_r1.sh` (server post-processing chain).
* `code/make_tex_r1.py`: all `\N{...}` keys and the tables `tab_forward, tab_K, tab_energy, tab_balance, tab_inverse (3 tables), tab_2d (2 tables), tab_orig`; prints keys used in the text but missing in the data. Values use `\ensuremath`, so `\N{}` works in and out of math mode.
* `code/make_figures_r1.py`: `fig_energy, fig_K, fig_balance, fig_inverse (3 panels), fig_b3`.
* `manuscript/build_r1.py` (nested `\input` now flattened), `manuscript/package_r1.py` (upload zips), `manuscript/r1_src/73_app_orig.tex` (Appendix D).
* `data_r1/runs_workstation.txt`: the 42 runs computed on the i5-9400F workstation (all others on ai.atu, Xeon Gold 6430).
* Lean: `lean/AllenCahnVerification/AllenCahnVerification/PerturbedFlow.lean` added and imported by the root file.

## 4. Reproduce

```bash
# on ai.atu (repo is private: copy it with `git bundle create spinn.bundle --all` + scp + `git clone spinn.bundle`)
cd SPINN_Allen_Khan/code
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=~/spectral_pinn_test/venv/bin/python          # torch 2.14 cpu; no system torch on ai.atu
nohup $PY r1_run.py all 36 > ../data_r1/all.log 2>&1 &   # finished jobs are skipped; ~95 min wall on 36 processes
./finalize_r1.sh                                 # waits for "all", then profiles, timing (1 process), classical timing, quadrature check
$PY r1_classical.py                              # classical solver + reaction-only ODE -> data_r1/classical.json (needs torch)
$PY r1_classical2.py b3                          # classical solver + ODE on B3 -> data_r1/classical2_b3.json (needs torch)
# anywhere with numpy/scipy/matplotlib (torch not needed):
python r1_classical2.py inverse; python r1_classical2.py floor
python make_tex_r1.py && python make_figures_r1.py
cd ../manuscript && python build_r1.py both && pdflatex response_r1.tex && pdflatex response_r1.tex && python package_r1.py
```

## 5. Remaining author actions

1. Read `manuscript/manuscript_r1_highlighted.pdf` and `manuscript/response_r1.pdf`; all co-authors approve. The response is signed by the corresponding author.
2. Check the GenAI statement in the Acknowledgments (rewritten to describe the actual AI assistance in the revision).
3. Lean: the parallel session on the writing PC reported `lake build` OK with the pinned Lean 4.33.1 project including `PerturbedFlow.lean` (commit ae7911f, 8710 jobs). If you want to repeat it, run `lake build` in `lean/AllenCahnVerification` on a machine with a fast GitHub connection. `PerturbedFlow.lean`, `Basic.lean` and `GradientFlow.lean` were re-checked here only against Lean 4.35.0-rc2 + Mathlib `b0df0c4` (no errors, no warnings, standard axioms); cloning Mathlib for the pinned 4.33.1 build failed (GitHub git transfer ~50-100 KB/s, then disconnects) on both the PC and ai.atu.
4. Upload at https://susy.mdpi.com/user/manuscripts/upload?pre_hash_key=760fc201bdbd8cda5fc7fe42d548ad3d&resubmit=1 (keep the manuscript ID): highlighted + clean PDF, response letter, `submission/symmetry-4596052_R1_manuscript.zip`, Supporting Information `submission/symmetry-4596052_R1_supporting_data.zip`.

## 6. Error audit (2026-10-13/15)

A multi-agent audit of the rebuilt manuscript (14 independent reviewers by section and by error type, every finding checked against the sources, data and code by a skeptical verifier, then a second round on the applied fixes) found 280 raw findings = 120 unique issues, of which 110 were confirmed and fixed, plus 31 in the check round (several introduced by the first fixes). Examples of what was wrong: |u_0|<=1/8 "in the interior" (false; u_0 reaches -1 at the walls); the IC term of the loss is lambda_IC/(K+1)|c(0)-c_0|^2 in the code; "the penalty is inactive" although it is active on B3 and for soft K=64 (it has no substantial effect there); counts of runs and the "revised" set in the maximum-bound sentence; the Lean claim (only the abstract form is machine-checked); several misattributed citations (Raissi Allen-Cahn example is discrete-time, S-PINNs are equivariant not penalty-based, Invariant PINNs, SAV genealogy, Tsachouridis); the ICCMA abstract has five authors (Nauryzbay first, Ertai last); a TAB that had eaten `\t` of `\texttt`; lambda_E=50 stated as default though the default is no penalty; British/American spelling mix (now American). A third round (same day) fixed 76 further issues and 12 more found by its check round, among them the Data S1 and data-availability references to the 2D results (Section 3.6, not 2.7), an orphaned "Abbreviations" heading, a stale `c' = -G + r` in the `PerturbedFlow.lean` header, the README's job-group and script notes, and upload ZIPs that predated the fixes. Audit material: `%TEMP%`-only, not in the repo; the fixes are in the git history.

A fourth round (2026-10-14; 19 reviewers, 118 raw findings = 100 issues, each checked by two independent verifiers with a tie-breaker) confirmed 71, and its check round 12 more. Besides wording, numbers and captions (e.g. variable projection is accurate on noise-free data at eps_true = 0.05 and 0.1 but 22 % low at 0.01; the a-posteriori bound includes the boundary defect; Liu et al. 2024 misdescribed; Table 6/7 layout and captions), it found that the highlighted version marked only the 'principal changes': the new title, many corrected sentences, every table and all new references were black. `manuscript/build_r1.py` now builds the highlighted version mechanically: everything is blue except paragraphs, headings, title-page data and references that are identical to the first version (`manuscript/v1_submitted/manuscript.tex`, needed for the build); `\rev{}` no longer matters for the colours. The response letter describes this.

A fifth round (2026-10-15; 22 reviewers including two that re-read the paper as the original referees and one that recompiled the Lean files; 68 raw findings = 52 issues, 39 confirmed) found, besides wording in code comments, docs and appendices: the Neumann basis phi_k = T_k - k^2/(k+2)^2 T_{k+2} is the classical Chebyshev-Galerkin basis of Shen (1995), which had no citation and was listed as a contribution (now cited, and the contribution is its use as the output basis); the first-version reference had first-order (not inconsistent) wall rows; Table A1 lacked the 1/4 ||r||^2 bound; tables were cited out of order; `code/stress_time_sweep.py` no longer passed the horizon T (fixed). The GenAI statement now also mentions the analysis of the results, as the first version did; the authors should confirm it.

Points for the authors to confirm: the ICCMA 2026 abstract author list as printed in the Abstract Book (Nauryzbay, Abdiramanov, Zhumazhanov, Urnek, Ertai; Issa is not an author of the abstract, Ertai is not an author of the paper), and the model versions in the GenAI statement (only the models used for the revision are named).

Machine notes: everything is float64 with tiny networks, CPU only. ai.atu has 32 physical cores / 64 hardware threads and is shared (another user's job uses ~24 threads); 36 processes was a good size. The Windows bash here mangles backslashes in heredocs: write LaTeX/Python with an editor (or build strings with `chr(92)`), not `cat <<EOF`.

## 6. Update 2026-10-19: Lean removed from the paper, text condensed (37 -> 27 pages)

* User decision: no Lean in the manuscript, pure mathematics instead. The energy balance is Theorem 2 with a four-step proof in Section 2.5 (`manuscript/r1_src/20_methods.tex`); the former Appendices A (Lean basics) and B (perturbed flow, Lean table) are deleted and replaced by a short Appendix A "Energy Balance: Evaluation and Scope" (quadrature exactness, finite-difference check, what the theorem does not say). Appendices are now A balance, B symmetry diagnostic, C original IC, D certified reference. Lean mentions removed from keywords, abstract, Introduction, Methods (old Section 2.8), Discussion (old Section 4.3 title), Conclusions, Data Availability and the GenAI statement. The Lean project (`lean/`, `revision1/lean/`) stays in the repository but is not referenced by the paper.
* Introduction, Methods, Results, Discussion were condensed without changing any number or conclusion (all `\N{...}` keys are unchanged). New related work: Yu, Qi, Oseledets, Chen, Spectral informed neural networks, J. Comput. Appl. Math. 477 (2026) 117178 (bib key `ref_Yu2026SINN`). `references/3138_Spectral_Informed_Neural__Originally Submitted PDF.pdf` is an anonymous, unpublished ICML submission (Modified SINNs, high-dimensional PDEs): deliberately NOT cited. The OpenReview page is behind a browser check that must not be bypassed.
* `manuscript/response_r1.tex` updated (Lean statements replaced, label `app:balance`, SINN mentioned); `\rev` is defined in the letter because the manuscript .aux carries captions with `\rev`.
* Rebuild: `python manuscript/build_r1.py both`, then `pdflatex response_r1.tex` twice in `manuscript/`, then `python manuscript/package_r1.py`.

## 7. Update 2026-10-19 (later): short sentences

* User request: the text was hard to read. Introduction, Methods, Results, Discussion, Conclusions, abstract, appendices, all table captions (generator `code/make_tex_r1.py`) and the response letter were rewritten with short sentences (mean sentence length of the body 49 -> 20 words, 62% -> 12% of sentences over 30 words). No number, claim or `\N{...}` key changed. Manuscript: 26 pages. The forward table floats (`place="t"`) to avoid a blank half page.
* Style rule for further edits: one idea per sentence, no nested parentheses, details into tables or appendices.
