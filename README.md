# Chebyshev-coefficient PINN for the Allen–Cahn equation (Symmetry, symmetry-4596052)

Research code, raw results, manuscript and Lean 4 proofs for:

> **A Chebyshev-Coefficient Physics-Informed Neural Network for the Allen–Cahn Equation**
> Zhuldyz Issa, Adil Nauryzbay, Zhanars Abdiramanov (corresponding), Yessimkhan Zhumazhanov, Dastan Urnek.
> *Symmetry* (MDPI), manuscript symmetry-4596052, revision 1 (resubmission). Extended version of an ICCMA 2026 conference abstract.

The network outputs time-dependent Chebyshev (1D) or Fourier (2D) coefficients, so spatial derivatives are exact linear maps in
coefficient space. Revision 1:

- evaluates the network on three 1D benchmarks against a certified reference solution, a space–time PINN and a classical Chebyshev
  solver, with an exactly Neumann-compatible basis;
- tests the optional one-sided energy penalty (every 1D forward-problem Spectral-PINN has monotone energy without it, so the penalty
  is unnecessary; the maximum bound is checked separately);
- adds an a-posteriori energy balance that bounds the energy growth of any trained trajectory by its measured residual and boundary defect (abstract form
  machine-checked in Lean 4);
- revisits the inverse problem for ε (joint Adam, variable projection, classical least squares, identifiability).

The 2D reflection-symmetric Fourier results are those of the first version (`run_experiments.py 2d`, `server_sweep_2d.py`,
`stress_2d_sweep.py`; Tables 8 and 9). The 1D stress tests and sweeps of the first version (`stress_k_sweep.py`,
`stress_eps_sweep.py`, `stress_time_sweep.py`, `server_sweep_forward.py`, `server_sweep_inverse.py`) are superseded and not used
in revision 1. `HANDOFF.md` records the state of the revision; where its counts differ, the manuscript and
`data_r1/numbers.json` are authoritative.

## Repository layout

```
code/                              PyTorch, double precision
  common.py, chebyshev_spectral.py, vanilla_pinn.py, spectral_pinn.py
                                   shared utilities, classical solver, space-time PINN, Spectral-PINN
  refsol.py                        certified reference solutions (Chebyshev collocation + Radau, N-refinement)
  r1_core.py                       evaluation protocol of revision 1 (200-point grid, max-bound check, energy, a-posteriori balance)
  r1_inverse.py                    inverse problem (joint Adam, variable projection, data-fitting phase)
  r1_jobs.py, r1_run.py            job groups;  python r1_run.py <group> [nproc] -> data_r1/runs/<tag>.json, data_r1/models/<tag>.pt
  r1_classical.py                  classical baselines -> data_r1/classical.json
  r1_classical2.py                 python r1_classical2.py {b3|inverse|floor|timing} -> data_r1/classical2_<mode>.json
  r1_profiles.py, r1_check_quad.py B3 profiles, quadrature check
  check_identity.py                finite-difference sanity check of the energy identity (Appendix B)
  r1_tables.py                     aggregation -> data_r1/summary.json
  make_tex_r1.py                   the revision-1 results of the training runs and classical baselines, and the data tables of the manuscript -> manuscript/r1_src/numbers_r1.tex, tab_*.tex, data_r1/numbers.json
  make_figures_r1.py               figures -> manuscript/figures_r1/*.png
  finalize_r1.sh                   server post-processing chain (timing runs, profiles, quadrature check)
  experiment_2d_symmetry.py        2D Fourier network and symmetry diagnostic (first version)
  run_experiments.py, inverse_problem.py, make_figures.py, stress_*.py, server_sweep_*.py
                                   first-version pipeline, kept for reference (the scripts write to data/,
                                   make_figures.py to figures/)
data_r1/       raw results of revision 1: runs/ (one JSON per training run), models/ (weights of the 1D forward networks),
               ref/ (certified reference solutions), classical*.json, check_quad.json, profiles_b3.json,
               runs_workstation.txt, numbers.json
data/          raw results of the first version; experiment_2d.json, symmetry_2d_multiseed.json and stress_2d_sweep.json
               feed the 2D tables of revision 1
manuscript/    r1_src/ (LaTeX sources), build_r1.py, make_docx_mdpi.py, package_r1.py, manuscript_r1_mdpi.docx (the article, MDPI Word
               template), response_r1.{tex,pdf}, Definitions/ (MDPI class files); build/ = LaTeX intermediates (git-ignored)
lean/AllenCahnVerification/AllenCahnVerification/
  Basic.lean          polarization identity, discrete Gronwall step, step-size existence (Appendix A)
  GradientFlow.lean   abstract gradient-flow energy-dissipation theorem (Appendix A)
  PerturbedFlow.lean  a-posteriori energy balance for an arbitrary trajectory (Appendix B)
submission/    upload packages and submission notes
```

## Reproducing revision 1

```bash
cd code
pip install torch numpy scipy matplotlib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python r1_run.py all 36                  # every training run used in revision 1 except the timings; finished jobs are skipped
python r1_run.py timing 1                # clean timings, one job at a time
python r1_classical.py
python r1_classical2.py b3               # likewise: inverse, floor, timing
python r1_profiles.py
python r1_check_quad.py
python make_tex_r1.py                    # numbers and tables of the manuscript
python make_figures_r1.py                # figures of the manuscript
cd ../manuscript
python build_r1.py                       # LaTeX intermediates in manuscript/build/ (aux file with the numbers)
python make_docx_mdpi.py                 # manuscript_r1_mdpi.docx (the single version of the article)
pdflatex response_r1.tex && pdflatex response_r1.tex
python package_r1.py                     # upload packages in submission/
```

Do not run `run_experiments.py 2d`, `server_sweep_2d.py` or `stress_2d_sweep.py`: they overwrite `data/experiment_2d.json`,
`data/symmetry_2d_multiseed.json` and `data/stress_2d_sweep.json`, which `make_tex_r1.py` reads for Tables 8 and 9 (the other
first-version scripts overwrite other files in `data/`). Run only the job groups `all` and `timing`: the other groups in
`r1_jobs.py` (e.g. `variants`, `forward`, `ablation_lamE`, `penalty_regimes`, `stress_K`, `stress_T`) also contain jobs that
are not used in the article, and `make_tex_r1.py` would pick up their records.

## Checking the Lean proofs

```bash
cd lean/AllenCahnVerification
lake build
```

Requires [`elan`](https://github.com/leanprover/elan); the toolchain version is pinned in `lean-toolchain`. The `.lake/` build cache
(Mathlib source and compiled `.olean` files, several GB) is not tracked and is fetched/built by `lake build`. `PerturbedFlow.lean` ends
with `#print axioms` for each of its 39 theorems; the output must list only `propext`, `Classical.choice` and `Quot.sound`. For the
other files, use `#print axioms <theorem_name>` in an editor with the Lean 4 extension.
