# Spectral-PINN for the Allen–Cahn Equation

Research code, experimental data, manuscript, and formal proofs for:

> **A Spectral–Physics-Informed Neural Network Hybrid with Energy-Constrained
> Training for the Allen–Cahn Equation: Forward and Inverse Problems**
> Zhuldyz Issa, Adil Nauryzbay, Zhanars Abdiramanov (corresponding), Esimkhan
> Zhumazhanov, Dastan Urnek
> Submitted to *Symmetry* (MDPI), Special Issue "Selected Papers from the 4th
> International Conference: Constructive Mathematical Analysis" (ICCMA 2026).

The paper introduces a **Spectral-PINN**: instead of a network that outputs
the field value directly, the network outputs time-dependent Chebyshev (1D)
or Fourier (2D) spectral coefficients, and spatial derivatives are computed
exactly via fixed linear operators in coefficient space, avoiding spatial
automatic differentiation. Training adds a one-sided **energy-constrained**
penalty that enforces discrete Ginzburg–Landau energy dissipation. The method
is evaluated on the forward problem, an inverse problem (recovering the
interface parameter ε), a 2D Fourier variant with architecturally exact
reflection symmetry, and stress-tested at scale (large mode counts, small ε,
long time horizons, larger 2D grids). The continuous energy-dissipation
theorem underlying the method is additionally verified formally in
**Lean 4 / Mathlib**.

## Repository layout

```
code/          Experiment implementation (PyTorch)
  common.py                  shared utilities (grids, IC, PDE residual, energy)
  chebyshev_spectral.py      classical Chebyshev spectral-collocation baseline
  vanilla_pinn.py            standard PINN baseline
  spectral_pinn.py           Spectral-PINN (Chebyshev-coefficient network)
  inverse_problem.py         inverse recovery of epsilon
  experiment_2d_symmetry.py  2D Fourier-Spectral-PINN, symmetry diagnostics
  stress_*.py                large-scale stress tests (K, epsilon, T, 2D)
  server_sweep_*.py          multiprocessing sweep drivers for server runs
  run_experiments.py         entry point tying the above together
  make_figures.py            regenerates the paper's figures from data/

data/          Raw JSON results for every experiment reported in the paper
figures/       Generated figures (PNG)

manuscript/    LaTeX source (MDPI Symmetry template)
  manuscript.tex             the manuscript
  manuscript.pdf             compiled PDF
  Definitions/               MDPI class/style files (mdpi.cls, logos, etc.)

lean/AllenCahnVerification/  Lean 4 + Mathlib formal verification
  AllenCahnVerification/Basic.lean         polarization identity, discrete
                                            Gronwall step, step-size existence
  AllenCahnVerification/GradientFlow.lean  abstract Hilbert-space gradient-flow
                                            energy-dissipation theorem

submission/    Cover-letter / submission-checklist materials for MDPI
```

## Reproducing the experiments

```bash
cd code
pip install torch numpy scipy matplotlib
python run_experiments.py        # forward, inverse, 2D experiments -> data/
python make_figures.py           # data/*.json -> figures/*.png
```

The `stress_*.py` and `server_sweep_*.py` scripts run the large-scale sweeps
(up to K=256 Chebyshev modes) and are intended for a multi-core machine; they
set `OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS`/`OMP_NUM_THREADS` explicitly to
avoid BLAS thread oversubscription under `multiprocessing`.

## Rebuilding the manuscript

```bash
cd manuscript
pdflatex -interaction=nonstopmode manuscript.tex
pdflatex -interaction=nonstopmode manuscript.tex   # twice, for cross-references
```

Requires a TeX distribution with the MDPI `symmetry` document class (provided
in `Definitions/`, or install via `mdpi.cls` on CTAN).

## Checking the Lean proofs

```bash
cd lean/AllenCahnVerification
lake build
```

Requires [`elan`](https://github.com/leanprover/elan) (Lean 4 toolchain
manager); the exact toolchain version is pinned in `lean-toolchain`. The
`.lake/` build cache (Mathlib source and compiled `.olean` files, several GB)
is not tracked in this repository and is fetched/built automatically by
`lake build`. Verify no `sorry` was used with:

```bash
lake env lean --run <(echo '#print axioms AllenCahnVerification.gradientFlow_deriv_energy_nonpos')
```

or open the relevant theorem in an editor with the Lean 4 extension and check
via `#print axioms <theorem_name>` — it should depend only on `propext`,
`Classical.choice`, and `Quot.sound`.
