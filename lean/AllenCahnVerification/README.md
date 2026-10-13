# AllenCahnVerification

Lean 4 / Mathlib proofs for the article "A Chebyshev-Coefficient Physics-Informed Neural Network for the Allen–Cahn Equation"
(Symmetry, symmetry-4596052):

- `AllenCahnVerification/Basic.lean`: polarization identity, discrete Gronwall step, existence of an admissible step size (Appendix A).
- `AllenCahnVerification/GradientFlow.lean`: energy dissipation along an exact gradient flow in a real Hilbert space (Appendix A).
- `AllenCahnVerification/PerturbedFlow.lean`: a-posteriori energy balance for an arbitrary differentiable trajectory with residual and
  boundary defect, integrated and discrete-time bounds, integration by parts for polynomials (Appendix B; 39 theorems).

Build with `lake build` (toolchain pinned in `lean-toolchain`).
