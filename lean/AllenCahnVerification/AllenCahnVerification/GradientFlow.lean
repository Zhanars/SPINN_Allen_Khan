import Mathlib

set_option linter.style.header false

/-!
# Theorem 1 (energy dissipation) as an abstract gradient-flow theorem

Theorem 1 of the manuscript states that a smooth solution of the Allen-Cahn
equation `u_t = eps^2 * u_xx - f(u)` with homogeneous Neumann boundary
conditions satisfies `dE/dt = -‖u_t‖^2_{L^2} <= 0`, where `E` is the
Ginzburg-Landau free energy. The paper's proof (integration by parts) shows
that the right-hand side of the PDE is exactly the negative *L^2 gradient* of
`E`, i.e. the Allen-Cahn equation is the `L^2`-gradient flow `u_t = -grad E(u)`
of `E`. This file formalizes, completely and for a general real Hilbert space
`F` (the instantiation `F = L^2(Omega)` for the PDE is only formal, because the Ginzburg-Landau
energy is not Fréchet differentiable on `L^2`; see Remark A1(iii) and Theorem A3 of the
article), the theorem that
*any* gradient flow dissipates its energy at exactly the rate `-‖u_t‖^2`:
this is the actual ODE/calculus content of Theorem 1, proved here via
Mathlib's `HasGradientAt` / `HasFDerivAt` / chain-rule machinery, rather than
only the finite-dimensional real-arithmetic lemmas of `Basic.lean`.

Unlike `Basic.lean`, this file is not about a numerical-scheme inequality: it
is a genuine derivative computation `d/dt E(c(t)) = ...`, i.e. the calculus
content of Theorem 1. The one part
of Theorem 1 *not* re-derived here is the identification, for the specific
Ginzburg-Landau functional `E[u] = int (eps^2/2 |u_x|^2 + F(u)) dx`, of its
`L^2` gradient with `-eps^2 u_xx + f(u)` (a classical calculus-of-variations
computation via integration by parts, already given in the manuscript's proof
of Theorem 1 and requiring Sobolev-space machinery beyond what is needed for
the abstract statement here).
-/

open scoped Gradient
open InnerProductSpace

variable {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F] [CompleteSpace F]

/-- **Gradient-flow energy dissipation** (abstract form of Theorem 1).
If `E : F -> R` has gradient `gradient E (c t)` at `c t`, and the path `c`
solves the gradient-flow ODE `c'(t) = -gradient E (c t)` at time `t`, then
`s -> E (c s)` has derivative `-‖c'(t)‖^2` at `t`. -/
theorem gradientFlow_hasDerivAt_energy
    (E : F → ℝ) (c : ℝ → F) (t : ℝ)
    (hE : HasGradientAt E (gradient E (c t)) (c t))
    (hflow : HasDerivAt c (-(gradient E (c t))) t) :
    HasDerivAt (fun s => E (c s)) (-‖deriv c t‖ ^ 2) t := by
  have hderiv_c : deriv c t = -gradient E (c t) := hflow.deriv
  have hFD : HasFDerivAt E (toDual ℝ F (gradient E (c t))) (c t) := hE.hasFDerivAt
  have hcomp : HasDerivAt (E ∘ c)
      (toDual ℝ F (gradient E (c t)) (-(gradient E (c t)))) t :=
    HasFDerivAt.comp_hasDerivAt_of_eq t hFD hflow rfl
  have hval : toDual ℝ F (gradient E (c t)) (-(gradient E (c t))) = -‖deriv c t‖ ^ 2 := by
    rw [toDual_apply_apply, inner_neg_right, real_inner_self_eq_norm_sq, hderiv_c, norm_neg]
  rw [hval] at hcomp
  exact hcomp

/-- **Corollary: the energy is instantaneously non-increasing along the flow.**
This is Theorem 1's conclusion `dE/dt <= 0`, in the abstract gradient-flow
setting: with `F = L^2(Omega)`, `E` the Ginzburg-Landau functional, and the
gradient-flow hypothesis instantiated via the manuscript's calculation
`grad E(u) = -eps^2 u_xx + f(u)`, this formally recovers Theorem 1 (Corollary A1 of the article). -/
theorem gradientFlow_deriv_energy_nonpos
    (E : F → ℝ) (c : ℝ → F) (t : ℝ)
    (hE : HasGradientAt E (gradient E (c t)) (c t))
    (hflow : HasDerivAt c (-(gradient E (c t))) t) :
    deriv (fun s => E (c s)) t ≤ 0 := by
  rw [(gradientFlow_hasDerivAt_energy E c t hE hflow).deriv]
  exact neg_nonpos.mpr (sq_nonneg _)

/-- **Strict form**: the energy strictly decreases unless the flow has already
reached a critical point (`c'(t) = 0`), matching the manuscript's remark that
equality `dE/dt = 0` holds only where `u_t = 0`. -/
theorem gradientFlow_deriv_energy_eq_zero_iff
    (E : F → ℝ) (c : ℝ → F) (t : ℝ)
    (hE : HasGradientAt E (gradient E (c t)) (c t))
    (hflow : HasDerivAt c (-(gradient E (c t))) t) :
    deriv (fun s => E (c s)) t = 0 ↔ deriv c t = 0 := by
  rw [(gradientFlow_hasDerivAt_energy E c t hE hflow).deriv]
  constructor
  · intro h
    have : ‖deriv c t‖ ^ 2 = 0 := by linarith
    have hnorm : ‖deriv c t‖ = 0 := by nlinarith [norm_nonneg (deriv c t), sq_nonneg (‖deriv c t‖)]
    exact norm_eq_zero.mp hnorm
  · intro h
    simp [h]
