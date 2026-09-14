import Mathlib

set_option linter.style.header false

/-!
# Energy-stability lemmas for the Allen-Cahn Spectral-PINN manuscript

Formal verification of the two real-arithmetic identities/inequalities
underlying the energy-stability arguments used in the manuscript
"A Spectral-PINN Hybrid with Energy-Constrained Training for the
Allen-Cahn Equation" (Section: energy dissipation / discrete stability).

This file is deliberately *not* an attempt to formalize the full PDE
theory (weak derivatives, Sobolev spaces, ...) in Mathlib. Instead it
isolates the finite-dimensional algebraic steps that, in the numerical
PDE literature, are the most common source of silent errors -- e.g. a
peer review of an earlier, related Allen-Cahn manuscript from the same
research group flagged exactly this kind of slip: forgetting a term
when moving a squared-norm quantity from one side of an inequality to
the other before invoking a discrete Gronwall argument. Both lemmas
below are checked by Lean's kernel, so no such slip is possible here.
-/

/-- **Polarization identity.**
This is the elementary identity that converts a mixed product term
`2 a (a - b)` (as it appears when one takes the inner product of a
semi-implicit update with the increment `u^{n+1} - u^n`) into a
telescoping difference of squares plus a nonnegative remainder. It is
the algebraic core of every discrete energy-dissipation proof for
gradient-flow schemes (e.g. SAV-type methods), including the
`r^{n+1}(r^{n+1}-r^n)` step used to close the energy estimate for the
auxiliary variable. -/
theorem polarization_identity (a b : ℝ) :
    2 * a * (a - b) = a ^ 2 - b ^ 2 + (a - b) ^ 2 := by
  ring

/-- **Discrete Grönwall step.**
Suppose a one-step error estimate of the form

  `a^2 ≤ b^2 + τ(L+1)(a^2 + b^2)`

has been derived for a semi-implicit scheme, where `a = ‖e^{n+1}‖`,
`b = ‖e^{n}‖`, `L` is the Lipschitz constant of the nonlinearity, and
`τ` is the time step. (This is exactly the point at which the review
of an earlier draft in this project flagged a missing step: the term
`τ(L+1) a^2` on the right-hand side must be moved to the left and the
inequality divided by `1 - τ(L+1)` -- *not* dropped.) Under the
stability restriction `τ(L+1) < 1`, this lemma proves the sharpened,
explicit one-step bound
  `a^2 ≤ ((1+τ(L+1)) / (1-τ(L+1))) * b^2`,
from which the standard discrete Grönwall lemma (iterating over `n`)
yields exponential-in-time stability. -/
theorem discrete_gronwall_step (a b L τ : ℝ) (hL : 0 ≤ L) (hτ : 0 < τ)
    (hτL : τ * (L + 1) < 1)
    (h : a ^ 2 ≤ b ^ 2 + τ * (L + 1) * (a ^ 2 + b ^ 2)) :
    a ^ 2 ≤ (1 + τ * (L + 1)) / (1 - τ * (L + 1)) * b ^ 2 := by
  set c := τ * (L + 1) with hc_def
  have hc_nonneg : 0 ≤ c := by
    have hL1 : (0:ℝ) ≤ L + 1 := by linarith
    exact mul_nonneg hτ.le hL1
  have hc_lt_one : c < 1 := hτL
  have hpos : (0:ℝ) < 1 - c := by linarith
  -- Rearrange `h` into `(1 - c) * a^2 ≤ (1 + c) * b^2`.
  have key : (1 - c) * a ^ 2 ≤ (1 + c) * b ^ 2 := by nlinarith [h]
  rw [div_mul_eq_mul_div, le_div_iff₀ hpos]
  nlinarith [key]

/-- **Sanity check.** The stability restriction `τ(L+1) < 1` is not
vacuous: for a fixed Lipschitz constant `L`, it simply bounds the time
step `τ < 1/(L+1)`, which is always achievable by taking `τ` small
enough. This closes the loop between the abstract lemma above and its
use as a genuine, satisfiable hypothesis in the numerical scheme. -/
theorem stability_step_size_exists (L : ℝ) (hL : 0 ≤ L) :
    ∃ τ : ℝ, 0 < τ ∧ τ * (L + 1) < 1 := by
  have hL1 : (0:ℝ) < L + 1 := by linarith
  have hden : (0:ℝ) < 2 * (L + 1) := by linarith
  refine ⟨1 / (2 * (L + 1)), div_pos one_pos hden, ?_⟩
  rw [div_mul_eq_mul_div, div_lt_one hden]
  linarith
