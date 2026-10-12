import Mathlib

set_option linter.style.header false

/-!
# A-posteriori energy balance for an arbitrary trajectory (perturbed gradient flow)

`GradientFlow.lean` proves that an *exact* gradient flow `c'(t) = -grad E(c(t))`
in a real Hilbert space dissipates the energy at the rate `-‖c'‖^2`.  A trained
Spectral-PINN does **not** satisfy the PDE exactly: its coefficient trajectory
`c(t)` satisfies
```
  c'(t) = -grad E(c(t)) + r(t)
```
where `r(t)` is the PDE residual (the quantity that the training loss
minimises).  This file proves, for an *arbitrary* differentiable trajectory,
the resulting a-posteriori energy balance, so that the conclusion no longer
depends on the trajectory being an exact gradient flow:

1.  `perturbedFlow_hasDerivAt_energy`:  `d/dt E(c(t)) = -‖c'‖^2 + <r, c'>` (exact identity);
2.  `perturbedFlow_deriv_energy_le`:  `d/dt E(c(t)) <= (1/2)‖r‖^2 - (1/2)‖c'‖^2`;
3.  `defectFlow_*`: the version in which the differential of `E` is only
    `DE(c)[v] = <G, v> + beta(v)` with a defect functional `beta`
    (e.g. `eps^2 [u_x v]_{-1}^{1}` when the Neumann condition is not exact), and
    the PDE holds as `c' = -G + r` in a Hilbert space `H` into which the space
    `V` carrying `E` is embedded (`V = ` coefficient space / `H^1`, `H = L^2`);
4.  `*_energy_growth_bound`: if `‖r(t)‖ <= delta` on `[0,T]` then
    `E(c(t)) <= E(c(0)) + (delta^2/2) t`, and with `|beta(c'(t))| <= b`,
    `E(c(t)) <= E(c(0)) + (delta^2/2 + b) t`;
5.  `discrete_energy_balance`: the discrete-time analogue;
6.  `first_variation_split`: the spatial integration by parts for polynomials that
    identifies `G = -eps^2 u_xx + f(u)` and `β(w) = eps^2 [u_x w]_{-1}^{1}` for Allen-Cahn.

Remark on constants: the inequality (2) keeps the dissipation term and has constant `1/2`.
If the dissipation term is discarded, the sharp constant is `1/4`
(`inner_balance_quarter`, `*_quarter` theorems); it is attained when `c' = r/2 = grad E`.

Everything is proved for abstract Hilbert spaces by Mathlib's chain rule and
mean-value theorem, with no unproven placeholders and only the three standard axioms
(see the `#print axioms` lines at the end).
-/

open scoped Gradient RealInnerProductSpace
open InnerProductSpace

/-! ## Section 1. The algebra of one instant (pure inner-product-space facts) -/

section Algebra

variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- If `w = -G + r` then `<G, w> = -‖w‖^2 + <r, w>`.  (Here `w = c'`, `G = grad E`,
`r` = residual.)  For `r = 0` this is the exact-gradient-flow identity
`<G, w> = -‖w‖^2` used in `GradientFlow.lean`. -/
theorem inner_balance (G r w : H) (hw : w = -G + r) :
    ⟪G, w⟫_ℝ = -‖w‖ ^ 2 + ⟪r, w⟫_ℝ := by
  have hG : G = r - w := by rw [hw]; abel
  rw [hG, inner_sub_left, real_inner_self_eq_norm_sq]
  ring

/-- Sharper, exact form of the balance: `<G, w> = (1/2)‖r‖^2 - (1/2)‖w‖^2 - (1/2)‖G‖^2`
whenever `w = -G + r`.  (Follows from `‖r‖^2 = ‖w + G‖^2 = ‖w‖^2 + 2<w,G> + ‖G‖^2`.) -/
theorem inner_balance_exact (G r w : H) (hw : w = -G + r) :
    ⟪G, w⟫_ℝ = (1 / 2) * ‖r‖ ^ 2 - (1 / 2) * ‖w‖ ^ 2 - (1 / 2) * ‖G‖ ^ 2 := by
  have hr : r = w + G := by rw [hw]; abel
  have h2 : ‖r‖ ^ 2 = ‖w‖ ^ 2 + 2 * ⟪w, G⟫_ℝ + ‖G‖ ^ 2 := by
    rw [hr]; exact norm_add_sq_real w G
  rw [real_inner_comm] at h2
  linarith

/-- Young-type bound used for the a-posteriori estimate:
`<G, w> <= (1/2)‖r‖^2 - (1/2)‖w‖^2` whenever `w = -G + r`. -/
theorem inner_balance_le (G r w : H) (hw : w = -G + r) :
    ⟪G, w⟫_ℝ ≤ (1 / 2) * ‖r‖ ^ 2 - (1 / 2) * ‖w‖ ^ 2 := by
  have := inner_balance_exact G r w hw
  nlinarith [sq_nonneg ‖G‖]

/-- Completing the square: `<G, w> = (1/4)‖r‖^2 - ‖G - (1/2) r‖^2` whenever `w = -G + r`.
Hence the best constant that bounds `<G, w>` by a multiple of `‖r‖^2` alone is `1/4`
(attained for `G = r/2`, i.e. `w = r/2`), whereas `inner_balance_le` keeps the
dissipation `-(1/2)‖w‖^2` at the price of the constant `1/2`. -/
theorem inner_balance_quarter (G r w : H) (hw : w = -G + r) :
    ⟪G, w⟫_ℝ = (1 / 4) * ‖r‖ ^ 2 - ‖G - (1 / 2 : ℝ) • r‖ ^ 2 := by
  have hw' : w = r - G := by rw [hw]; abel
  rw [hw', inner_sub_right, real_inner_self_eq_norm_sq, norm_sub_sq_real, inner_smul_right,
    norm_smul, Real.norm_eq_abs, abs_of_pos (by norm_num : (0 : ℝ) < 1 / 2)]
  rw [real_inner_comm]
  ring

/-- Sharp bound: `<G, w> <= (1/4)‖r‖^2` whenever `w = -G + r`. -/
theorem inner_balance_le_quarter (G r w : H) (hw : w = -G + r) :
    ⟪G, w⟫_ℝ ≤ (1 / 4) * ‖r‖ ^ 2 := by
  have := inner_balance_quarter G r w hw
  nlinarith [sq_nonneg ‖G - (1 / 2 : ℝ) • r‖]

end Algebra

/-! ## Section 2. Perturbed gradient flow in a Hilbert space `F` (items (1) and (2)) -/

section Perturbed

variable {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F] [CompleteSpace F]

/-- Chain rule: if `E` has gradient `G` at `c t` and `c` has derivative `c'` at `t`,
then `s ↦ E (c s)` has derivative `<G, c'>` at `t` (no flow hypothesis). -/
theorem hasDerivAt_energy_inner
    (E : F → ℝ) (c : ℝ → F) (c' G : F) (t : ℝ)
    (hE : HasGradientAt E G (c t)) (hc : HasDerivAt c c' t) :
    HasDerivAt (fun s => E (c s)) ⟪G, c'⟫_ℝ t := by
  have hFD : HasFDerivAt E (toDual ℝ F G) (c t) := hE.hasFDerivAt
  have hcomp : HasDerivAt (E ∘ c) (toDual ℝ F G c') t :=
    HasFDerivAt.comp_hasDerivAt_of_eq t hFD hc rfl
  rw [toDual_apply_apply] at hcomp
  exact hcomp

/-- **(1) Pointwise energy balance for the perturbed gradient flow** (exact identity).
Let `E` have gradient `G` at `c t`, let `c` be differentiable at `t` with derivative `c'`,
and suppose `c' = -G + r` (the flow holds up to the residual `r`).  Then
`d/dt E(c(t)) = -‖c'‖^2 + <r, c'>`. -/
theorem perturbedFlow_hasDerivAt_energy
    (E : F → ℝ) (c : ℝ → F) (c' G r : F) (t : ℝ)
    (hE : HasGradientAt E G (c t)) (hc : HasDerivAt c c' t)
    (hflow : c' = -G + r) :
    HasDerivAt (fun s => E (c s)) (-‖c'‖ ^ 2 + ⟪r, c'⟫_ℝ) t := by
  have h := hasDerivAt_energy_inner E c c' G t hE hc
  rwa [inner_balance G r c' hflow] at h

/-- Same identity in the format of `GradientFlow.lean` (`gradient E (c t)`; the flow
hypothesis is `HasDerivAt c (-(gradient E (c t)) + r) t`). -/
theorem perturbedFlow_hasDerivAt_energy_gradient
    (E : F → ℝ) (c : ℝ → F) (r : F) (t : ℝ)
    (hE : HasGradientAt E (gradient E (c t)) (c t))
    (hflow : HasDerivAt c (-(gradient E (c t)) + r) t) :
    HasDerivAt (fun s => E (c s)) (-‖deriv c t‖ ^ 2 + ⟪r, deriv c t⟫_ℝ) t := by
  have hd : deriv c t = -(gradient E (c t)) + r := hflow.deriv
  rw [hd]
  exact perturbedFlow_hasDerivAt_energy E c _ (gradient E (c t)) r t hE hflow rfl

/-- **Exact refined balance**: `d/dt E(c(t)) = (1/2)‖r‖^2 - (1/2)‖c'‖^2 - (1/2)‖G‖^2`. -/
theorem perturbedFlow_hasDerivAt_energy_exact
    (E : F → ℝ) (c : ℝ → F) (c' G r : F) (t : ℝ)
    (hE : HasGradientAt E G (c t)) (hc : HasDerivAt c c' t)
    (hflow : c' = -G + r) :
    HasDerivAt (fun s => E (c s))
      ((1 / 2) * ‖r‖ ^ 2 - (1 / 2) * ‖c'‖ ^ 2 - (1 / 2) * ‖G‖ ^ 2) t := by
  have h := hasDerivAt_energy_inner E c c' G t hE hc
  rwa [inner_balance_exact G r c' hflow] at h

/-- **(2) A-posteriori energy inequality**:
`d/dt E(c(t)) <= (1/2)‖r‖^2 - (1/2)‖c'‖^2`. -/
theorem perturbedFlow_deriv_energy_le
    (E : F → ℝ) (c : ℝ → F) (c' G r : F) (t : ℝ)
    (hE : HasGradientAt E G (c t)) (hc : HasDerivAt c c' t)
    (hflow : c' = -G + r) :
    deriv (fun s => E (c s)) t ≤ (1 / 2) * ‖r‖ ^ 2 - (1 / 2) * ‖c'‖ ^ 2 := by
  rw [(hasDerivAt_energy_inner E c c' G t hE hc).deriv]
  exact inner_balance_le G r c' hflow

/-- The same inequality with the derivative written as `deriv c t`. -/
theorem perturbedFlow_deriv_energy_le_deriv
    (E : F → ℝ) (c : ℝ → F) (G r : F) (t : ℝ)
    (hE : HasGradientAt E G (c t)) (hc : DifferentiableAt ℝ c t)
    (hflow : deriv c t = -G + r) :
    deriv (fun s => E (c s)) t ≤ (1 / 2) * ‖r‖ ^ 2 - (1 / 2) * ‖deriv c t‖ ^ 2 :=
  perturbedFlow_deriv_energy_le E c (deriv c t) G r t hE hc.hasDerivAt hflow

/-- Consistency with `GradientFlow.lean`: for `r = 0` the identity (1) is exactly
`d/dt E(c(t)) = -‖c'‖^2`. -/
theorem perturbedFlow_exact_case
    (E : F → ℝ) (c : ℝ → F) (c' G : F) (t : ℝ)
    (hE : HasGradientAt E G (c t)) (hc : HasDerivAt c c' t)
    (hflow : c' = -G) :
    HasDerivAt (fun s => E (c s)) (-‖c'‖ ^ 2) t := by
  have h := perturbedFlow_hasDerivAt_energy E c c' G 0 t hE hc (by simpa using hflow)
  simpa using h

/-- If the residual is no larger than the speed, `‖r‖ <= ‖c'‖`, the energy is still
instantaneously non-increasing. -/
theorem perturbedFlow_deriv_energy_nonpos_of_residual_le
    (E : F → ℝ) (c : ℝ → F) (c' G r : F) (t : ℝ)
    (hE : HasGradientAt E G (c t)) (hc : HasDerivAt c c' t)
    (hflow : c' = -G + r) (hrc : ‖r‖ ≤ ‖c'‖) :
    deriv (fun s => E (c s)) t ≤ 0 := by
  have h := perturbedFlow_deriv_energy_le E c c' G r t hE hc hflow
  have h2 : ‖r‖ ^ 2 ≤ ‖c'‖ ^ 2 := pow_le_pow_left₀ (norm_nonneg _) hrc 2
  linarith

/-- Relative-residual version: if `‖r‖ <= θ ‖c'‖` with `0 <= θ`, then
`d/dt E(c(t)) <= -(1 - θ) ‖c'‖^2` (strict dissipation for `θ < 1`). -/
theorem perturbedFlow_deriv_energy_le_of_relative
    (E : F → ℝ) (c : ℝ → F) (c' G r : F) (θ : ℝ) (t : ℝ)
    (hE : HasGradientAt E G (c t)) (hc : HasDerivAt c c' t)
    (hflow : c' = -G + r) (hrel : ‖r‖ ≤ θ * ‖c'‖) :
    deriv (fun s => E (c s)) t ≤ -(1 - θ) * ‖c'‖ ^ 2 := by
  rw [(perturbedFlow_hasDerivAt_energy E c c' G r t hE hc hflow).deriv]
  have h1 : ⟪r, c'⟫_ℝ ≤ ‖r‖ * ‖c'‖ := real_inner_le_norm r c'
  have h2 : ‖r‖ * ‖c'‖ ≤ θ * ‖c'‖ * ‖c'‖ :=
    mul_le_mul_of_nonneg_right hrel (norm_nonneg _)
  nlinarith

end Perturbed

/-! ## Section 3. Defect functional (item (3))

The energy `E` is differentiable on a normed space `V` (coefficient space, `H^1`, ...)
which is embedded by a continuous linear map `ι : V → H` into the Hilbert space `H`
(`L^2`) in which the PDE residual is measured.  The differential `DE(c) : V → ℝ` is *not*
assumed to be an `H`-inner product with a gradient: it only splits as
`DE(c)[v] = <G, ι v> + β v`, where `G ∈ H` is the "strong-form gradient"
(for Allen-Cahn: `G = -eps^2 u_xx + f(u)`) and `β` is a defect functional
(for Allen-Cahn: `β(v) = eps^2 [u_x v]_{-1}^{1}`, which vanishes iff the Neumann
condition holds).  The PDE is satisfied up to a residual `r ∈ H`:
`ι c' = -G + r`.  The case `V = H`, `ι = id`, `β = 0`, `G = grad E` is Section 2. -/

section Defect

variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]
variable {H : Type*} [NormedAddCommGroup H] [InnerProductSpace ℝ H]

/-- **(3) Energy balance with a defect.**  If `DE(c')= <G, ι c'> + β` (with `β` a real
number, the value of the defect functional at the velocity) and `ι c' = -G + r`, then
`d/dt E(c(t)) = -‖ι c'‖^2 + <r, ι c'> + β`. -/
theorem defectFlow_hasDerivAt_energy
    (E : V → ℝ) (DE : V →L[ℝ] ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (c' : V) (G r : H)
    (β : ℝ) (t : ℝ)
    (hE : HasFDerivAt E DE (c t)) (hc : HasDerivAt c c' t)
    (hDE : DE c' = ⟪G, ι c'⟫_ℝ + β)
    (hflow : ι c' = -G + r) :
    HasDerivAt (fun s => E (c s)) (-‖ι c'‖ ^ 2 + ⟪r, ι c'⟫_ℝ + β) t := by
  have hcomp : HasDerivAt (E ∘ c) (DE c') t :=
    HasFDerivAt.comp_hasDerivAt_of_eq t hE hc rfl
  rw [hDE, inner_balance G r (ι c') hflow] at hcomp
  exact hcomp

/-- Version of (3) with a defect *functional* `β : V →L[ℝ] ℝ`: if
`DE v = <G, ι v> + β v` for all `v`, then `d/dt E(c(t)) = -‖ι c'‖^2 + <r, ι c'> + β c'`. -/
theorem defectFlow_hasDerivAt_energy_functional
    (E : V → ℝ) (DE β : V →L[ℝ] ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (c' : V) (G r : H)
    (t : ℝ)
    (hE : HasFDerivAt E DE (c t)) (hc : HasDerivAt c c' t)
    (hDE : ∀ v, DE v = ⟪G, ι v⟫_ℝ + β v)
    (hflow : ι c' = -G + r) :
    HasDerivAt (fun s => E (c s)) (-‖ι c'‖ ^ 2 + ⟪r, ι c'⟫_ℝ + β c') t :=
  defectFlow_hasDerivAt_energy E DE ι c c' G r (β c') t hE hc (hDE c') hflow

/-- Exact refined balance with defect:
`d/dt E = (1/2)‖r‖^2 - (1/2)‖ι c'‖^2 - (1/2)‖G‖^2 + β`. -/
theorem defectFlow_hasDerivAt_energy_exact
    (E : V → ℝ) (DE : V →L[ℝ] ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (c' : V) (G r : H)
    (β : ℝ) (t : ℝ)
    (hE : HasFDerivAt E DE (c t)) (hc : HasDerivAt c c' t)
    (hDE : DE c' = ⟪G, ι c'⟫_ℝ + β)
    (hflow : ι c' = -G + r) :
    HasDerivAt (fun s => E (c s))
      ((1 / 2) * ‖r‖ ^ 2 - (1 / 2) * ‖ι c'‖ ^ 2 - (1 / 2) * ‖G‖ ^ 2 + β) t := by
  have hcomp : HasDerivAt (E ∘ c) (DE c') t :=
    HasFDerivAt.comp_hasDerivAt_of_eq t hE hc rfl
  rw [hDE, inner_balance_exact G r (ι c') hflow] at hcomp
  exact hcomp

/-- **(3') Inequality with defect**: `d/dt E(c(t)) <= (1/2)‖r‖^2 - (1/2)‖ι c'‖^2 + β`. -/
theorem defectFlow_deriv_energy_le
    (E : V → ℝ) (DE : V →L[ℝ] ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (c' : V) (G r : H)
    (β : ℝ) (t : ℝ)
    (hE : HasFDerivAt E DE (c t)) (hc : HasDerivAt c c' t)
    (hDE : DE c' = ⟪G, ι c'⟫_ℝ + β)
    (hflow : ι c' = -G + r) :
    deriv (fun s => E (c s)) t ≤ (1 / 2) * ‖r‖ ^ 2 - (1 / 2) * ‖ι c'‖ ^ 2 + β := by
  have hcomp : HasDerivAt (E ∘ c) (DE c') t :=
    HasFDerivAt.comp_hasDerivAt_of_eq t hE hc rfl
  have h : deriv (fun s => E (c s)) t = DE c' := hcomp.deriv
  rw [h, hDE]
  have := inner_balance_le G r (ι c') hflow
  linarith

/-- Sharp pointwise bound with defect, discarding the dissipation term:
`d/dt E(c(t)) <= (1/4)‖r‖^2 + β`. -/
theorem defectFlow_deriv_energy_le_quarter
    (E : V → ℝ) (DE : V →L[ℝ] ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (c' : V) (G r : H)
    (β : ℝ) (t : ℝ)
    (hE : HasFDerivAt E DE (c t)) (hc : HasDerivAt c c' t)
    (hDE : DE c' = ⟪G, ι c'⟫_ℝ + β)
    (hflow : ι c' = -G + r) :
    deriv (fun s => E (c s)) t ≤ (1 / 4) * ‖r‖ ^ 2 + β := by
  have hcomp : HasDerivAt (E ∘ c) (DE c') t :=
    HasFDerivAt.comp_hasDerivAt_of_eq t hE hc rfl
  have h : deriv (fun s => E (c s)) t = DE c' := hcomp.deriv
  rw [h, hDE]
  have := inner_balance_le_quarter G r (ι c') hflow
  linarith

/-- Sanity check: with `β = 0`, `V = H`, `ι = id` the defect balance (3) contains the
balance (1) (here stated for a general inner-product space `H` and the derivative
`DE` of `E`). -/
theorem defectFlow_zero_defect
    (E : H → ℝ) (DE : H →L[ℝ] ℝ) (c : ℝ → H) (c' G r : H) (t : ℝ)
    (hE : HasFDerivAt E DE (c t)) (hc : HasDerivAt c c' t)
    (hDE : DE c' = ⟪G, c'⟫_ℝ)
    (hflow : c' = -G + r) :
    HasDerivAt (fun s => E (c s)) (-‖c'‖ ^ 2 + ⟪r, c'⟫_ℝ) t := by
  have h := defectFlow_hasDerivAt_energy E DE (ContinuousLinearMap.id ℝ H) c c' G r 0 t hE hc
    (by simpa using hDE) (by simpa using hflow)
  simpa using h

/-! ## Section 4. Integrated consequences (item (4)) -/

/-- Mean-value step: if `ψ` is differentiable on `[0,T]` and `ψ' <= C` on `(0,T)`, then
`ψ t <= ψ s + C (t - s)` for `0 <= s <= t <= T`. -/
theorem energy_le_of_deriv_le (ψ : ℝ → ℝ) (T C : ℝ)
    (hψ : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ ψ t)
    (hd : ∀ t ∈ Set.Ioo 0 T, deriv ψ t ≤ C) :
    ∀ s ∈ Set.Icc 0 T, ∀ t ∈ Set.Icc 0 T, s ≤ t → ψ t ≤ ψ s + C * (t - s) := by
  intro s hs t ht hst
  have hcont : ContinuousOn ψ (Set.Icc 0 T) :=
    fun u hu => (hψ u hu).continuousAt.continuousWithinAt
  have hdiff : DifferentiableOn ℝ ψ (interior (Set.Icc 0 T)) := by
    rw [interior_Icc]
    exact fun u hu => (hψ u (Set.Ioo_subset_Icc_self hu)).differentiableWithinAt
  have hd' : ∀ u ∈ interior (Set.Icc 0 T), deriv ψ u ≤ C := by
    rw [interior_Icc]; exact hd
  have := (convex_Icc (0 : ℝ) T).image_sub_le_mul_sub_of_deriv_le hcont hdiff hd' s hs t ht hst
  linarith

/-- **(4, with defect) Integrated energy bound.**  Along a trajectory `c` on `[0,T]` with
`ι c' = -G + r`, `DE(c)[c'] = <G, ι c'> + β(t)`, residual bound `‖r(t)‖ <= δ` and
defect bound `|β(t)| <= b`, we have for `0 <= s <= t <= T`
`E(c(t)) <= E(c(s)) + (δ^2/2 + b) (t - s)`. -/
theorem defectFlow_energy_growth_bound
    (E : V → ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (G r : ℝ → H) (β : ℝ → ℝ) (T δ b : ℝ)
    (hEd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ E (c t))
    (hcd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ c t)
    (hDE : ∀ t ∈ Set.Icc 0 T, fderiv ℝ E (c t) (deriv c t) = ⟪G t, ι (deriv c t)⟫_ℝ + β t)
    (hflow : ∀ t ∈ Set.Icc 0 T, ι (deriv c t) = -G t + r t)
    (hr : ∀ t ∈ Set.Icc 0 T, ‖r t‖ ≤ δ)
    (hβ : ∀ t ∈ Set.Icc 0 T, |β t| ≤ b) :
    ∀ s ∈ Set.Icc 0 T, ∀ t ∈ Set.Icc 0 T, s ≤ t →
      E (c t) ≤ E (c s) + (δ ^ 2 / 2 + b) * (t - s) := by
  have hcomp : ∀ t ∈ Set.Icc 0 T,
      HasDerivAt (fun u => E (c u)) (fderiv ℝ E (c t) (deriv c t)) t := by
    intro t ht
    exact HasFDerivAt.comp_hasDerivAt_of_eq t (hEd t ht).hasFDerivAt (hcd t ht).hasDerivAt rfl
  refine energy_le_of_deriv_le (fun u => E (c u)) T (δ ^ 2 / 2 + b)
    (fun t ht => (hcomp t ht).differentiableAt) ?_
  intro t ht
  have ht' : t ∈ Set.Icc 0 T := Set.Ioo_subset_Icc_self ht
  have h1 := defectFlow_deriv_energy_le E (fderiv ℝ E (c t)) ι c (deriv c t) (G t) (r t)
    (β t) t (hEd t ht').hasFDerivAt (hcd t ht').hasDerivAt (hDE t ht') (hflow t ht')
  have h2 : ‖r t‖ ^ 2 ≤ δ ^ 2 := pow_le_pow_left₀ (norm_nonneg _) (hr t ht') 2
  have h3 : β t ≤ b := (le_abs_self _).trans (hβ t ht')
  have h4 : 0 ≤ ‖ι (deriv c t)‖ ^ 2 := sq_nonneg _
  linarith

/-- **(4, with defect), from time zero**:
`E(c(t)) <= E(c(0)) + (δ^2/2 + b) t` for `t ∈ [0,T]`. -/
theorem defectFlow_energy_bound
    (E : V → ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (G r : ℝ → H) (β : ℝ → ℝ) (T δ b : ℝ)
    (hEd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ E (c t))
    (hcd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ c t)
    (hDE : ∀ t ∈ Set.Icc 0 T, fderiv ℝ E (c t) (deriv c t) = ⟪G t, ι (deriv c t)⟫_ℝ + β t)
    (hflow : ∀ t ∈ Set.Icc 0 T, ι (deriv c t) = -G t + r t)
    (hr : ∀ t ∈ Set.Icc 0 T, ‖r t‖ ≤ δ)
    (hβ : ∀ t ∈ Set.Icc 0 T, |β t| ≤ b) :
    ∀ t ∈ Set.Icc 0 T, E (c t) ≤ E (c 0) + (δ ^ 2 / 2 + b) * t := by
  intro t ht
  have h0 : (0 : ℝ) ∈ Set.Icc 0 T := ⟨le_rfl, ht.1.trans ht.2⟩
  have := defectFlow_energy_growth_bound E ι c G r β T δ b hEd hcd hDE hflow hr hβ 0 h0 t ht ht.1
  simpa using this

/-- **(4'), energy estimate in integral form.**  If `ψ` has a continuous derivative `ψ'` on
`[0,T]` and `ψ' <= (1/2) ρ - (1/2) d + b` pointwise, with `ρ, d, b` continuous, then
`ψ(t) + (1/2) ∫_0^t d <= ψ(0) + (1/2) ∫_0^t ρ + ∫_0^t b`. -/
theorem energy_integral_le (ψ ψ' ρ d b : ℝ → ℝ) (T : ℝ)
    (hψ : ∀ t ∈ Set.Icc 0 T, HasDerivAt ψ (ψ' t) t)
    (hψ' : ContinuousOn ψ' (Set.Icc 0 T))
    (hρ : ContinuousOn ρ (Set.Icc 0 T)) (hd : ContinuousOn d (Set.Icc 0 T))
    (hb : ContinuousOn b (Set.Icc 0 T))
    (hle : ∀ t ∈ Set.Icc 0 T, ψ' t ≤ (1 / 2) * ρ t - (1 / 2) * d t + b t) :
    ∀ t ∈ Set.Icc 0 T, ψ t + (1 / 2) * ∫ s in (0 : ℝ)..t, d s
      ≤ ψ 0 + (1 / 2) * (∫ s in (0 : ℝ)..t, ρ s) + ∫ s in (0 : ℝ)..t, b s := by
  intro t ht
  have ht0 : (0 : ℝ) ≤ t := ht.1
  have hsub : Set.uIcc (0 : ℝ) t ⊆ Set.Icc 0 T := by
    rw [Set.uIcc_of_le ht0]; exact Set.Icc_subset_Icc_right ht.2
  have hsub' : Set.Icc (0 : ℝ) t ⊆ Set.Icc 0 T := Set.Icc_subset_Icc_right ht.2
  have iψ' : IntervalIntegrable ψ' MeasureTheory.volume 0 t :=
    (hψ'.mono hsub).intervalIntegrable
  have iρ : IntervalIntegrable ρ MeasureTheory.volume 0 t := (hρ.mono hsub).intervalIntegrable
  have id' : IntervalIntegrable d MeasureTheory.volume 0 t := (hd.mono hsub).intervalIntegrable
  have ib : IntervalIntegrable b MeasureTheory.volume 0 t := (hb.mono hsub).intervalIntegrable
  have hFTC : ∫ s in (0 : ℝ)..t, ψ' s = ψ t - ψ 0 :=
    intervalIntegral.integral_eq_sub_of_hasDerivAt (fun s hs => hψ s (hsub hs)) iψ'
  have iR : IntervalIntegrable (fun s => (1 / 2) * ρ s - (1 / 2) * d s + b s)
      MeasureTheory.volume 0 t :=
    ((iρ.const_mul _).sub (id'.const_mul _)).add ib
  have hmono : ∫ s in (0 : ℝ)..t, ψ' s
      ≤ ∫ s in (0 : ℝ)..t, ((1 / 2) * ρ s - (1 / 2) * d s + b s) :=
    intervalIntegral.integral_mono_on ht0 iψ' iR (fun s hs => hle s (hsub' hs))
  have hlin : ∫ s in (0 : ℝ)..t, ((1 / 2) * ρ s - (1 / 2) * d s + b s)
      = (1 / 2) * (∫ s in (0 : ℝ)..t, ρ s) - (1 / 2) * (∫ s in (0 : ℝ)..t, d s)
        + ∫ s in (0 : ℝ)..t, b s := by
    rw [intervalIntegral.integral_add ((iρ.const_mul _).sub (id'.const_mul _)) ib,
      intervalIntegral.integral_sub (iρ.const_mul _) (id'.const_mul _),
      intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul]
  linarith

/-- **(4', with defect), integral form along a trajectory.**
`E(c(t)) + (1/2) ∫_0^t ‖ι c'‖^2 <= E(c(0)) + (1/2) ∫_0^t ‖r‖^2 + ∫_0^t β`. -/
theorem defectFlow_energy_integral_bound
    (E : V → ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (G r : ℝ → H) (β : ℝ → ℝ) (T : ℝ)
    (hEd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ E (c t))
    (hcd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ c t)
    (hDE : ∀ t ∈ Set.Icc 0 T, fderiv ℝ E (c t) (deriv c t) = ⟪G t, ι (deriv c t)⟫_ℝ + β t)
    (hflow : ∀ t ∈ Set.Icc 0 T, ι (deriv c t) = -G t + r t)
    (hcont : ContinuousOn (fun t => fderiv ℝ E (c t) (deriv c t)) (Set.Icc 0 T))
    (hρ : ContinuousOn (fun t => ‖r t‖ ^ 2) (Set.Icc 0 T))
    (hd : ContinuousOn (fun t => ‖ι (deriv c t)‖ ^ 2) (Set.Icc 0 T))
    (hb : ContinuousOn β (Set.Icc 0 T)) :
    ∀ t ∈ Set.Icc 0 T,
      E (c t) + (1 / 2) * ∫ s in (0 : ℝ)..t, ‖ι (deriv c s)‖ ^ 2
        ≤ E (c 0) + (1 / 2) * (∫ s in (0 : ℝ)..t, ‖r s‖ ^ 2) + ∫ s in (0 : ℝ)..t, β s := by
  refine energy_integral_le (fun u => E (c u)) (fun t => fderiv ℝ E (c t) (deriv c t))
    (fun t => ‖r t‖ ^ 2) (fun t => ‖ι (deriv c t)‖ ^ 2) β T ?_ hcont hρ hd hb ?_
  · intro t ht
    exact HasFDerivAt.comp_hasDerivAt_of_eq t (hEd t ht).hasFDerivAt (hcd t ht).hasDerivAt rfl
  · intro t ht
    rw [hDE t ht]
    have := inner_balance_le (G t) (r t) (ι (deriv c t)) (hflow t ht)
    linarith

/-- **(4', zero defect)**: if `DE(c)[c'] = <G, ι c'>` (no boundary defect) and
`ι c' = -G + r`, then
`E(c(t)) + (1/2) ∫_0^t ‖ι c'‖^2 <= E(c(0)) + (1/2) ∫_0^t ‖r‖^2`. -/
theorem defectFree_energy_integral_bound
    (E : V → ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (G r : ℝ → H) (T : ℝ)
    (hEd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ E (c t))
    (hcd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ c t)
    (hDE : ∀ t ∈ Set.Icc 0 T, fderiv ℝ E (c t) (deriv c t) = ⟪G t, ι (deriv c t)⟫_ℝ)
    (hflow : ∀ t ∈ Set.Icc 0 T, ι (deriv c t) = -G t + r t)
    (hcont : ContinuousOn (fun t => fderiv ℝ E (c t) (deriv c t)) (Set.Icc 0 T))
    (hρ : ContinuousOn (fun t => ‖r t‖ ^ 2) (Set.Icc 0 T))
    (hd : ContinuousOn (fun t => ‖ι (deriv c t)‖ ^ 2) (Set.Icc 0 T)) :
    ∀ t ∈ Set.Icc 0 T,
      E (c t) + (1 / 2) * ∫ s in (0 : ℝ)..t, ‖ι (deriv c s)‖ ^ 2
        ≤ E (c 0) + (1 / 2) * ∫ s in (0 : ℝ)..t, ‖r s‖ ^ 2 := by
  intro t ht
  have := defectFlow_energy_integral_bound E ι c G r (fun _ => 0) T hEd hcd
    (by intro u hu; simpa using hDE u hu) hflow hcont hρ hd continuousOn_const t ht
  simpa using this

/-- **(4), sharp constant, with defect**, from time zero:
`E(c(t)) <= E(c(0)) + (δ^2/4 + b) t`. -/
theorem defectFlow_energy_bound_quarter
    (E : V → ℝ) (ι : V →L[ℝ] H) (c : ℝ → V) (G r : ℝ → H) (β : ℝ → ℝ) (T δ b : ℝ)
    (hEd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ E (c t))
    (hcd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ c t)
    (hDE : ∀ t ∈ Set.Icc 0 T, fderiv ℝ E (c t) (deriv c t) = ⟪G t, ι (deriv c t)⟫_ℝ + β t)
    (hflow : ∀ t ∈ Set.Icc 0 T, ι (deriv c t) = -G t + r t)
    (hr : ∀ t ∈ Set.Icc 0 T, ‖r t‖ ≤ δ)
    (hβ : ∀ t ∈ Set.Icc 0 T, |β t| ≤ b) :
    ∀ t ∈ Set.Icc 0 T, E (c t) ≤ E (c 0) + (δ ^ 2 / 4 + b) * t := by
  have hcomp : ∀ t ∈ Set.Icc 0 T,
      HasDerivAt (fun u => E (c u)) (fderiv ℝ E (c t) (deriv c t)) t := by
    intro t ht
    exact HasFDerivAt.comp_hasDerivAt_of_eq t (hEd t ht).hasFDerivAt (hcd t ht).hasDerivAt rfl
  have key := energy_le_of_deriv_le (fun u => E (c u)) T (δ ^ 2 / 4 + b)
    (fun t ht => (hcomp t ht).differentiableAt) (by
      intro t ht
      have ht' : t ∈ Set.Icc 0 T := Set.Ioo_subset_Icc_self ht
      have h1 := defectFlow_deriv_energy_le_quarter E (fderiv ℝ E (c t)) ι c (deriv c t) (G t)
        (r t) (β t) t (hEd t ht').hasFDerivAt (hcd t ht').hasDerivAt (hDE t ht') (hflow t ht')
      have h2 : ‖r t‖ ^ 2 ≤ δ ^ 2 := pow_le_pow_left₀ (norm_nonneg _) (hr t ht') 2
      have h3 : β t ≤ b := (le_abs_self _).trans (hβ t ht')
      linarith)
  intro t ht
  have h0 : (0 : ℝ) ∈ Set.Icc 0 T := ⟨le_rfl, ht.1.trans ht.2⟩
  have := key 0 h0 t ht ht.1
  simpa using this

end Defect

section PerturbedIntegrated

variable {F : Type*} [NormedAddCommGroup F] [InnerProductSpace ℝ F] [CompleteSpace F]

/-- **(4) Integrated energy bound for the perturbed gradient flow**:
if `c'(t) = -grad E(c(t)) + r(t)` and `‖r(t)‖ <= δ` on `[0,T]`, then for `0 <= s <= t <= T`
`E(c(t)) <= E(c(s)) + (δ^2/2)(t - s)`. -/
theorem perturbedFlow_energy_growth_bound
    (E : F → ℝ) (c : ℝ → F) (r : ℝ → F) (T δ : ℝ)
    (hEd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ E (c t))
    (hcd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ c t)
    (hflow : ∀ t ∈ Set.Icc 0 T, deriv c t = -gradient E (c t) + r t)
    (hr : ∀ t ∈ Set.Icc 0 T, ‖r t‖ ≤ δ) :
    ∀ s ∈ Set.Icc 0 T, ∀ t ∈ Set.Icc 0 T, s ≤ t →
      E (c t) ≤ E (c s) + δ ^ 2 / 2 * (t - s) := by
  intro s hs t ht hst
  have hDE : ∀ u ∈ Set.Icc 0 T,
      fderiv ℝ E (c u) (deriv c u)
        = ⟪gradient E (c u), (ContinuousLinearMap.id ℝ F) (deriv c u)⟫_ℝ + (fun _ ↦ (0 : ℝ)) u := by
    intro u _
    rw [← toDual_gradient, toDual_apply_apply]
    simp
  have := defectFlow_energy_growth_bound E (ContinuousLinearMap.id ℝ F) c
    (fun u => gradient E (c u)) r (fun _ => 0) T δ 0 hEd hcd hDE
    (by intro u hu; simpa using hflow u hu) hr (by intro u _; simp) s hs t ht hst
  simpa using this

/-- **(4) from time zero**: `E(c(t)) <= E(c(0)) + (δ^2/2) t` for `t ∈ [0,T]`
(the statement requested for the Spectral-PINN trajectory). -/
theorem perturbedFlow_energy_bound
    (E : F → ℝ) (c : ℝ → F) (r : ℝ → F) (T δ : ℝ)
    (hEd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ E (c t))
    (hcd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ c t)
    (hflow : ∀ t ∈ Set.Icc 0 T, deriv c t = -gradient E (c t) + r t)
    (hr : ∀ t ∈ Set.Icc 0 T, ‖r t‖ ≤ δ) :
    ∀ t ∈ Set.Icc 0 T, E (c t) ≤ E (c 0) + δ ^ 2 / 2 * t := by
  intro t ht
  have h0 : (0 : ℝ) ∈ Set.Icc 0 T := ⟨le_rfl, ht.1.trans ht.2⟩
  have := perturbedFlow_energy_growth_bound E c r T δ hEd hcd hflow hr 0 h0 t ht ht.1
  simpa using this

/-- **(4), sharp constant**: `E(c(t)) <= E(c(0)) + (δ^2/4) t` for `t ∈ [0,T]`. -/
theorem perturbedFlow_energy_bound_quarter
    (E : F → ℝ) (c : ℝ → F) (r : ℝ → F) (T δ : ℝ)
    (hEd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ E (c t))
    (hcd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ c t)
    (hflow : ∀ t ∈ Set.Icc 0 T, deriv c t = -gradient E (c t) + r t)
    (hr : ∀ t ∈ Set.Icc 0 T, ‖r t‖ ≤ δ) :
    ∀ t ∈ Set.Icc 0 T, E (c t) ≤ E (c 0) + δ ^ 2 / 4 * t := by
  have hDE : ∀ u ∈ Set.Icc 0 T,
      fderiv ℝ E (c u) (deriv c u)
        = ⟪gradient E (c u), (ContinuousLinearMap.id ℝ F) (deriv c u)⟫_ℝ + (fun _ ↦ (0 : ℝ)) u := by
    intro u _
    rw [← toDual_gradient, toDual_apply_apply]
    simp
  have := defectFlow_energy_bound_quarter E (ContinuousLinearMap.id ℝ F) c
    (fun u => gradient E (c u)) r (fun _ => 0) T δ 0 hEd hcd hDE
    (by intro u hu; simpa using hflow u hu) hr (by intro u _; simp)
  intro t ht
  simpa using this t ht

/-- Exact gradient flow (`δ = 0`): the integrated energy is non-increasing on `[0,T]`
(this recovers the integrated form of `GradientFlow.lean`). -/
theorem gradientFlow_energy_le_initial
    (E : F → ℝ) (c : ℝ → F) (T : ℝ)
    (hEd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ E (c t))
    (hcd : ∀ t ∈ Set.Icc 0 T, DifferentiableAt ℝ c t)
    (hflow : ∀ t ∈ Set.Icc 0 T, deriv c t = -gradient E (c t)) :
    ∀ t ∈ Set.Icc 0 T, E (c t) ≤ E (c 0) := by
  intro t ht
  have := perturbedFlow_energy_bound E c (fun _ => 0) T 0 hEd hcd
    (by intro u hu; simpa using hflow u hu) (by intro u _; simp) t ht
  simpa using this

end PerturbedIntegrated

/-! ## Section 5. Discrete-time statement (item (5)) -/

section Discrete

variable {H : Type*} [NormedAddCommGroup H]

/-- Telescoping: if `E_{n+1} - E_n <= a_n` for all `n`, then `E_N <= E_0 + Σ_{n<N} a_n`. -/
theorem energy_le_of_increment_le (E a : ℕ → ℝ) (h : ∀ n, E (n + 1) - E n ≤ a n) (N : ℕ) :
    E N ≤ E 0 + ∑ n ∈ Finset.range N, a n := by
  induction N with
  | zero => simp
  | succ N ih =>
    rw [Finset.sum_range_succ]
    linarith [h N]

/-- **(5) Discrete energy balance.**  If `E_{n+1} - E_n <= (1/2)‖r_n‖^2 - (1/2)‖d_n‖^2` for
all `n`, then `E_N + (1/2) Σ_{n<N} ‖d_n‖^2 <= E_0 + (1/2) Σ_{n<N} ‖r_n‖^2`. -/
theorem discrete_energy_balance (E : ℕ → ℝ) (r d : ℕ → H)
    (h : ∀ n, E (n + 1) - E n ≤ (1 / 2) * ‖r n‖ ^ 2 - (1 / 2) * ‖d n‖ ^ 2) (N : ℕ) :
    E N + (1 / 2) * ∑ n ∈ Finset.range N, ‖d n‖ ^ 2
      ≤ E 0 + (1 / 2) * ∑ n ∈ Finset.range N, ‖r n‖ ^ 2 := by
  have h1 := energy_le_of_increment_le E (fun n => (1 / 2) * ‖r n‖ ^ 2 - (1 / 2) * ‖d n‖ ^ 2) h N
  rw [Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum] at h1
  linarith

/-- **(5) Discrete energy bound**: `E_N <= E_0 + (1/2) Σ_{n<N} ‖r_n‖^2`. -/
theorem discrete_energy_bound (E : ℕ → ℝ) (r d : ℕ → H)
    (h : ∀ n, E (n + 1) - E n ≤ (1 / 2) * ‖r n‖ ^ 2 - (1 / 2) * ‖d n‖ ^ 2) (N : ℕ) :
    E N ≤ E 0 + (1 / 2) * ∑ n ∈ Finset.range N, ‖r n‖ ^ 2 := by
  have h1 := discrete_energy_balance E r d h N
  have h2 : 0 ≤ ∑ n ∈ Finset.range N, ‖d n‖ ^ 2 :=
    Finset.sum_nonneg (fun n _ => sq_nonneg _)
  linarith

/-- Uniform residual bound: if `‖r_n‖ <= δ` for `n < N` then `E_N <= E_0 + N δ^2 / 2`. -/
theorem discrete_energy_bound_uniform (E : ℕ → ℝ) (r d : ℕ → H) (δ : ℝ)
    (h : ∀ n, E (n + 1) - E n ≤ (1 / 2) * ‖r n‖ ^ 2 - (1 / 2) * ‖d n‖ ^ 2)
    (hr : ∀ n, ‖r n‖ ≤ δ) (N : ℕ) :
    E N ≤ E 0 + (N : ℝ) * δ ^ 2 / 2 := by
  have h1 := discrete_energy_bound E r d h N
  have h2 : ∑ n ∈ Finset.range N, ‖r n‖ ^ 2 ≤ ∑ n ∈ Finset.range N, δ ^ 2 :=
    Finset.sum_le_sum (fun n _ => pow_le_pow_left₀ (norm_nonneg _) (hr n) 2)
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul] at h2
  linarith

/-- Step-size weighted version (time steps `τ_n >= 0`):
if `E_{n+1} - E_n <= τ_n ((1/2)‖r_n‖^2 - (1/2)‖d_n‖^2)` then
`E_N + (1/2) Σ τ_n ‖d_n‖^2 <= E_0 + (1/2) Σ τ_n ‖r_n‖^2`.  With `Σ_{n<N} τ_n = t_N` and
`‖r_n‖ <= δ` this gives `E_N <= E_0 + (δ^2/2) t_N`, the discrete counterpart of (4). -/
theorem discrete_energy_balance_weighted (E τ : ℕ → ℝ) (r d : ℕ → H)
    (h : ∀ n, E (n + 1) - E n ≤ τ n * ((1 / 2) * ‖r n‖ ^ 2 - (1 / 2) * ‖d n‖ ^ 2)) (N : ℕ) :
    E N + (1 / 2) * ∑ n ∈ Finset.range N, τ n * ‖d n‖ ^ 2
      ≤ E 0 + (1 / 2) * ∑ n ∈ Finset.range N, τ n * ‖r n‖ ^ 2 := by
  have h1 := energy_le_of_increment_le E
    (fun n => τ n * ((1 / 2) * ‖r n‖ ^ 2 - (1 / 2) * ‖d n‖ ^ 2)) h N
  have h2 : ∑ n ∈ Finset.range N, τ n * ((1 / 2) * ‖r n‖ ^ 2 - (1 / 2) * ‖d n‖ ^ 2)
      = (1 / 2) * ∑ n ∈ Finset.range N, τ n * ‖r n‖ ^ 2
        - (1 / 2) * ∑ n ∈ Finset.range N, τ n * ‖d n‖ ^ 2 := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl (fun n _ => by ring)
  rw [h2] at h1
  linarith

/-- Weighted discrete bound with uniform residual bound and `Σ τ_n = t_N`:
`E_N <= E_0 + (δ^2/2) Σ_{n<N} τ_n`. -/
theorem discrete_energy_bound_weighted (E τ : ℕ → ℝ) (r d : ℕ → H) (δ : ℝ)
    (hτ : ∀ n, 0 ≤ τ n)
    (h : ∀ n, E (n + 1) - E n ≤ τ n * ((1 / 2) * ‖r n‖ ^ 2 - (1 / 2) * ‖d n‖ ^ 2))
    (hr : ∀ n, ‖r n‖ ≤ δ) (N : ℕ) :
    E N ≤ E 0 + δ ^ 2 / 2 * ∑ n ∈ Finset.range N, τ n := by
  have h1 := discrete_energy_balance_weighted E τ r d h N
  have h2 : 0 ≤ ∑ n ∈ Finset.range N, τ n * ‖d n‖ ^ 2 :=
    Finset.sum_nonneg (fun n _ => mul_nonneg (hτ n) (sq_nonneg _))
  have h3 : ∑ n ∈ Finset.range N, τ n * ‖r n‖ ^ 2 ≤ ∑ n ∈ Finset.range N, τ n * δ ^ 2 :=
    Finset.sum_le_sum (fun n _ =>
      mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (norm_nonneg _) (hr n) 2) (hτ n))
  rw [← Finset.sum_mul] at h3
  nlinarith

end Discrete

/-! ## Section 6. The spatial step for Allen-Cahn: integration by parts for polynomials

For the Spectral-PINN the state `u(., t) = Σ_k c_k(t) T_k` is a polynomial in `x`.  The
following two lemmas machine-check the spatial integration by parts that identifies, in the
splitting `DE(u)[w] = <G, w>_{L^2} + β(w)` of Section 3, the strong-form gradient
`G = -eps^2 u_xx + f(u)` (with `f(u) = u^3 - u`) and the defect functional
`β(w) = eps^2 [u_x w]_{x=-1}^{x=1}` (which vanishes iff `u_x(±1) = 0`).
(The differentiation of `E` under the integral sign is *not* formalized here; see
Steps 1-2 of the proof of Proposition A2 (Appendix B.2) of the accompanying article.) -/

section PolynomialIBP

open Polynomial

/-- Integration by parts on `[-1,1]` for real polynomials `p`, `w`:
`∫ p' w' = [p' w]_{-1}^{1} - ∫ p'' w`. -/
theorem ibp_polynomial (p w : ℝ[X]) :
    ∫ x in (-1:ℝ)..1, (derivative p).eval x * (derivative w).eval x
      = (derivative p).eval 1 * w.eval 1 - (derivative p).eval (-1) * w.eval (-1)
        - ∫ x in (-1:ℝ)..1, (derivative (derivative p)).eval x * w.eval x := by
  exact intervalIntegral.integral_mul_deriv_eq_deriv_mul (a := -1) (b := 1)
    (u := fun x => (derivative p).eval x) (v := fun x => w.eval x)
    (u' := fun x => (derivative (derivative p)).eval x) (v' := fun x => (derivative w).eval x)
    (fun x _ => (derivative p).hasDerivAt x) (fun x _ => w.hasDerivAt x)
    ((derivative (derivative p)).continuous.intervalIntegrable _ _)
    ((derivative w).continuous.intervalIntegrable _ _)

/-- **First variation of the Ginzburg-Landau energy, split into bulk and boundary parts.**
For polynomials `p` (state) and `w` (direction),
`∫ (eps^2 p' w' + (p^3 - p) w) = ∫ (-(eps^2 p'') + (p^3 - p)) w + eps^2 [p' w]_{-1}^{1}`.
Left side `= DE(p)[w]`; the first term on the right is `<G, w>_{L^2}` with
`G = -eps^2 p'' + f(p)`, the second is the defect `β(w)`. -/
theorem first_variation_split (ε : ℝ) (p w : ℝ[X]) :
    (∫ x in (-1:ℝ)..1, (ε ^ 2 * ((derivative p).eval x * (derivative w).eval x)
        + (p.eval x ^ 3 - p.eval x) * w.eval x))
      = (∫ x in (-1:ℝ)..1, (-(ε ^ 2 * (derivative (derivative p)).eval x)
        + (p.eval x ^ 3 - p.eval x)) * w.eval x)
        + ε ^ 2 * ((derivative p).eval 1 * w.eval 1 - (derivative p).eval (-1) * w.eval (-1)) := by
  have i1 : IntervalIntegrable (fun x : ℝ => (derivative p).eval x * (derivative w).eval x)
      MeasureTheory.volume (-1) 1 :=
    (((derivative p).continuous).mul ((derivative w).continuous)).intervalIntegrable _ _
  have i2 : IntervalIntegrable (fun x : ℝ => (p.eval x ^ 3 - p.eval x) * w.eval x)
      MeasureTheory.volume (-1) 1 :=
    ((((p.continuous).pow 3).sub p.continuous).mul w.continuous).intervalIntegrable _ _
  have i3 : IntervalIntegrable (fun x : ℝ => (derivative (derivative p)).eval x * w.eval x)
      MeasureTheory.volume (-1) 1 :=
    (((derivative (derivative p)).continuous).mul w.continuous).intervalIntegrable _ _
  have hibp := ibp_polynomial p w
  have e1 : (∫ x in (-1:ℝ)..1, (ε ^ 2 * ((derivative p).eval x * (derivative w).eval x)
        + (p.eval x ^ 3 - p.eval x) * w.eval x))
      = ε ^ 2 * (∫ x in (-1:ℝ)..1, (derivative p).eval x * (derivative w).eval x)
        + (∫ x in (-1:ℝ)..1, (p.eval x ^ 3 - p.eval x) * w.eval x) := by
    rw [intervalIntegral.integral_add (i1.const_mul (ε ^ 2)) i2,
      intervalIntegral.integral_const_mul]
  have e2 : (∫ x in (-1:ℝ)..1, (-(ε ^ 2 * (derivative (derivative p)).eval x)
        + (p.eval x ^ 3 - p.eval x)) * w.eval x)
      = (-(ε ^ 2)) * (∫ x in (-1:ℝ)..1, (derivative (derivative p)).eval x * w.eval x)
        + (∫ x in (-1:ℝ)..1, (p.eval x ^ 3 - p.eval x) * w.eval x) := by
    have h : ∀ x : ℝ, (-(ε ^ 2 * (derivative (derivative p)).eval x)
        + (p.eval x ^ 3 - p.eval x)) * w.eval x
        = (-(ε ^ 2)) * ((derivative (derivative p)).eval x * w.eval x)
          + (p.eval x ^ 3 - p.eval x) * w.eval x := fun x => by ring
    simp_rw [h]
    rw [intervalIntegral.integral_add (i3.const_mul (-(ε ^ 2))) i2,
      intervalIntegral.integral_const_mul]
  rw [e1, e2, hibp]
  ring

end PolynomialIBP

/-! ## Axiom audit -/

#print axioms inner_balance
#print axioms inner_balance_exact
#print axioms inner_balance_le
#print axioms hasDerivAt_energy_inner
#print axioms perturbedFlow_hasDerivAt_energy
#print axioms perturbedFlow_hasDerivAt_energy_gradient
#print axioms perturbedFlow_hasDerivAt_energy_exact
#print axioms perturbedFlow_deriv_energy_le
#print axioms perturbedFlow_deriv_energy_le_deriv
#print axioms perturbedFlow_exact_case
#print axioms perturbedFlow_deriv_energy_nonpos_of_residual_le
#print axioms perturbedFlow_deriv_energy_le_of_relative
#print axioms defectFlow_hasDerivAt_energy
#print axioms defectFlow_hasDerivAt_energy_functional
#print axioms defectFlow_hasDerivAt_energy_exact
#print axioms defectFlow_deriv_energy_le
#print axioms defectFlow_zero_defect
#print axioms energy_le_of_deriv_le
#print axioms defectFlow_energy_growth_bound
#print axioms defectFlow_energy_bound
#print axioms perturbedFlow_energy_growth_bound
#print axioms perturbedFlow_energy_bound
#print axioms gradientFlow_energy_le_initial
#print axioms energy_le_of_increment_le
#print axioms discrete_energy_balance
#print axioms discrete_energy_bound
#print axioms discrete_energy_bound_uniform
#print axioms discrete_energy_balance_weighted
#print axioms discrete_energy_bound_weighted
#print axioms ibp_polynomial
#print axioms first_variation_split
#print axioms inner_balance_quarter
#print axioms inner_balance_le_quarter
#print axioms defectFlow_deriv_energy_le_quarter
#print axioms defectFlow_energy_bound_quarter
#print axioms perturbedFlow_energy_bound_quarter
#print axioms energy_integral_le
#print axioms defectFlow_energy_integral_bound
#print axioms defectFree_energy_integral_bound
