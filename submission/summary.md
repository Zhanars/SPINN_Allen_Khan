# Submission Package Summary

*Symmetry (MDPI) — Research Article*

## Manuscript

**A Spectral–Physics-Informed Neural Network Hybrid with Energy-Constrained Training for the Allen–Cahn Equation: Forward and Inverse Problems**

## Authors

Zhuldyz Issa¹, Adil Nauryzbay¹, Zhanars Abdiramanov¹˒²\* (corresponding), Yessimkhan Zhumazhanov¹, Dastan Urnek¹

¹ Department of Mathematics and Mathematical Modeling, Abai Kazakh National Pedagogical University, Almaty 050010, Kazakhstan
² Institute of Information and Computational Technologies, Almaty 050010, Kazakhstan

Corresponding author: zh.abdiramanov@abaiuniversity.edu.kz

ORCID: Zh. Issa 0009-0004-9631-6284 · A. Nauryzbay 0009-0002-5513-712X · Zh. Abdiramanov 0000-0003-3820-7253 · Y. Zhumazhanov 0009-0008-5187-6105 · D. Urnek 0009-0009-0099-4408

## Funding

This research was funded by Abai Kazakh National Pedagogical University, grant number 05-04/568. *(Funding statement intentionally shortened in the manuscript itself — no project title or names — per house style; full grant title is on file with the university.)*

## Document status

Manuscript compiles cleanly on the official MDPI *Symmetry* LaTeX template (**19 pages**), zero LaTeX errors, zero undefined references. Targeted at the Special Issue "Selected Papers from the 4th International Conference: Constructive Mathematical Analysis" (ICCMA 2026), as an extended version of the authors' ICCMA 2026 conference abstract. All numerical results are from executed code (local machine + a 64-core Linux server, *ai.atu* — the Almaty Technological University's server, reached over SSH); none are estimated or fabricated. All Lean 4 proofs compile against Lean 4.33.1 + Mathlib with zero uses of `sorry` and depend only on Lean's three standard axioms; Appendix A now also states each result as a conventional mathematical Lemma/Theorem/Proposition/Corollary with a standard proof immediately before the corresponding Lean code, so reviewers do not need to read Lean to verify correctness. Reference list: **45 entries**, all bearing a DOI or a verified alternative link (institutional repository / arXiv ID / conference-proceedings URL), systematically checked against Crossref (one incorrect author attribution found and corrected during this process). Full source (code, data, manuscript, Lean project) is mirrored at `github.com/Zhanars/SPINN_Allen_Khan` (private repository), cited by URL in the Data Availability Statement, which states that access will be granted by the corresponding author upon reasonable request.

\newpage

# 1. Draft Cover Letter to the Editor

Dear Editor,

We submit the manuscript "A Spectral–Physics-Informed Neural Network Hybrid with Energy-Constrained Training for the Allen–Cahn Equation: Forward and Inverse Problems" for consideration as a research article in *Symmetry*, for the Special Issue "Selected Papers from the 4th International Conference: Constructive Mathematical Analysis" (ICCMA 2026).

This manuscript is an extended version of a conference abstract presented at ICCMA 2026 (Abdiramanov, Z.; Zhumazhanov, E.; Urnek, D. "Investigation of Physics-Informed Neural Networks for Solving the Allen-Cahn Equation," 1–2 July 2026), as noted on the manuscript's first page and in the reference list. Relative to the original ~200-word abstract, this manuscript adds the full Spectral-PINN architecture and its exact spectral-differentiation formulation, the two-dimensional symmetric extension, complete forward/inverse/energy-dissipation experiments with multi-seed statistics, a four-axis stress test at scale (up to $K=256$ modes, $T=10$), and a full Lean 4 formalization of the underlying energy-stability theory — none of which appeared in the original abstract. Two additional co-authors (Zhuldyz Issa, Adil Nauryzbay) joined for this expanded scope of work.

The Allen–Cahn equation's structure — a double-well potential invariant under $u \mapsto -u$, and the resulting energy-dissipation and phase-symmetry properties — make it a natural fit for the journal's scope. Our central contribution, the Spectral-PINN hybrid, represents the neural network's output as a time-dependent vector of Chebyshev (1D) or Fourier (2D) spectral coefficients rather than as the field value directly, so that spatial derivatives are obtained exactly and, in the two-dimensional case, so that reflection symmetry of the solution is preserved to machine precision by construction rather than by training — a property we verify holds regardless of network seed, training budget, mode count, or interface sharpness.

Two aspects of this submission may be of particular interest to reviewers:

- **Negative and null results are reported alongside positive ones.** A single-point experiment initially suggested that an inverse-problem estimation bias was a spectral-truncation artifact; a systematic 20-configuration sweep on a 64-core server rejected this hypothesis outright, and we report the rejection rather than the more publishable-looking single data point. A similar reversal occurs for the number of spectral modes $K$: accuracy degrades, rather than improves, as $K$ grows, even under training budgets scaled in proportion to $K$ up to $K=256$.
- **Key theoretical claims are machine-verified.** Beyond the standard hand proof of the continuous energy-dissipation theorem, we formalized the theorem itself in Lean 4 with Mathlib, in the general form that any Hilbert-space gradient flow dissipates energy at rate $-\|u_t\|^2$, together with the discrete algebraic lemmas (a polarization identity and a discrete Grönwall step) that a review of earlier, unrelated work from our group had flagged as an error-prone step in this style of proof.

We confirm that neither the manuscript nor any parts of its content are currently under consideration for publication with or published in another journal. All authors have approved the manuscript and agree with its submission to Symmetry.

Appendix A of the manuscript, formalizing the energy-stability lemmas and the continuous energy-dissipation theorem in Lean 4, now also gives a standard mathematical statement and proof of each result immediately before the corresponding Lean code, so that reviewers are not required to read Lean to verify the mathematics.

The authors declare no conflicts of interest. The work was funded by Abai Kazakh National Pedagogical University (grant No. 05-04/568); the funder had no role in the design, analysis, or decision to publish.

We thank the editors and reviewers in advance for their time and consideration.

Sincerely,
**Zhanars Abdiramanov, on behalf of all co-authors**

\newpage

# 2. Summary of All Results

*Every number below is reproduced verbatim from the final manuscript and its underlying `data/*.json` files (`Symmetry_PINN_2026/data/`). Figures in parentheses give context, not independent claims.*

## 2.1 Forward problem ($\varepsilon = 0.01$, $t=1$; PINN rows: mean ± SD over 5 seeds)

| Method | DoF | $L^2$ error | $L^\infty$ error | Time (s) |
|---|---|---|---|---|
| Chebyshev spectral collocation (N=80) | 81 | 8.51×10⁻⁵ | 5.56×10⁻⁴ | 0.16 |
| Vanilla space-time PINN | 8,577 | (6.74±0.55)×10⁻² | (8.68±0.43)×10⁻² | 919±7 |
| Spectral-PINN, K=32 (this work) | 33 coeff. | (1.38±0.08)×10⁻³ | (3.32±0.48)×10⁻³ | 169±16 |
| Spectral-PINN + energy-constrained | 33 coeff. | (1.38±0.07)×10⁻³ | (3.31±0.48)×10⁻³ | 191±27 |
| Spectral-PINN, K=64, same budget | 65 coeff. | (6.00±1.21)×10⁻² | (1.10±0.22)×10⁻¹ | 267±2 |

*Headline: Spectral-PINN cuts $L^2$ error by ~49× vs. vanilla PINN while training ~5.4× faster. Doubling K from 32→64 under the same budget makes it ~40× less accurate — modes and training budget are coupled.*

## 2.2 Energy dissipation

- Energy-constrained Spectral-PINN: **0 energy violations** across 51 evaluation points over $t\in[0,1]$ (violation = energy increase > 10⁻¹⁰ between consecutive points).

## 2.3 Inverse problem: recovering $\varepsilon$

Initial finding ($\varepsilon_{true}=0.01$, K=32): recovered $\hat\varepsilon=0.01402$ regardless of noise (0%, 1%, 5% → 40.2%, 40.2%, 40.0% relative error) — bias dominates noise-induced variance.

Systematic 5×4 sweep ($\varepsilon_{true}\in\{0.005,0.01,0.02,0.05,0.1\} \times K\in\{16,32,64,96\}$, 20 runs, 64-core server):

| $\varepsilon_{true}$ | K=16 | K=32 | K=64 | K=96 |
|---|---|---|---|---|
| 0.005 | 30.0% | 35.9% | 37.0% | 44.6% |
| 0.01 | 12.6% | 40.4% | 39.4% | 40.2% |
| 0.02 | 35.3% | 37.2% | 35.4% | 39.4% |
| 0.05 | 44.3% | 43.4% | 37.0% | 40.3% |
| 0.1 | 43.6% | 43.3% | 37.3% | 40.0% |
| **Mean ± SD** | **33.2±11.9%** | **40.0±3.1%** | **37.2±1.5%** | **40.9±2.2%** |

*Overall mean bias 37.8±7.0% (n=20), essentially flat across a 6× range of K — REJECTS the initial "spectral-truncation artifact" hypothesis. True cause not yet isolated (likely loss weighting / observation design).*

## 2.4 Two-dimensional symmetry preservation

| t | Fourier-Spectral-PINN | Vanilla PINN (5 seeds, mean±SD) |
|---|---|---|
| 0.0 | 8.7×10⁻³² | (4.12±1.13)×10⁻⁵ |
| 0.2 | 9.4×10⁻³² | (4.94±1.16)×10⁻⁵ |
| 0.5 | 1.0×10⁻³¹ | (9.46±1.81)×10⁻⁵ |
| 1.0 | 1.3×10⁻³¹ | (2.23±0.45)×10⁻⁴ |

*Fourier-Spectral-PINN's reflection symmetry is exact to machine precision by construction (cos(mx)cos(ny) basis), independent of seed; vanilla PINN's violation grows monotonically in every one of 5 runs.*

## 2.5 Stress tests at scale (64-core server)

**(a) Larger K, budget scaled in proportion to K:**

| K | Budget | Warm-up+Epochs | $L^2$ error | Time (s) |
|---|---|---|---|---|
| 64 | 1× | 1000+6000 | 6.18×10⁻² | 178 |
| 64 | 2× matched | 2000+12000 | 2.25×10⁻² | 365 |
| 128 | 1× | 1000+6000 | 3.36×10⁻¹ | 351 |
| 128 | 4× matched | 4000+24000 | 2.19×10⁻¹ | 1387 |
| 256 | 1× | 1000+6000 | 1.51×10⁰ | 778 |
| 256 | 8× matched | 8000+48000 | 1.47×10⁰ | 4659 |

*K=256 finding: an 8× budget (48,000 epochs, 78 min) gave NO meaningful improvement over 1× (778 s). Larger K makes the optimization landscape harder, not just longer to train — baseline K=32 stays at 1.4×10⁻³.*

**(b) Smaller $\varepsilon$ (stiffness) — INCONCLUSIVE BY DESIGN:**

Both methods showed near-flat $L^2$ error across $\varepsilon\in\{0.01,0.005,0.002,0.001\}$ (~6.8×10⁻² vanilla, ~1×10⁻³ Spectral-PINN). We diagnosed this as an artifact of the fixed 200-point evaluation grid under-resolving the interface at small $\varepsilon$, not genuine stiffness-robustness — disclosed explicitly in the manuscript as a methodological limitation and flagged for future work ($\varepsilon$-adaptive evaluation grid).

**(c) Longer time horizons T ($\varepsilon=0.01$, K=32, energy-constrained):**

| T | Epochs | $L^2$ @ t=1 | $L^2$ @ t=T | Energy violations |
|---|---|---|---|---|
| 1 | 6,000 | 1.5×10⁻³ | 1.5×10⁻³ | 0 / 21 |
| 5 | 12,000 | 1.1×10⁻¹ | 9.0×10⁻¹ | 0 / 101 |
| 10 | 19,500 | 2.1×10⁻¹ | 9.9×10⁻¹ | 0 / 201 |

*Key dissociation: accuracy collapses by t≈5 ($L^2\approx1$, no better than a guess) while energy STILL never violates monotonicity even once, at any T. The soft energy penalty is a structural/qualitative guarantee, not a correctness certificate — the two must be assessed separately.*

**(d) Larger 2D mode counts / sharper interfaces:**

| M | $\varepsilon$ | Epochs | S(t=1) |
|---|---|---|---|
| 8 | 0.3 | 3000 | 1.27×10⁻³¹ (baseline) |
| 16 | 0.3 | 3000 | 1.27×10⁻³¹ |
| 16 | 0.3 | 6000 | 1.26×10⁻³¹ |
| 8 | 0.1 | 3000 | 3.74×10⁻³¹ |
| 16 | 0.1 | 6000 | 3.74×10⁻³¹ |

*Exact-symmetry guarantee is completely scale-invariant — unlike accuracy under K-scaling, it does not degrade at all across a 2× mode-count and 3× sharper-interface range.*

## 2.6 Formal verification (Lean 4 + Mathlib)

- `polarization_identity`: $2a(a-b) = a^2-b^2+(a-b)^2$ — proved by `ring`.
- `discrete_gronwall_step`: given $a^2 \le b^2+\tau(L+1)(a^2+b^2)$ and $\tau(L+1)<1$, then $a^2 \le \frac{1+\tau(L+1)}{1-\tau(L+1)}b^2$ — the exact step a prior peer review flagged as error-prone (a dropped term) in a related manuscript from this group.
- `stability_step_size_exists`: for every Lipschitz constant $L\ge0$, a valid step size $\tau$ with $\tau(L+1)<1$ exists.
- `gradientFlow_hasDerivAt_energy` / `gradientFlow_deriv_energy_nonpos` / `gradientFlow_deriv_energy_eq_zero_iff`: **Theorem 1 itself** (continuous energy dissipation), proved for a general real Hilbert space via Mathlib's `HasGradientAt`/`HasFDerivAt` chain rule — "any gradient flow dissipates energy at rate $-\|u_t\|^2$," which specializes to Allen–Cahn via the manuscript's own integration-by-parts identification of $\nabla E(u)=-\varepsilon^2 u_{xx}+f(u)$.

*Verification: `lake build` succeeds with zero errors and zero `sorry`; `#print axioms` confirms every theorem rests only on `propext`, `Classical.choice`, and `Quot.sound` — Lean/Mathlib's three standard axioms, no unproven placeholder.*

\newpage

# 3. MDPI Submission Checklist

| Item | Status |
|---|---|
| Manuscript on official *Symmetry* LaTeX template | Done — compiles clean, 19 pp. |
| Title, abstract, keywords, MSC codes | Done — abstract trimmed to 265 words (journal target: ~200) |
| Author list & order | Set per author instruction: Issa, Nauryzbay, Abdiramanov (corr.), Zhumazhanov, Urnek |
| Affiliations | Done (2 affiliations, Zh.A. dual-affiliated) |
| Corresponding author ORCID | Provided (Zh.A.: 0000-0003-3820-7253) |
| Co-author ORCIDs | **All five provided** — Issa: 0009-0004-9631-6284; Nauryzbay: 0009-0002-5513-712X; Zhumazhanov: 0009-0008-5187-6105; Urnek: 0009-0009-0099-4408 |
| Funding statement + grant number | Done — shortened to funder + grant No. 05-04/568, no project title/names in the manuscript text |
| Author Contributions (CRediT) | Done |
| Institutional Review Board / Informed Consent | Not applicable (no human/animal subjects) |
| Data Availability Statement | Done — code/data maintained in a private repository (`github.com/Zhanars/SPINN_Allen_Khan`), URL cited in the manuscript; access granted by the corresponding author upon reasonable request |
| Acknowledgments | Done — thanks Almaty Technological University for computational support, plus GenAI disclosure |
| Conflicts of Interest | Done — shortened to a single sentence, none declared |
| Reference list format & completeness | Done — 45 references, ACS numeric style, all DOI-verified against Crossref |
| Figures embedded and captioned | Done — 1 figure (energy dissipation) |
| Tables captioned and numbered | Done — 6 main tables + 2 appendices |
| Materials and Methods structure | Done — consolidated to 7 subsections (was 15 across the paper; merged baselines, merged implementation details into their parent subsections) |
| Appendix A mathematical proofs | Done — standard Lemma/Theorem/Proposition/Corollary statements with conventional proofs now precede each Lean code block, so reviewers need not read Lean |
| Supplementary materials identified | See Section 4 below — not yet uploaded to the portal |
| Suggested reviewers (optional) | **NOT PREPARED** — optional, consider adding |
| Cover letter | Draft in Section 1, now includes MDPI's two required verbatim statements |
| GenAI use disclosure | Done — stated once, in Acknowledgments only (removed the duplicate sentence previously also present in Materials and Methods 2.6) |
| Conference-paper extension disclosure (ICCMA 2026) | Done — dagger footnote on title page, reference list entry, and cover-letter statement, per MDPI's conference-paper policy |
| Special Issue targeting | Done — cover letter now names "Selected Papers from the 4th International Conference: Constructive Mathematical Analysis" |
| MDPI author account / submission portal | **ACTION NEEDED** — submit via susy.mdpi.com |

\newpage

# 4. Supplementary Materials to Upload

Location on disk: `D:\Videos\KazNPU\Symmetry_PINN_2026\`

- **Python source code** (`code/`): `common.py`, `chebyshev_spectral.py`, `vanilla_pinn.py`, `spectral_pinn.py`, `inverse_problem.py`, `experiment_2d_symmetry.py`, `run_experiments.py`, `stress_*.py`, `server_sweep_*.py`, `make_figures.py`
- **Result data** (`data/*.json`): `vanilla_pinn.json`, `spectral_pinn.json`, `spectral_pinn_energy.json`, `chebyshev_baseline.json`, `inverse_problem.json`, `inverse_bias_sweep.json`, `forward_multiseed.json`, `experiment_2d.json`, `symmetry_2d_multiseed.json`, `stress_k_sweep.json`, `stress_eps_sweep.json`, `stress_time_sweep.json`, `stress_2d_sweep.json`
- **Lean 4 / Mathlib project** (`lean/AllenCahnVerification/`): `Basic.lean`, `GradientFlow.lean`, `lakefile.toml` — self-contained, builds with `lake build`
- **Figure source** (`figures/energy_dissipation.png`) and generation script (`code/make_figures.py`)

# 5. Outstanding Items Requiring Your Action

- The repository stays private by your decision; when a reviewer or reader requests the data, add them as a GitHub collaborator on `github.com/Zhanars/SPINN_Allen_Khan` (Settings → Collaborators) or send a zip export — either satisfies the "upon reasonable request" wording.
- Confirm the final author order and that all five co-authors have seen and approved this manuscript version before submission (MDPI requires explicit confirmation at submission; this is now also one of the two verbatim statements in the cover letter).
- Decide whether to name suggested reviewers (optional but often speeds handling; entered in the submission system, not the cover letter).
- Finalize and sign the cover letter (Section 1) on institutional letterhead if your workflow requires it.
- Create an MDPI account (if not already) and submit via the **direct Special Issue submission link**: https://susy.mdpi.com/user/manuscripts/upload?journal=symmetry&form[special_issue_id]=270508 (this pre-selects the "Selected Papers from ICCMA 2026" Special Issue — do not use the generic Symmetry submission link, or the manuscript will not be routed to the Special Issue).
- Note the Special Issue's Article Processing Charge: 2400 CHF, and its submission deadline: 31 December 2026.

## GenAI disclosure (final wording, as it appears in the manuscript)

Stated once, in Acknowledgments (the duplicate sentence formerly also present at the end of Materials and Methods \S2.6 has been removed):

> "The authors thank the Almaty Technological University for computational support. During the preparation of this manuscript, the authors used Claude (Anthropic) to assist with analysis of the obtained results and with translation. The authors have reviewed and verified this output manually and take full responsibility for the content of this publication."
