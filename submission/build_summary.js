const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, AlignmentType, BorderStyle, LevelFormat, ExternalHyperlink,
  PageBreak, VerticalAlign,
} = require("docx");
const fs = require("fs");

const FULL_WIDTH = 9360; // DXA, ~6.5in usable width on A4 with 1in margins approx

function p(text, opts = {}) {
  const { bold, italic, size, spacingAfter = 160, alignment, color, heading } = opts;
  return new Paragraph({
    heading,
    alignment,
    spacing: { after: spacingAfter },
    children: [new TextRun({ text, bold, italics: italic, size, color })],
  });
}

function bullet(text, opts = {}) {
  return new Paragraph({
    numbering: { reference: "bullets", level: 0 },
    spacing: { after: 100 },
    children: [new TextRun({ text, ...opts })],
  });
}

function heading1(text) {
  return new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 320, after: 160 }, children: [new TextRun({ text })] });
}
function heading2(text) {
  return new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 240, after: 120 }, children: [new TextRun({ text })] });
}

function cell(text, opts = {}) {
  const { bold, width, shade, align = AlignmentType.LEFT } = opts;
  return new TableCell({
    width: { size: width, type: WidthType.DXA },
    shading: shade ? { type: ShadingType.CLEAR, fill: shade } : undefined,
    verticalAlign: VerticalAlign.CENTER,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: [new Paragraph({ alignment: align, children: [new TextRun({ text: String(text), bold, size: 20 })] })],
  });
}

function makeTable(headerRow, rows, colWidths) {
  const totalWidth = colWidths.reduce((a, b) => a + b, 0);
  return new Table({
    width: { size: totalWidth, type: WidthType.DXA },
    columnWidths: colWidths,
    rows: [
      new TableRow({
        tableHeader: true,
        children: headerRow.map((h, i) => cell(h, { bold: true, width: colWidths[i], shade: "D9E2F3", align: AlignmentType.CENTER })),
      }),
      ...rows.map(
        (r) =>
          new TableRow({
            children: r.map((v, i) => cell(v, { width: colWidths[i], align: i === 0 ? AlignmentType.LEFT : AlignmentType.CENTER })),
          })
      ),
    ],
  });
}

const doc = new Document({
  numbering: {
    config: [
      {
        reference: "bullets",
        levels: [{ level: 0, format: LevelFormat.BULLET, text: "\u2022", alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 420, hanging: 260 } } } }],
      },
    ],
  },
  sections: [
    {
      properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
      children: [
        // ===== Title block =====
        p("Submission Package Summary", { bold: true, size: 40, spacingAfter: 60 }),
        p("Symmetry (MDPI) — Research Article", { italic: true, size: 24, color: "555555", spacingAfter: 300 }),

        heading2("Manuscript"),
        p("A Spectral–Physics-Informed Neural Network Hybrid with Energy-Constrained Training for the Allen–Cahn Equation: Forward and Inverse Problems", { bold: true }),

        heading2("Authors"),
        p("Zhuldyz Issa \u00b9, Adil Nauryzbay \u00b9, Zhanars Abdiramanov \u00b9\u00b7\u00b2\u002A (corresponding), Esimkhan Zhumazhanov \u00b9, Dastan Urnek \u00b9"),
        p("\u00b9 Department of Mathematics and Mathematical Modeling, Abai Kazakh National Pedagogical University, Almaty 050010, Kazakhstan", { size: 20 }),
        p("\u00b2 Institute of Information and Computational Technologies, Almaty 050010, Kazakhstan", { size: 20 }),
        p("Corresponding author: zh.abdiramanov@abaiuniversity.edu.kz \u00b7 ORCID (Zh.A.): 0000-0003-3820-7253", { size: 20 }),

        heading2("Funding"),
        p("Abai Kazakh National Pedagogical University, intra-university grant competition for young researchers' scientific projects (2026), grant No. 05-04/568 — \u00abInvestigation of Physics-Informed Neural Network (PINN) methods for solving the Allen–Cahn equation\u00bb."),

        heading2("Document status"),
        p("Manuscript compiles cleanly on the official MDPI Symmetry LaTeX template (16 pages), zero LaTeX errors, zero undefined references. All numerical results are from executed code (local machine + a 64-core Linux server, ai.atu, reached over SSH); none are estimated or fabricated. All Lean 4 proofs compile against Lean 4.33.1 + Mathlib with zero uses of sorry and depend only on Lean's three standard axioms.", { italic: true }),

        new Paragraph({ children: [new PageBreak()] }),

        // ===== Cover letter =====
        heading1("1. Draft Cover Letter to the Editor"),
        p("Dear Editor,", { spacingAfter: 200 }),
        p("We submit the manuscript \u201cA Spectral\u2013Physics-Informed Neural Network Hybrid with Energy-Constrained Training for the Allen\u2013Cahn Equation: Forward and Inverse Problems\u201d for consideration as a research article in Symmetry.", { spacingAfter: 200 }),
        p("The Allen\u2013Cahn equation's structure \u2014 a double-well potential invariant under u \u21a6 \u2212u, and the resulting energy-dissipation and phase-symmetry properties \u2014 make it a natural fit for the journal's scope. Our central contribution, the Spectral-PINN hybrid, represents the neural network's output as a time-dependent vector of Chebyshev (1D) or Fourier (2D) spectral coefficients rather than as the field value directly, so that spatial derivatives are obtained exactly and, in the two-dimensional case, so that reflection symmetry of the solution is preserved to machine precision by construction rather than by training \u2014 a property we verify holds regardless of network seed, training budget, mode count, or interface sharpness.", { spacingAfter: 200 }),
        p("Two aspects of this submission may be of particular interest to reviewers:", { spacingAfter: 120 }),
        bullet("Negative and null results are reported alongside positive ones. A single-point experiment initially suggested that an inverse-problem estimation bias was a spectral-truncation artifact; a systematic 20-configuration sweep on a 64-core server rejected this hypothesis outright, and we report the rejection rather than the more publishable-looking single data point. A similar reversal occurs for the number of spectral modes K: accuracy degrades, rather than improves, as K grows, even under training budgets scaled in proportion to K up to K=256."),
        bullet("Key theoretical claims are machine-verified. Beyond the standard hand proof of the continuous energy-dissipation theorem, we formalized the theorem itself in Lean 4 with Mathlib, in the general form that any Hilbert-space gradient flow dissipates energy at rate \u2212\u2016u_t\u2016\u00b2, together with the discrete algebraic lemmas (a polarization identity and a discrete Gr\u00f6nwall step) that a review of earlier, unrelated work from our group had flagged as an error-prone step in this style of proof."),
        p("All authors have approved the manuscript and declare no conflicts of interest. The work was funded by Abai Kazakh National Pedagogical University (grant No. 05-04/568); the funder had no role in the design, analysis, or decision to publish. We confirm the manuscript is not under consideration elsewhere.", { spacingAfter: 200 }),
        p("We thank the editors and reviewers in advance for their time and consideration.", { spacingAfter: 200 }),
        p("Sincerely,", { spacingAfter: 60 }),
        p("Zhanars Abdiramanov, on behalf of all co-authors", { bold: true }),

        new Paragraph({ children: [new PageBreak()] }),

        // ===== Results summary =====
        heading1("2. Summary of All Results"),
        p("Every number below is reproduced verbatim from the final manuscript and its underlying data/*.json files (Symmetry_PINN_2026/data/). Figures in parentheses give context, not independent claims.", { italic: true, size: 20 }),

        heading2("2.1 Forward problem (\u03b5 = 0.01, t = 1; PINN rows: mean \u00b1 SD over 5 seeds)"),
        makeTable(
          ["Method", "DoF", "L2 error", "L\u221e error", "Time (s)"],
          [
            ["Chebyshev spectral collocation (N=80)", "81", "8.51\u00d710\u207b\u2075", "5.56\u00d710\u207b\u2074", "0.16"],
            ["Vanilla space-time PINN", "8,577", "(6.74\u00b10.55)\u00d710\u207b\u00b2", "(8.68\u00b10.43)\u00d710\u207b\u00b2", "919\u00b17"],
            ["Spectral-PINN, K=32 (this work)", "33 coeff.", "(1.38\u00b10.08)\u00d710\u207b\u00b3", "(3.32\u00b10.48)\u00d710\u207b\u00b3", "169\u00b116"],
            ["Spectral-PINN + energy-constrained", "33 coeff.", "(1.38\u00b10.07)\u00d710\u207b\u00b3", "(3.31\u00b10.48)\u00d710\u207b\u00b3", "191\u00b127"],
            ["Spectral-PINN, K=64, same budget", "65 coeff.", "(6.00\u00b11.21)\u00d710\u207b\u00b2", "(1.10\u00b10.22)\u00d710\u207b\u00b9", "267\u00b12"],
          ],
          [3200, 1400, 1800, 1800, 1360]
        ),
        p("Headline: Spectral-PINN cuts L2 error by ~49\u00d7 vs. vanilla PINN while training ~5.4\u00d7 faster. Doubling K from 32\u219264 under the same budget makes it ~40\u00d7 less accurate \u2014 modes and training budget are coupled.", { italic: true, size: 20, spacingAfter: 260 }),

        heading2("2.2 Energy dissipation"),
        bullet("Energy-constrained Spectral-PINN: 0 energy violations across 51 evaluation points over t\u2208[0,1] (violation = energy increase > 10\u207b\u00b9\u2070 between consecutive points)."),
        p("", { spacingAfter: 100 }),

        heading2("2.3 Inverse problem: recovering \u03b5"),
        p("Initial finding (\u03b5_true=0.01, K=32): recovered \u03b5\u0302=0.01402 regardless of noise (0%, 1%, 5% \u2192 40.2%, 40.2%, 40.0% relative error) \u2014 bias dominates noise-induced variance.", { size: 20 }),
        p("Systematic 5\u00d74 sweep (\u03b5_true\u2208{0.005,0.01,0.02,0.05,0.1} \u00d7 K\u2208{16,32,64,96}, 20 runs, 64-core server):", { size: 20, spacingAfter: 120 }),
        makeTable(
          ["\u03b5_true", "K=16", "K=32", "K=64", "K=96"],
          [
            ["0.005", "30.0%", "35.9%", "37.0%", "44.6%"],
            ["0.01", "12.6%", "40.4%", "39.4%", "40.2%"],
            ["0.02", "35.3%", "37.2%", "35.4%", "39.4%"],
            ["0.05", "44.3%", "43.4%", "37.0%", "40.3%"],
            ["0.1", "43.6%", "43.3%", "37.3%", "40.0%"],
            ["Mean \u00b1 SD", "33.2\u00b111.9%", "40.0\u00b13.1%", "37.2\u00b11.5%", "40.9\u00b12.2%"],
          ],
          [1800, 1890, 1890, 1890, 1890]
        ),
        p("Overall mean bias 37.8\u00b17.0% (n=20), essentially flat across a 6\u00d7 range of K \u2014 REJECTS the initial \u201cspectral-truncation artifact\u201d hypothesis. True cause not yet isolated (likely loss weighting / observation design).", { italic: true, size: 20, spacingAfter: 260 }),

        heading2("2.4 Two-dimensional symmetry preservation"),
        makeTable(
          ["t", "Fourier-Spectral-PINN", "Vanilla PINN (5 seeds, mean\u00b1SD)"],
          [
            ["0.0", "8.7\u00d710\u207b\u00b3\u00b2", "(4.12\u00b11.13)\u00d710\u207b\u2075"],
            ["0.2", "9.4\u00d710\u207b\u00b3\u00b2", "(4.94\u00b11.16)\u00d710\u207b\u2075"],
            ["0.5", "1.0\u00d710\u207b\u00b3\u00b9", "(9.46\u00b11.81)\u00d710\u207b\u2075"],
            ["1.0", "1.3\u00d710\u207b\u00b3\u00b9", "(2.23\u00b10.45)\u00d710\u207b\u2074"],
          ],
          [1200, 2600, 3200]
        ),
        p("Fourier-Spectral-PINN's reflection symmetry is exact to machine precision by construction (cos(mx)cos(ny) basis), independent of seed; vanilla PINN's violation grows monotonically in every one of 5 runs.", { italic: true, size: 20, spacingAfter: 260 }),

        heading2("2.5 Stress tests at scale (64-core server)"),
        p("(a) Larger K, budget scaled in proportion to K:", { bold: true, size: 22, spacingAfter: 100 }),
        makeTable(
          ["K", "Budget", "Warm-up+Epochs", "L2 error", "Time (s)"],
          [
            ["64", "1\u00d7", "1000+6000", "6.18\u00d710\u207b\u00b2", "178"],
            ["64", "2\u00d7 matched", "2000+12000", "2.25\u00d710\u207b\u00b2", "365"],
            ["128", "1\u00d7", "1000+6000", "3.36\u00d710\u207b\u00b9", "351"],
            ["128", "4\u00d7 matched", "4000+24000", "2.19\u00d710\u207b\u00b9", "1387"],
            ["256", "1\u00d7", "1000+6000", "1.51\u00d710\u2070", "778"],
            ["256", "8\u00d7 matched", "8000+48000", "1.47\u00d710\u2070", "4659"],
          ],
          [1400, 1800, 2200, 1800, 1560]
        ),
        p("K=256 finding: an 8\u00d7 budget (48,000 epochs, 78 min) gave NO meaningful improvement over 1\u00d7 (778 s). Larger K makes the optimization landscape harder, not just longer to train \u2014 baseline K=32 stays at 1.4\u00d710\u207b\u00b3.", { italic: true, size: 20, spacingAfter: 200 }),

        p("(b) Smaller \u03b5 (stiffness) \u2014 INCONCLUSIVE BY DESIGN:", { bold: true, size: 22, spacingAfter: 100 }),
        p("Both methods showed near-flat L2 error across \u03b5\u2208{0.01,0.005,0.002,0.001} (~6.8\u00d710\u207b\u00b2 vanilla, ~1\u00d710\u207b\u00b3 Spectral-PINN). We diagnosed this as an artifact of the fixed 200-point evaluation grid under-resolving the interface at small \u03b5, not genuine stiffness-robustness \u2014 disclosed explicitly in the manuscript as a methodological limitation and flagged for future work (\u03b5-adaptive evaluation grid).", { size: 20, spacingAfter: 200 }),

        p("(c) Longer time horizons T (\u03b5=0.01, K=32, energy-constrained):", { bold: true, size: 22, spacingAfter: 100 }),
        makeTable(
          ["T", "Epochs", "L2 @ t=1", "L2 @ t=T", "Energy violations"],
          [
            ["1", "6,000", "1.5\u00d710\u207b\u00b3", "1.5\u00d710\u207b\u00b3", "0 / 21"],
            ["5", "12,000", "1.1\u00d710\u207b\u00b9", "9.0\u00d710\u207b\u00b9", "0 / 101"],
            ["10", "19,500", "2.1\u00d710\u207b\u00b9", "9.9\u00d710\u207b\u00b9", "0 / 201"],
          ],
          [1400, 1800, 2000, 2000, 2160]
        ),
        p("Key dissociation: accuracy collapses by t\u22485 (L2\u22481, no better than a guess) while energy STILL never violates monotonicity even once, at any T. The soft energy penalty is a structural/qualitative guarantee, not a correctness certificate \u2014 the two must be assessed separately.", { italic: true, size: 20, spacingAfter: 200 }),

        p("(d) Larger 2D mode counts / sharper interfaces:", { bold: true, size: 22, spacingAfter: 100 }),
        makeTable(
          ["M", "\u03b5", "Epochs", "S(t=1)"],
          [
            ["8", "0.3", "3000", "1.27\u00d710\u207b\u00b3\u00b9 (baseline)"],
            ["16", "0.3", "3000", "1.27\u00d710\u207b\u00b3\u00b9"],
            ["16", "0.3", "6000", "1.26\u00d710\u207b\u00b3\u00b9"],
            ["8", "0.1", "3000", "3.74\u00d710\u207b\u00b3\u00b9"],
            ["16", "0.1", "6000", "3.74\u00d710\u207b\u00b3\u00b9"],
          ],
          [1400, 1400, 1800, 4760]
        ),
        p("Exact-symmetry guarantee is completely scale-invariant \u2014 unlike accuracy under K-scaling, it does not degrade at all across a 2\u00d7 mode-count and 3\u00d7 sharper-interface range.", { italic: true, size: 20, spacingAfter: 260 }),

        heading2("2.6 Formal verification (Lean 4 + Mathlib)"),
        bullet("polarization_identity: 2a(a\u2212b) = a\u00b2\u2212b\u00b2+(a\u2212b)\u00b2 \u2014 proved by ring."),
        bullet("discrete_gronwall_step: given a\u00b2 \u2264 b\u00b2+\u03c4(L+1)(a\u00b2+b\u00b2) and \u03c4(L+1)<1, then a\u00b2 \u2264 [(1+\u03c4(L+1))/(1\u2212\u03c4(L+1))]\u00b7b\u00b2 \u2014 the exact step a prior peer review flagged as error-prone (a dropped term) in a related manuscript from this group."),
        bullet("stability_step_size_exists: for every Lipschitz constant L\u22650, a valid step size \u03c4 with \u03c4(L+1)<1 exists."),
        bullet("gradientFlow_hasDerivAt_energy / gradientFlow_deriv_energy_nonpos / gradientFlow_deriv_energy_eq_zero_iff: Theorem 1 itself (continuous energy dissipation), proved for a general real Hilbert space via Mathlib's HasGradientAt/HasFDerivAt chain rule \u2014 \u201cany gradient flow dissipates energy at rate \u2212\u2016u_t\u2016\u00b2,\u201d which specializes to Allen\u2013Cahn via the manuscript's own integration-by-parts identification of \u2207E(u)=\u2212\u03b5\u00b2u_xx+f(u)."),
        p("Verification: lake build succeeds with zero errors and zero sorry; #print axioms confirms every theorem rests only on propext, Classical.choice, and Quot.sound \u2014 Lean/Mathlib's three standard axioms, no unproven placeholder.", { italic: true, size: 20, spacingAfter: 260 }),

        new Paragraph({ children: [new PageBreak()] }),

        // ===== Checklist =====
        heading1("3. MDPI Submission Checklist"),
        makeTable(
          ["Item", "Status"],
          [
            ["Manuscript on official Symmetry LaTeX template", "Done \u2014 compiles clean, 16 pp."],
            ["Title, abstract, keywords, MSC codes", "Done"],
            ["Author list & order", "Set per author instruction: Issa, Nauryzbay, Abdiramanov (corr.), Zhumazhanov, Urnek"],
            ["Affiliations", "Done (2 affiliations, Zh.A. dual-affiliated)"],
            ["Corresponding author ORCID", "Provided (Zh.A.: 0000-0003-3820-7253)"],
            ["Co-author ORCIDs", "NOT YET PROVIDED \u2014 action needed (see Section 4)"],
            ["Funding statement + grant number", "Done \u2014 No. 05-04/568"],
            ["Author Contributions (CRediT)", "Done"],
            ["Institutional Review Board / Informed Consent", "Not applicable (no human/animal subjects)"],
            ["Data Availability Statement", "Done \u2014 code/data \u201cupon reasonable request\u201d"],
            ["Conflicts of Interest", "Done \u2014 none declared"],
            ["Reference list format & completeness", "Done \u2014 20 references, ACS numeric style"],
            ["Figures embedded and captioned", "Done \u2014 1 figure (energy dissipation)"],
            ["Tables captioned and numbered", "Done \u2014 6 main tables + 2 appendices"],
            ["Supplementary materials identified", "See Section 5 below \u2014 not yet uploaded to the portal"],
            ["Suggested reviewers (optional)", "NOT PREPARED \u2014 optional, consider adding"],
            ["Cover letter", "Draft provided in Section 1 \u2014 needs your signature/finalization"],
            ["MDPI author account / submission portal", "ACTION NEEDED \u2014 submit via susy.mdpi.com"],
          ],
          [5200, 4160]
        ),

        new Paragraph({ children: [new PageBreak()] }),

        // ===== Supplementary materials =====
        heading1("4. Supplementary Materials to Upload"),
        p("Location on disk: D:\\Videos\\KazNPU\\Symmetry_PINN_2026\\", { size: 20, spacingAfter: 160 }),
        bullet("Python source code (code/): common.py, chebyshev_spectral.py, vanilla_pinn.py, spectral_pinn.py, inverse_problem.py, experiment_2d_symmetry.py, run_experiments.py, stress_*.py, server_sweep_*.py, make_figures.py"),
        bullet("Result data (data/*.json): vanilla_pinn.json, spectral_pinn.json, spectral_pinn_energy.json, chebyshev_baseline.json, inverse_problem.json, inverse_bias_sweep.json, forward_multiseed.json, experiment_2d.json, symmetry_2d_multiseed.json, stress_k_sweep.json, stress_eps_sweep.json, stress_time_sweep.json, stress_2d_sweep.json"),
        bullet("Lean 4 / Mathlib project (lean/AllenCahnVerification/): Basic.lean, GradientFlow.lean, lakefile.toml \u2014 self-contained, builds with `lake build`"),
        bullet("Figure source (figures/energy_dissipation.png) and generation script (code/make_figures.py)"),

        heading1("5. Outstanding Items Requiring Your Action"),
        bullet("Provide ORCID iDs for Zhuldyz Issa, Adil Nauryzbay, Esimkhan Zhumazhanov, and Dastan Urnek (MDPI requests these for all authors where available)."),
        bullet("Confirm the final author order and that all four co-authors have seen and approved this manuscript version before submission (MDPI requires explicit confirmation at submission)."),
        bullet("Decide whether to name suggested reviewers (optional but often speeds handling)."),
        bullet("Finalize and sign the cover letter (Section 1) on institutional letterhead if your workflow requires it."),
        bullet("Decide the actual hosting/access mechanism for the Data Availability Statement (currently \u201cupon reasonable request\u201d) \u2014 consider a public repository (e.g., Zenodo, GitHub) if you would prefer open code/data, which MDPI encourages."),
        bullet("Create an MDPI account (if not already) and submit via https://susy.mdpi.com/user/manuscripts/upload, selecting Symmetry as the target journal."),
      ],
    },
  ],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync("Symmetry_Submission_Summary.docx", buf);
  console.log("written");
});
