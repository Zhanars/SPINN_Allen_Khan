# Verification of the first 16 new (`ref_`) bibliography entries

Source file: `D:\Videos\KazNPU\Symmetry_PINN_2026\manuscript\r1_src\80_refs.tex` (entries `ref_Xia2023` ... `ref_Li2026`, lines 49-64).
Date of check: 2026-10-12. Sources: Crossref REST API (`api.crossref.org/works/<DOI>`), arXiv API (`export.arxiv.org`), DataCite API (arXiv DOIs), iclr.cc (conference venue/dates).

Result: all 16 entries are consistent with the registry metadata (authors incl. order, initials, diacritics; title; venue; year; volume; pages/article number; DOI). No FIX needed. Two optional remarks are listed below the table (they are not errors).

| # | key | status | what was wrong | corrected line |
|---|-----|--------|----------------|----------------|
| 1 | ref_Xia2023 | OK | nothing (Crossref: Xia, Mingtao; Bottcher, Lucas; Chou, Tom; Mach. Learn.: Sci. Technol. 4(2), 025024, 2023) | -- |
| 2 | ref_Kutuk2025 | OK | nothing (Crossref: Kutuk, Mustafa; Yucel, Hamdullah; J. Comput. Sci. 87, 102577, 2025; diacritics match) | -- |
| 3 | ref_He2019 | OK | nothing (Crossref: He, Dongdong; Pan, Kejia; Numer. Methods Partial Differ. Equ. 35(3), 955-975; print year 2019, online 2018-12) | -- |
| 4 | ref_Du2024 | OK | nothing (arXiv 2312.05225 / DataCite: Du, Yiheng; Chalapathi, Nithin; Krishnapriyan, Aditi; title matches; "Accepted to ICLR 2024"; iclr.cc: 12th ICLR, Vienna, Austria, 7-11 May 2024; DOI is the arXiv DataCite DOI, not in Crossref by design) | -- |
| 5 | ref_Bizzi2025 | OK | nothing wrong; all 12 authors present in the right order (Marques, M. = Marcio; Mendon\c{c}a, L. = Leonardo Mendonca, cedilla correct). Note: Crossref gives only first names for Moreira/de Oliveira/Pereira; the middle initials L.M., C.J., J.M. in the .tex come from arXiv 2509.04966 (Leonardo M. Moreira, Christian Junior de Oliveira, Joao M. Pereira), so they are consistent but not present in the NeurIPS/Crossref record. Pages 60600-60630, NeurIPS 38, 2025, DOI verified | -- |
| 6 | ref_Zhang2024 | OK | nothing (Crossref: Zhang, Jiahao; Zhang, Shiheng; Shen, Jie; Lin, Guang; J. Comput. Phys. 498, 112638, 2024) | -- |
| 7 | ref_Lin2025 | OK | nothing (Crossref: Lin, Guang; Mou, Changhong; Zhang, Jiahao; J. Comput. Phys. 541, 114326, 2025-11) | -- |
| 8 | ref_Du2021b | OK | nothing (Crossref: Du, Qiang; Ju, Lili; Li, Xiao; Qiao, Zhonghua; SIAM Rev. 63(2), 317-359, 2021) | -- |
| 9 | ref_Tang2016 | OK | nothing (Crossref: Tang, Tao; Yang, Jiang; J. Comput. Math. 34(5), 451-461, 2016) | -- |
| 10 | ref_Grossmann2024 | OK | nothing (Crossref: Grossmann, Tamara G; Komorowska, Urszula Julia; Latz, Jonas; Schonlieb, Carola-Bibiane; IMA J. Appl. Math. 89(1), 143-174, 2024) | -- |
| 11 | ref_Tekriwal2021 | OK | nothing (Crossref: Tekriwal, Mohit; Duraisamy, Karthik; Jeannin, Jean-Baptiste; NASA Formal Methods, LNCS, pp. 322-339, 2021) | -- |
| 12 | ref_Zhao2021 | OK | nothing (Crossref: Zhao, Jia; Commun. Appl. Math. Comput. 3(2), 357-369, 2021) | -- |
| 13 | ref_Bruna2024 | OK | nothing (Crossref: Bruna, Joan; Peherstorfer, Benjamin; Vanden-Eijnden, Eric; J. Comput. Phys. 496, 112588, 2024) | -- |
| 14 | ref_Du2021a | OK | nothing (Crossref: Du, Yifan; Zaki, Tamer A.; Phys. Rev. E 104(4), 045303, 2021) | -- |
| 15 | ref_Cooley2025 | OK | nothing (arXiv 2410.03496 / DataCite: Cooley, Madison; Shankar, Varun; Kirby, Robert M.; Zhe, Shandian; arXiv journal_ref "Transactions on Machine Learning Research, 2025"; title matches; DOI is the arXiv DataCite DOI) | -- |
| 16 | ref_Li2026 | OK | nothing (Crossref: Li, Chunyan; Yu, Wenkai; Wang, Qi; Int. J. Numer. Anal. Model. 23(4), 576-604, published online 2026-05-29) | -- |

## Optional remarks (not errors, no change required)

- ref_Du2024, ref_Cooley2025: the DOIs are arXiv (DataCite) DOIs, so they do not resolve through Crossref (HTTP 404 there is expected); they do resolve in DataCite and the metadata matches. TMLR itself does not issue DOIs; if a direct link to the journal version is preferred for ref_Cooley2025, the arXiv record lists https://openreview.net/forum?id=KqRnsEMYLx (not added, to keep the file style).
- ref_Bizzi2025: if the middle initials are to follow the NeurIPS/Crossref record exactly, drop them (L., C., J.); the current form follows arXiv and is not wrong. The publisher (Neural Information Processing Systems Foundation) is not given in the entry; the other NeurIPS entry (ref_Krishnapriyan2021) uses the same format, so the style is consistent.
- ref_Tekriwal2021: Crossref lists no volume/editors for the chapter; the entry matches the style of the existing DeMouraUllrich2021 entry.

## Verified fields per entry (compact evidence)

- Xia2023: DOI 10.1088/2632-2153/acd0a1 -> title "Spectrally adapted physics-informed neural networks for solving unbounded domain problems", vol 4, issue 2, art. 025024, issued 2023-05-26.
- Kutuk2025: 10.1016/j.jocs.2025.102577 -> "Energy dissipation preserving physics informed neural network for Allen-Cahn equations", vol 87, art. 102577, 2025-05.
- He2019: 10.1002/num.22333 -> "Maximum norm error analysis of an unconditionally stable semi-implicit scheme for multi-dimensional Allen-Cahn equations", vol 35, issue 3, pp. 955-975, print 2019-05.
- Du2024: arXiv 2312.05225v2 -> "Neural Spectral Methods: Self-supervised learning in the spectral domain" (Yiheng Du, Nithin Chalapathi, Aditi Krishnapriyan), ICLR 2024.
- Bizzi2025: 10.52202/085713-1818 -> "Neuro-Spectral Architectures for Causal Physics-Informed Networks", Advances in Neural Information Processing Systems 38, pp. 60600-60630, 2025; arXiv 2509.04966 ("Accepted at NeurIPS 2025").
- Zhang2024: 10.1016/j.jcp.2023.112638 -> "Energy-dissipative evolutionary deep operator neural networks", vol 498, art. 112638, 2024-02.
- Lin2025: 10.1016/j.jcp.2025.114326 -> "Energy-dissipative evolutionary Kolmogorov-Arnold networks for complex PDE systems", vol 541, art. 114326, 2025-11.
- Du2021b: 10.1137/19M1243750 -> "Maximum Bound Principles for a Class of Semilinear Parabolic Equations and Exponential Time-Differencing Schemes", vol 63, issue 2, pp. 317-359, 2021.
- Tang2016: 10.4208/jcm.1603-m2014-0017 -> "Implicit-Explicit Scheme for the Allen-Cahn Equation Preserves the Maximum Principle", vol 34, issue 5, pp. 451-461, 2016-06.
- Grossmann2024: 10.1093/imamat/hxae011 -> "Can physics-informed neural networks beat the finite element method?", vol 89, issue 1, pp. 143-174, 2024.
- Tekriwal2021: 10.1007/978-3-030-76384-8_20 -> "A Formal Proof of the Lax Equivalence Theorem for Finite Difference Schemes", NASA Formal Methods (LNCS), pp. 322-339, 2021.
- Zhao2021: 10.1007/s42967-020-00105-2 -> "Discovering Phase Field Models from Image Data with the Pseudo-Spectral Physics Informed Neural Networks", vol 3, issue 2, pp. 357-369, 2021.
- Bruna2024: 10.1016/j.jcp.2023.112588 -> "Neural Galerkin schemes with active learning for high-dimensional evolution equations", vol 496, art. 112588, 2024-01.
- Du2021a: 10.1103/PhysRevE.104.045303 -> "Evolutional deep neural network", vol 104, issue 4, art. 045303, 2021-10-04.
- Cooley2025: arXiv 2410.03496v2 -> "Fourier PINNs: From Strong Boundary Conditions to Adaptive Fourier Bases", journal_ref TMLR 2025.
- Li2026: 10.4208/ijnam2026-1025 -> "Energy Dissipation Rate Guided Adaptive Sampling for Physics-Informed Neural Networks: Resolving Surface-Bulk Dynamics in Allen-Cahn Systems", vol 23, issue 4, pp. 576-604, 2026-05-29.
