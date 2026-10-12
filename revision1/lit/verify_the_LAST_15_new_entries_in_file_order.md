# Bibliography verification: last 15 new (ref_*) entries of 80_refs.tex

Source file: `D:\Videos\KazNPU\Symmetry_PINN_2026\manuscript\r1_src\80_refs.tex` (31 `ref_` entries on lines 49-79; this report covers lines 65-79, i.e. the LAST 15 in file order).

Method: Crossref REST API (`https://api.crossref.org/works/<DOI>`, User-Agent with mailto) for the 12 journal/proceedings DOIs; arXiv API (`export.arxiv.org/api/query?id_list=...`) plus DataCite (`api.datacite.org/dois/10.48550/arXiv...`) for the 3 arXiv-DOI items; for the 3 ML-proceedings items the PMLR / NeurIPS proceedings landing pages were additionally read to confirm page ranges. Compared: author list (family names, initials, order, diacritics), title, venue, year, volume, issue, pages/article number, DOI.

Result: 15 of 15 entries agree with the registry metadata. No FIX needed. One optional completeness improvement (pages for ref_Krishnapriyan2021), see below.

| # | key | status | what was wrong | corrected line |
|---|-----|--------|----------------|----------------|
| 1 | ref_Geng2025 | OK | nothing. Authors Geng Y., Burkovska O., Ju L., Zhang G., Gunzburger M.; title; CMAME 436 (Mar 2025), article 117721; DOI all match | unchanged |
| 2 | ref_Shen2016 | OK | nothing. Shen J., Tang T., Yang J.; Commun. Math. Sci. 14(6), 1517-1534, 2016; DOI 10.4310/CMS.2016.v14.n6.a3 match | unchanged |
| 3 | ref_Dong2021 | OK | nothing. Dong S., Ni N.; J. Comput. Phys. 435, 110242 (Jun 2021) match | unchanged |
| 4 | ref_Sukumar2022 | OK | nothing. Sukumar N., Srivastava A.; CMAME 389, 114333 (Feb 2022) match | unchanged |
| 5 | ref_Lu2021 | OK | nothing. Lu L., Pestourie R. (Raphael with diaeresis, initial only used), Yao W., Wang Z., Verdugo F., Johnson S.G.; SIAM J. Sci. Comput. 43(6), B1105-B1132, 2021 match | unchanged |
| 6 | ref_CardosoBihlo2025 | OK | nothing. Cardoso-Bihlo E., Bihlo A.; Neural Netw. 181, 106826 (Jan 2025) match | unchanged |
| 7 | ref_Krishnapriyan2021 | OK (optional addition) | all stated data correct (arXiv 2109.01050 authors/title match; arXiv journal-ref "NeurIPS 2021"; NeurIPS vol. 34). Only incomplete: no page range. NeurIPS proceedings page gives pp. 26548-26560 (citation_firstpage/lastpage meta). Not an error; the other proceedings entries in this block carry pages | unchanged is acceptable; optional line below |
| 8 | ref_Rathore2024 | OK | nothing. arXiv 2402.01868 and PMLR page: Rathore P., Lei W., Frangella Z., Lu L., Udell M.; ICML 2024 (41st), PMLR 235, 42159-42191 match | unchanged |
| 9 | ref_Wang2021a | OK | nothing. Wang S., Teng Y., Perdikaris P.; SIAM J. Sci. Comput. 43(5), A3055-A3081, 2021 match | unchanged |
| 10 | ref_Rahaman2019 | OK | nothing. arXiv 1806.08734 and PMLR page: 8 authors in order (Rahaman, Baratin, Arpit, Draxler, Lin, Hamprecht, Bengio, Courville); ICML 2019 (36th), PMLR 97, 5301-5310 match | unchanged |
| 11 | ref_McGreivy2024 | OK | nothing. McGreivy N., Hakim A.; Nat. Mach. Intell. 6(10), 1256-1269 (Sep 2024) match | unchanged |
| 12 | ref_Mishra2023 | OK | nothing. Mishra S., Molinaro R.; IMA J. Numer. Anal. 43(1), 1-43; print year 2023 (online Jan 2022) so year 2023 is the print/volume year, consistent | unchanged |
| 13 | ref_Boldo2013 | OK | nothing. Boldo S., Cl\'ement F., Filli\^atre J.-C., Mayero M., Melquiond G., Weis P. (diacritics correct); J. Autom. Reason. 50(4), 423-456, print Apr 2013 (online Aug 2012) match | unchanged |
| 14 | ref_Boldo2017 | OK | nothing. Boldo S., Cl\'ement F., Faissole F., Martin V., Mayero M.; Proc. 6th ACM SIGPLAN CPP, ACM, 2017 (Paris, Jan 2017), pp. 79-89 match | unchanged |
| 15 | ref_Immler2018 | OK | nothing. Immler F. (sole author); J. Autom. Reason. 61(1-4), 73-111, print Jun 2018 (online Jan 2018) match | unchanged |

## Optional improvement (not an error)

ref_Krishnapriyan2021: add the verified page range, in the same style as the ref_Bizzi2025 NeurIPS entry:

```
\bibitem{ref_Krishnapriyan2021} Krishnapriyan, A.S.; Gholami, A.; Zhe, S.; Kirby, R.M.; Mahoney, M.W. Characterizing possible failure modes in physics-informed neural networks. In \textit{Advances in Neural Information Processing Systems 34 (NeurIPS 2021)}; 2021; pp.\ 26548--26560. https://doi.org/10.48550/arXiv.2109.01050.
```

## Notes on what could not be verified

- NeurIPS 2021 has no publisher DOI in Crossref (the query.bibliographic search returned only different papers), so the arXiv DOI is the only resolvable DOI; it is registered with DataCite and its metadata (title, 5 authors, 2021) matches.
- ref_Boldo2017: Crossref gives the venue "Proceedings of the 6th ACM SIGPLAN Conference on Certified Programs and Proofs", event "CPP '17", Paris, France; start date 2017-01-16. The end date of the conference is not in the Crossref record, so the entry correctly omits conference dates.
- For the three arXiv-DOI items the DOI points to the arXiv preprint, not the proceedings version; this is how the file cites them and it is accurate for what the DOI resolves to.
