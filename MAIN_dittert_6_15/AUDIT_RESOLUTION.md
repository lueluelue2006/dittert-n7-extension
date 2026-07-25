# Resolution of the independent 74+1-item audit

This matrix records the disposition of the independent audit dated 25 July 2026. “Resolved” means the manuscript or executable package was changed and regression-tested. “Author action” means the anonymous review package cannot truthfully supply personal or externally assigned metadata.

## Mathematical P0 items

| ID | Status | Resolution |
|---|---|---|
| M1 | Resolved | The tangent first derivative `D Phi(U_n)[X]=0` is written before the negative-Hessian argument. |
| M2 | Resolved | All three product-energy lemmas are non-strict in general; strictness is recovered in applications from `D>0`, hence `y != 0`. |
| M3 | Resolved | The dimension-seven text explicitly states that `(1-t)^6` is a stronger sufficient condition than the naturally obtained `(1-t)^5`. |
| M4 | Resolved | The dimension-six side functional is defined as a genuine polynomial on the closed product of simplices. |
| M5 | Resolved | Both signs in the dimension-seven conditional-energy integral argument are written out. |
| M6 | Resolved | The two-level maximizer proof now gives `P_0`, the exact pair quadratic, relative-interior stationarity, feasible flat segments, strict decrease of the auxiliary functional, and the active-to-inactive one-sided derivative. |
| M7 | Resolved | `Cap(0)=0` is defined and every iterative extraction verifies positive relevant degree. |
| M8 | Resolved | The complementary-minor regularization treats both zero exponent and positive exponent/zero-minor possibilities. |
| M9 | Resolved | The historical statement for dimension four is weakened to full indecomposability; dimensions 16 and `n>=17` are identified as preprints; “only dimension 5 remains” is removed. |

## Program P0 items

| ID | Status | Resolution |
|---|---|---|
| C1 | Resolved | All validation `assert` statements in both author verifiers were replaced by explicit `require`/exception checks. Ordinary and `python -O` results are compared literally. |
| C2 | Resolved | All dimension-seven polynomial functions use immutable `GAMMA7`; a regression test checks invariance after the dimension-fifteen verifier runs. |
| C3 | Resolved | Both author verifiers are import-safe and use `main()` guards. |
| C4 | Resolved | The orchestration script records environment, optimization level, hashes, times, return codes, exact statistics, and prints success only after all stages pass. |

## Mathematical rigor and readability items

| ID | Status | Resolution |
|---|---|---|
| R1 | Resolved | A normalized block-variation identity is stated once and specialized in all stationarity lemmas. |
| R2 | Resolved | The dimension-six monotonicity derivative is displayed and bounded below by `pH>0`. |
| R3 | Resolved | Every one-sided branch states `1 <= b <= n-1` and explains why `b=0` gives zero deficit. |
| R4 | Resolved | The range `0<tau<E_0/p` is derived from `(p tau)^2 <= 15/319 < E_0^2`. |
| R5 | Resolved | The complementary-minor limiting argument is fully expanded. |
| R6 | Resolved | Column entropy uses the explicit regularization `B_epsilon=(1-epsilon)B+epsilon U_n` and `0 log 0=0`. |
| R7 | Resolved | The denominator-variable argument for zero capacity is written explicitly. |
| R8 | Resolved | The support-sensitive iteration checks positive degree at every stage. |
| R9 | Resolved | “Sharp” was replaced by “best-possible constant”; no classification of all equality cases is claimed. |
| R10 | Resolved | Core saturation and the `H_0` lower bound explicitly use equal total mass plus entrywise domination. |
| R11 | Resolved | Positive denominators are collected before every rational-to-polynomial reduction. |
| R12 | Resolved | The dimension-seven certificate is stated to hold on a larger closed square than the feasible set. |
| R13 | Resolved | The dimension-six triangular map is shown to cover exactly the enlarged closed line segment. |
| R14 | Resolved | Monotonic decrease of `M_n^(1)(L,c)` in `c` is stated before replacing `c_*` by `93/100`. |
| R15 | Resolved | The oriented binary-Pinsker integral explicitly records `q>=p` in the application. |
| R16 | Resolved | An exact Bernstein/de Casteljau certificate lemma is stated in the manuscript. |
| R17 | Resolved | The appendix now maps manuscript formulas to code functions, enumeration ranges, and case counts. |
| R18 | Resolved | “Unconditional” is defined as not assuming prior Dittert dimensions or unproved conjectures, while listing the standard mathematical tools still used. |

## Certificate architecture items

| ID | Status | Resolution |
|---|---|---|
| C5 | Resolved by independent implementations | The author generators are accompanied by separately written polynomial/Bernstein reconstructions. The release does not freeze a second copy of millions of coefficients; hashes and literal expected statistics bind the executable objects instead. |
| C6 | Resolved | All three independent checkers are included and run by default. |
| C7 | Resolved | Exact unit tests cover conversion, subdivision, degenerate degrees, endpoints through exact evaluation, import behavior, and immutable state. |
| C8 | Resolved | Failure messages include dimension/type/support or bridge pair, subdivision path, depth, and coefficient bounds where applicable. |
| C9 | Resolved | Depth, leaf, count, and completeness guards use explicit exceptions. |
| C10 | Resolved | Nonzero polynomials, rational domains, variable order, multidegrees, and family counts are checked. |
| C11 | Resolved | Dimension-dependent state is local or immutable; mutable global `g` is gone. |
| C12 | Resolved | Human and JSON outputs report all case counts, minima, and tree statistics. |
| C13 | Resolved | `verification_manifest.json` is generated automatically. |
| C14 | Resolved | Python `>=3.11` is declared in README and `pyproject.toml` and checked at runtime. |
| C15 | Resolved | A GitHub Actions workflow covers Python 3.11--3.13, ordinary/optimized equality, unit tests, full verification, and LaTeX. |
| C16 | Resolved in the trust model | The independent dimension-six univariate checker is standard-library-only; separate dimension-six and dimension-seven implementations reduce reliance on SymPy/NumPy. The author generator retains them only for exact expansion and object storage. |
| C17 | Resolved | `verification_output.txt` is generated by the orchestration script and includes environment, hashes, timing, statistics, and stage digests. |
| C18 | Resolved | Both author scripts carry a certificate version and an exact `main.tex` SHA-256 binding. |

## Writing and presentation items

| ID | Status | Resolution |
|---|---|---|
| W1 | Resolved | The title now says “An Exact Computer-Assisted Proof.” |
| W2 | Resolved | The abstract describes explicit rational checkers only after the optimization-mode defects were removed. |
| W3 | Resolved | The introduction contains a dimension/method/certificate-size overview table. |
| W4 | Resolved | The introduction gives an acyclic proof-dependency diagram. |
| W5 | Resolved | Very long exact fractions were moved to code/manifest; the text retains positivity and readable decimal scale. |
| W6 | Substantially resolved | The most confusing local group-size symbol was renamed, and the notation section declares section-local reuse. A wholesale renaming was avoided because it would create unnecessary transcription risk. |
| W7 | Resolved | The original core lemmas are expanded rather than delegated to “exactly as before.” |
| W8 | Resolved | The verification-status paragraph lists a representative, explicitly non-exhaustive set of standard tools. |
| W9 | Resolved | Publication status, historical context, and logical dependencies are separated. |
| W10 | Resolved | The redundant terminal restatement of the zero-block theorem was removed. |
| W11 | Resolved | Equation numbering uses `numberwithin`; manual counter resets are gone. |
| W12 | Resolved | Nonbreaking references, punctuation, and LaTeX warnings were cleaned. |
| W13 | Resolved | The appendix correctly says the second author verifier covers dimensions seven through fifteen. |
| W14 | Resolved | Knopp--Sinkhorn is cited in the historical discussion. |
| W15 | Resolved to available verified metadata | The 2026 items are uniformly labeled as dated preprints, and the audited DOI for the 2024 paper is included. Unverified DOI/MR/Zbl data were not invented. |
| W16 | Resolved | The final proof separates compactness, boundary exclusion, positive-maximizer classification, value, and equality uniqueness. |
| W17 | Resolved | Certificate magnitudes are given in readable decimal form where useful; proof decisions remain exact. |

## Submission and reproducibility items

| ID | Status | Resolution |
|---|---|---|
| P1 | Author action | Anonymous-review macros replace misleading `Author Name`. Real names, affiliations, email, ORCID, and corresponding author must be supplied by the authors when anonymity ends. |
| P2 | Resolved in anonymous form / author action | Code/data availability, conflicts, contributions, funding, and computational/language-assistance declarations are present; personal funding and acknowledgements remain for the authors. |
| P3 | Resolved | MIT code license and CC BY 4.0 manuscript-source notice are included. |
| P4 | Resolved | `SHA256SUMS` is generated automatically. |
| P5 | Author action | The manuscript and README require a tagged public repository and immutable DOI, which cannot be assigned inside this local anonymous revision. |
| P6 | Resolved in anonymous form | `CITATION.cff` is included with explicit fields to replace before public release. |
| P7 | Resolved | `Makefile` supplies `test`, `verify`, `pdf`, `hashes`, and `all`. |
| P8 | Resolved to reproducible version pins | `requirements.txt` and `pyproject.toml` pin exact versions and Python range. Platform-specific wheel hashes should be generated by the eventual public build service if desired. |
| P9 | Resolved | README records the verification workflow and the manifest records the measured time/memory for the actual release run. |
| P10 | Resolved | Verification output contains environment, source hashes, times, return codes, statistics, and output digests. |
| P11 | Resolved | The paper is A4, fonts are embedded, Unicode mapping is enabled, and the final PDF is rendered and inspected. Tagged-PDF conversion may still be applied if required by the target journal. |
| P12 | Resolved | The release ZIP excludes auxiliary TeX files, caches, and bytecode. |

## Additional image-reported item

| ID | Status | Resolution |
|---|---|---|
| +1 | Resolved | Both occurrences, and the exact verifier, now use the correct chain `15/319 < (7/32)^2 < (11/50)^2`. |

## Remaining non-mathematical author actions

The anonymous package is ready for referee circulation after regression. Before public or accepted release, the authors must disclose personal metadata, choose the final journal-compatible licensing language, create the public tagged repository/archive, and replace the repository/DOI placeholders.
