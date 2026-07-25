# Dittert's Conjecture in Dimensions Six through Fifteen

This archive accompanies the manuscript

> **Dittert's Conjecture in Dimensions Six through Fifteen: An Exact Computer-Assisted Proof**

It contains the LaTeX source, compiled PDF, exact author verifiers, independent reconstruction checkers, unit tests, and machine-readable verification records.

## Scope of the result

The manuscript proves that, for every integer `6 <= n <= 15`, the Dittert functional on the simplex of nonnegative `n x n` matrices with total entry sum `n` is uniquely maximized by the uniform matrix. The analytic proof reduces finitely many boundary configurations to exact rational inequalities. No floating-point comparison, numerical optimizer, random sampling, or tolerance is used in any certificate decision.

The proof is computer-assisted in the ordinary mathematical sense. It is not a Lean/Coq/Isabelle kernel formalization; human review of the analytic reductions remains necessary.

## Requirements

- Python **3.11 or newer**.
- For the dimension-six generator: `sympy==1.14.0` and `numpy==2.3.5`.
- For PDF compilation: a reasonably complete TeX Live installation with `pdflatex`.

Install the pinned Python dependencies with

```bash
python3 -m pip install -r requirements.txt
```

The dimensions 7--15 verifier and the independent dimension-six univariate checker use only the Python standard library. SymPy is used only to expand displayed rational formulas in the author dimension-six generator; NumPy is used only as an object-array container. All Bernstein coefficients and all correctness decisions are exact `fractions.Fraction` values.

## One-command verification

Run the complete regression suite:

```bash
python3 run_all_certificates.py
```

The default run performs all of the following:

1. author verifier for `n=6` in ordinary mode;
2. author verifier for `n=7,...,15` in ordinary mode;
3. both author verifiers again under `python -O` and literal comparison of the JSON results;
4. three independent reconstruction checkers;
5. exact unit tests for power-to-Bernstein conversion, de Casteljau subdivision, degenerate degrees, import behavior, and immutable dimension-seven state;
6. source hashing and manuscript-hash binding.

A successful run writes:

- `verification_manifest.json` -- machine-readable environment, hashes, statistics, and stage results;
- `verification_output.txt` -- human-readable verification record;
- `SHA256SUMS` -- hashes of the clean submission files.

It ends with

```text
All Dittert certificates for dimensions 6 through 15 passed.
```

All mathematical checks use explicit exceptions rather than Python `assert`, so optimized execution cannot silently remove a certificate check.

Useful shorter commands are

```bash
python3 run_all_certificates.py --skip-independent --skip-optimized
python3 dittert_n6_certificates.py --json
python3 dittert_n7_n15_certificates.py --json
python3 -m unittest discover -s tests -v
```

## Make targets

```bash
make test       # exact unit tests
make verify     # full author + optimized + independent verification
make pdf        # compile main.tex twice
make hashes     # refresh SHA256SUMS
make all        # tests, verification, and PDF
```

## Principal files

- `main.tex`, `main.pdf` -- manuscript source and compiled paper.
- `dittert_n6_certificates.py` -- exact dimension-six generator and checker: 98 four-variable and 10 univariate polynomials.
- `dittert_n7_n15_certificates.py` -- standard-library exact certificates for dimensions 7 through 15.
- `independent_n6_fast_audit.py` -- independent reconstruction of the 98 dimension-six four-variable certificates.
- `independent_n6_univariate_audit.py` -- independent standard-library reconstruction of the ten dimension-six univariate certificates.
- `independent_n7_certificate_audit.py` -- independent symbolic and Bernstein reconstruction of the three dimension-seven single-bridge trees.
- `tests/test_exact_certificates.py` -- exact unit tests.
- `run_all_certificates.py` -- orchestration, environment capture, cross-mode regression, hashing, and manifest generation.
- `AUDIT_RESOLUTION.md` -- response matrix for the independent 74+1-item audit.
- `.github/workflows/verify.yml` -- continuous-integration template for Python 3.11--3.13 and LaTeX.

## Reproducibility and trust boundary

The author dimension-six generator and its checker are in one script, but the package also contains independently written reconstruction programs. These use different sparse-polynomial and Bernstein implementations and check the published case counts, subdivision statistics, depths, and minimum terminal coefficients. This reduces the risk of a shared implementation error, though it does not replace mathematical refereeing or proof-assistant formalization.

The exact release-run statistics are stored in `verification_manifest.json`. In the reference environment (CPython 3.13.5, SymPy 1.14.0, NumPy 2.3.5, Linux x86-64), a full ordinary/optimized/independent/test run takes about **81 seconds** and uses about **158 MB** of child-process peak resident memory. Timing and memory are informational only and play no role in correctness.

## Building the PDF

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The archived manuscript is A4, uses embedded Type 1 fonts, and is checked for undefined references, overfull/underfull boxes, and visual clipping before release.

## Anonymous-review metadata

The current source uses anonymous-review macros. Before a non-anonymous submission or accepted version, replace the author, affiliation, email, ORCID, funding, acknowledgements, contribution statement, repository URL, release tag, and archival DOI. `CITATION.cff` and `cover_letter_template.md` contain corresponding placeholders.

## Licensing

- Supplemental code: MIT License (`LICENSE-CODE`).
- Manuscript source: CC BY 4.0 notice (`LICENSE-MANUSCRIPT`), subject to the target journal's policy.

## Archival release

For final submission, create a tagged public repository release and an immutable archive (for example, a DOI-bearing repository deposit), then record its URL, tag/commit, and DOI in the manuscript, `CITATION.cff`, and cover letter.
