# Submission package: Dittert's Conjecture in Dimensions Seven through Fifteen

## Contents

- `main.tex` - self-contained LaTeX manuscript.
- `main.pdf` - compiled manuscript.
- `dittert_n7_n15_certificates.py` - exact rational and Bernstein certificate.
- `cover_letter_template.md` - editable cover-letter template.

## Reproduce the proof certificates

```bash
python3 dittert_n7_n15_certificates.py
```

The expected final line is:

```text
All exact rational and Bernstein certificates passed.
```

The script uses only Python's standard library and performs no floating-point comparison.

## Compile the manuscript

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

## Before submission

1. Replace `Author Name` in `main.tex` and add affiliation, postal address, e-mail, funding, and acknowledgements as appropriate.
2. Verify every bibliographic entry against the journal or arXiv record. The bibliography was inherited from the supplied draft and was not independently web-checked while preparing this package.
3. Select the target journal and apply its class file only after preserving a compiling copy of this version.
4. Upload `dittert_n7_n15_certificates.py` as supplementary material and mention the exact computer-assisted component in the cover letter.
5. Obtain an independent expert reading of the new dimension-seven argument, especially Theorems 3.10-3.11 and Lemmas 8.5-8.7. The exact code has been executed successfully, but the manuscript is not a formally machine-checked proof.

## Main new ingredient

The dimension-seven extension uses a complementary-minor entropy inequality and its single-bridge zero-block consequence. The final three two-variable alternatives are certified by exact dyadic Bernstein subdivision with `Fraction` arithmetic.
