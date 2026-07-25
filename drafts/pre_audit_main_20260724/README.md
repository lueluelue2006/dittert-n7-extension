# Dittert's Conjecture in Dimensions Six through Fifteen

This submission package contains the manuscript and the complete exact-arithmetic verification code.

## Files

- `main.tex` - LaTeX source.
- `main.pdf` - compiled manuscript.
- `dittert_n6_certificates.py` - exact dimension-six Bernstein certificates (98 four-variable and 10 univariate polynomials).
- `dittert_n7_n15_certificates.py` - exact certificates for dimensions 7 through 15.
- `run_all_certificates.py` - runs both verifiers.
- `requirements.txt` - tested SymPy and NumPy versions for the dimension-six verifier.
- `cover_letter_template.md` - editable submission letter.

## Reproduction

With Python 3 available, install the two dependencies used only by the dimension-six polynomial generator:

```bash
python3 -m pip install -r requirements.txt
python3 run_all_certificates.py
```

A successful run ends with:

```text
All Dittert certificates for dimensions 6 through 15 passed.
```

The dimensions 7--15 script uses only the Python standard library. The dimension-six script uses SymPy only for exact rational polynomial expansion and NumPy only as an object-array container. All Bernstein coefficients, de Casteljau subdivisions, and comparisons are exact `fractions.Fraction` computations; neither script performs floating-point comparisons.

## Author metadata

Before submission, replace `Author Name` in `main.tex` and the PDF metadata. Add affiliation, email, ORCID, funding, acknowledgements, data/code repository information, and the target journal's required declarations.

## Bibliographic and independent review

The bibliography should be checked against MathSciNet, zbMATH, arXiv, or publisher records before submission. The new capacity, stationarity, smoothing, and certificate reductions should also receive independent expert review; successful code execution is not a substitute for mathematical refereeing.
