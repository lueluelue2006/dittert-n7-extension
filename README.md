# PINNED — Final main PDF (dimensions 6–15)

**Pinned main file (standalone PDF):**

**[Dittert_Conjecture_6_15_MAIN.pdf](./Dittert_Conjecture_6_15_MAIN.pdf)**

- Final unified manuscript: *Dittert’s Conjecture in Dimensions Six through Fifteen: An Exact Computer-Assisted Proof*
- Theorem range: **6 ≤ n ≤ 15**
- Status: **final post-audit package** (74+1-item independent audit resolved in manuscript + verifiers)
- Source tree: [`MAIN_dittert_6_15/main.pdf`](./MAIN_dittert_6_15/main.pdf)

Full package: [`MAIN_dittert_6_15/`](./MAIN_dittert_6_15/) · [zip](./Dittert_Conjecture_6_15_final_submission.zip) · [zip (alias)](./Dittert_Conjecture_6_15_submission.zip)

Audit matrix: [`AUDIT_RESOLUTION.md`](./AUDIT_RESOLUTION.md) · in-package [`MAIN_dittert_6_15/AUDIT_RESOLUTION.md`](./MAIN_dittert_6_15/AUDIT_RESOLUTION.md)

---

# Dittert conjecture — dimensions 6 through 15

## PRIMARY (final)

This repository’s primary materials are the **final** unified package for

> **Dittert’s Conjecture in Dimensions Six through Fifteen: An Exact Computer-Assisted Proof**  
> Main theorem: **\(6 \le n \le 15\)**.

| Path | What it is |
|------|------------|
| **`MAIN_dittert_6_15/`** | Full final submission tree (open this first) |
| **`Dittert_Conjecture_6_15_final_submission.zip`** | Same package as zip |
| `Dittert_Conjecture_6_15_submission.zip` | Alias of the final zip (keeps old path working) |
| `Dittert_Conjecture_6_15_MAIN.pdf` / `Dittert_Conjecture_6_15.pdf` | Pinned PDF copies at repo root |
| `Dittert_Conjecture_6_15.tex` | Final TeX copy at repo root |
| `AUDIT_RESOLUTION.md` | 74+1 audit disposition matrix |

### Contents of `MAIN_dittert_6_15/`

- `main.tex` / `main.pdf` — final paper (51 pages)
- `dittert_n6_certificates.py` — exact n=6 author verifier
- `dittert_n7_n15_certificates.py` — exact n=7–15 author verifier
- `independent_n6_fast_audit.py`, `independent_n6_univariate_audit.py`, `independent_n7_certificate_audit.py` — independent rebuilds
- `run_all_certificates.py` — full regression (author + independent + unit tests)
- `tests/test_exact_certificates.py` — exact arithmetic unit tests
- `verification_output.txt`, `verification_manifest.json`, `SHA256SUMS`
- `.github/workflows/verify.yml`, `Makefile`, licenses, `CITATION.cff`

### Run certificates

```bash
cd MAIN_dittert_6_15
python3 -m pip install -r requirements.txt
python3 run_all_certificates.py
# expected: status: passed
# expected: All Dittert certificates for dimensions 6 through 15 passed.
```

Recorded full run (package): ~78.6s, ordinary and `python -O` results identical, independent checkers and unit tests included.

### Source chat

Final revision package from ChatGPT conversation  
`https://chatgpt.com/c/6a63bcb4-5054-83eb-8fd9-a5198d9c4222`

---

## DRAFTS ONLY — everything else

**All paths under `drafts/` are earlier drafts / intermediate dumps.**  
They are **not** the primary version.

Includes:

- `drafts/pre_audit_main_20260724/` — previous primary before the final audit revision
- earlier n=6 approaches, separate n=7 notes, n=8–15-only manuscripts, older 7–15-only submission

---

## Status note

- This upload is the **final 6–15 main line** for public sharing.
- Author name / affiliation / ORCID / funding / acknowledgements may still be placeholders for anonymous or pre-submission use.
- Independent expert reading is still recommended before journal submission.
