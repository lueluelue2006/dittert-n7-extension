# Dittert intermediate drafts (n=6, 7, 8–15)

**Status: DRAFT only. Not a final version.**

This repository is a public dump of materials that are **mostly confirmed on the author’s side** (local notes + certificate scripts that run with exact rational arithmetic). It is **not** a camera-ready paper, **not** a claim of journal-ready verification, and **not** a final unified write-up.

A later pass will **unify formatting, structure, and packaging**. Treat everything here as **working drafts** for feedback and timestamping.

## What is here

### n=6 (two different reduction styles)

Both claim a route to Dittert for dimension 6. They share the same broad Hall / core-entropy framework but are **not** two fully independent theories; they are two technical branches.

| Folder | Note | Certificates |
|--------|------|----------------|
| `n6/approach_A_fourvar/` | `dittert_n6_proof.tex` — pairwise smoothing → **four-variable** Bernstein family | `dittert_n6_certificates.py` (uses SymPy + NumPy; comparisons are exact `Fraction`) |
| `n6/approach_B_penalized/` | `dittert_n6_extension_note.tex` — **penalized** core–entropy criterion → **bivariate** Bernstein | `dittert_n6_certificates.py` (stdlib + `Fraction` only) |

### n=7

| File | Description |
|------|-------------|
| `dittert_n7_extension_note.pdf` / `.tex` | Extension note |
| `dittert_n7_certificates.py` | Exact certificates |

### n=8–15

| File | Description |
|------|-------------|
| `dittert_n8_n15_manuscript.tex` | Companion manuscript source |
| `dittert_n8_n15_certificates.py` | Certificates |

The n=6 and n=7 notes inherit dimension-free infrastructure from the n=8–15 line; they are **not** self-contained from first principles without that package (or equivalent).

## How to run certificates

```bash
# n=7
python3 dittert_n7_certificates.py

# n=8–15
python3 dittert_n8_n15_certificates.py

# n=6, approach B (stdlib only)
python3 n6/approach_B_penalized/dittert_n6_certificates.py

# n=6, approach A (needs sympy, numpy)
python3 n6/approach_A_fourvar/dittert_n6_certificates.py
```

## Disclaimer

- Draft upload for sharing and backup.
- **No** assertion that every lemma has been independently refereed.
- **Final** manuscript layout, single proof tree, and submission package will be prepared separately.
- Comments welcome; please treat claims as provisional until the unified version appears.
