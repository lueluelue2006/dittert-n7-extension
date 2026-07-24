#!/usr/bin/env python3
"""Run all exact certificates accompanying the Dittert 6--15 paper."""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parent
for name in ("dittert_n6_certificates.py", "dittert_n7_n15_certificates.py"):
    print(f"\n=== Running {name} ===", flush=True)
    subprocess.run([sys.executable, str(root / name)], check=True)
print("\nAll Dittert certificates for dimensions 6 through 15 passed.")
