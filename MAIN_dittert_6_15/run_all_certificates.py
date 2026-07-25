#!/usr/bin/env python3
"""Run, cross-check, and document all exact Dittert certificates.

The default full run executes both author verifiers in ordinary and optimized
Python modes, runs the independent reconstruction checkers and unit tests,
checks the embedded manuscript hash, and writes ``verification_manifest.json``,
``verification_output.txt``, and ``SHA256SUMS``.  Success is printed only after
all subprocesses and all literal statistic comparisons have completed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "1.0.0"
MIN_PYTHON = (3, 11)
ROOT = Path(__file__).resolve().parent
MANIFEST_PATH = ROOT / "verification_manifest.json"
OUTPUT_PATH = ROOT / "verification_output.txt"
HASHES_PATH = ROOT / "SHA256SUMS"

AUTHOR_SCRIPTS = (
    "dittert_n6_certificates.py",
    "dittert_n7_n15_certificates.py",
)
INDEPENDENT_SCRIPTS = (
    "independent_n6_fast_audit.py",
    "independent_n6_univariate_audit.py",
    "independent_n7_certificate_audit.py",
)
SOURCE_FILES = (
    "main.tex",
    "dittert_n6_certificates.py",
    "dittert_n7_n15_certificates.py",
    "run_all_certificates.py",
    "independent_n6_fast_audit.py",
    "independent_n6_univariate_audit.py",
    "independent_n7_certificate_audit.py",
    "tests/test_exact_certificates.py",
    "README.md",
    "requirements.txt",
    "pyproject.toml",
    "Makefile",
)


class VerificationError(RuntimeError):
    """Raised when a subprocess or consistency check fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_version(name: str) -> str | None:
    try:
        module = __import__(name)
    except Exception:
        return None
    return str(getattr(module, "__version__", "unknown"))


def environment() -> dict[str, Any]:
    return {
        "python": sys.version.replace("\n", " "),
        "python_executable": sys.executable,
        "python_implementation": platform.python_implementation(),
        "optimize_flag": sys.flags.optimize,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "numpy": package_version("numpy"),
        "sympy": package_version("sympy"),
        "timezone": "UTC",
    }


def existing_source_hashes() -> dict[str, str]:
    hashes: dict[str, str] = {}
    for name in SOURCE_FILES:
        path = ROOT / name
        if path.exists():
            hashes[name] = sha256_file(path)
    return hashes


def run_stage(
    name: str,
    command: list[str],
    *,
    expect: str | None = None,
    display_output: bool = True,
) -> dict[str, Any]:
    print(f"\n=== {name} ===", flush=True)
    print("$ " + " ".join(command), flush=True)
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        env={**os.environ, "PYTHONHASHSEED": "0"},
    )
    elapsed = time.perf_counter() - started
    if display_output and completed.stdout:
        print(completed.stdout.rstrip(), flush=True)
    if completed.stderr:
        print(completed.stderr.rstrip(), file=sys.stderr, flush=True)
    require(
        completed.returncode == 0,
        f"stage {name!r} failed with exit code {completed.returncode}",
    )
    if expect is not None:
        combined_output = completed.stdout + "\n" + completed.stderr
        require(expect in combined_output, f"stage {name!r} missing sentinel {expect!r}")
    return {
        "name": name,
        "command": command,
        "elapsed_seconds": round(elapsed, 6),
        "returncode": completed.returncode,
        "stdout_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr.encode()).hexdigest(),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def parse_json_stage(stage: dict[str, Any]) -> dict[str, Any]:
    try:
        value = json.loads(stage["stdout"])
    except json.JSONDecodeError as exc:
        raise VerificationError(f"stage {stage['name']!r} did not return valid JSON") from exc
    require(value.get("status") == "passed", f"stage {stage['name']!r} status={value.get('status')}")
    return value


def concise_stage_record(stage: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in stage.items() if key not in {"stdout", "stderr"}}


def author_summary(n6: dict[str, Any], n7: dict[str, Any]) -> list[str]:
    lines = [
        "Author-verifier statistics:",
        "  n=6: 98 four-variable polynomials and 10 univariate polynomials.",
    ]
    for pair, stats in n6["two_sided"].items():
        lines.append(
            f"    type ({pair}): polynomials={stats['polynomials']}, "
            f"nodes={stats['nodes']}, depth={stats['depth']}, minimum={stats['minimum']}"
        )
    lines.append(
        f"  n=7: ordinary={n7['n7']['ordinary_cases']}, mixed={n7['n7']['mixed_cases']}."
    )
    for pair, stats in n7["n7"]["single_bridge_trees"].items():
        lines.append(
            f"    bridge ({pair}): nodes={stats['nodes']}, "
            f"leaves={stats['direct_leaves']}+{stats['entropy_leaves']}, "
            f"depth={stats['max_depth']}"
        )
    for n, stats in n7["n8_n10"].items():
        lines.append(f"  n={n}: {stats['cases']} second-order Hall cases.")
    for n in ("11", "12", "13", "14"):
        lines.append(f"  n={n}: {n7['n11_n14'][n]['cases']} first-order Hall cases.")
    lines.append("  n=15: two rational margins, six integer gaps, and two coarse gaps.")
    return lines


def write_output(
    *,
    env: dict[str, Any],
    source_hashes: dict[str, str],
    stages: list[dict[str, Any]],
    n6: dict[str, Any],
    n7: dict[str, Any],
    total_elapsed: float,
    maxrss_kb: int | None,
) -> None:
    lines = [
        "Dittert dimensions 6--15 exact verification record",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}",
        "",
        "== Environment ==",
    ]
    for key, value in env.items():
        lines.append(f"{key}: {value}")
    lines.extend(["", "== Source SHA-256 =="])
    for name, digest in sorted(source_hashes.items()):
        lines.append(f"{digest}  {name}")
    lines.extend(["", "== Exact certificate statistics =="])
    lines.extend(author_summary(n6, n7))
    lines.extend(["", "== Stages =="])
    for stage in stages:
        lines.append(
            f"{stage['name']}: exit={stage['returncode']} "
            f"elapsed={stage['elapsed_seconds']:.6f}s stdout_sha256={stage['stdout_sha256']}"
        )
        stdout_lines = [line for line in stage["stdout"].splitlines() if line.strip()]
        if stdout_lines and not stage["name"].startswith("author"):
            for line in stdout_lines[-12:]:
                lines.append(f"  {line}")
    lines.extend([
        "",
        f"total_elapsed_seconds: {total_elapsed:.6f}",
        f"children_maxrss_kb: {maxrss_kb}",
        "status: passed",
        "All Dittert certificates for dimensions 6 through 15 passed.",
    ])
    OUTPUT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_sha256sums() -> None:
    excluded = {
        HASHES_PATH.name,
        "main.aux",
        "main.log",
        "main.out",
        "main.toc",
    }
    entries: list[tuple[str, str]] = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT).as_posix()
        if any(part in {"__pycache__", ".git"} for part in path.parts):
            continue
        if path.name in excluded or path.suffix in {".pyc"}:
            continue
        entries.append((relative, sha256_file(path)))
    HASHES_PATH.write_text(
        "".join(f"{digest}  {name}\n" for name, digest in entries),
        encoding="utf-8",
    )


def verify_full(*, skip_independent: bool, skip_optimized: bool, skip_tests: bool) -> dict[str, Any]:
    require(sys.version_info >= MIN_PYTHON, f"Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]} or newer is required")
    for name in AUTHOR_SCRIPTS + INDEPENDENT_SCRIPTS:
        require((ROOT / name).is_file(), f"missing verifier: {name}")

    started = time.perf_counter()
    env = environment()
    source_hashes = existing_source_hashes()
    manuscript_hash = source_hashes.get("main.tex")
    require(manuscript_hash is not None, "main.tex hash unavailable")
    stages: list[dict[str, Any]] = []

    author_n6_normal_stage = run_stage(
        "author n=6 (ordinary mode)",
        [sys.executable, str(ROOT / AUTHOR_SCRIPTS[0]), "--json"],
        display_output=False,
    )
    stages.append(author_n6_normal_stage)
    author_n6 = parse_json_stage(author_n6_normal_stage)

    author_n7_normal_stage = run_stage(
        "author n=7--15 (ordinary mode)",
        [sys.executable, str(ROOT / AUTHOR_SCRIPTS[1]), "--json"],
        display_output=False,
    )
    stages.append(author_n7_normal_stage)
    author_n7 = parse_json_stage(author_n7_normal_stage)

    require(
        author_n6.get("manuscript_sha256") == manuscript_hash,
        "n=6 verifier is not bound to the current main.tex hash",
    )
    require(
        author_n7.get("manuscript_sha256") == manuscript_hash,
        "n=7--15 verifier is not bound to the current main.tex hash",
    )

    optimized_match = None
    if not skip_optimized:
        n6_opt_stage = run_stage(
            "author n=6 (optimized mode)",
            [sys.executable, "-O", str(ROOT / AUTHOR_SCRIPTS[0]), "--json"],
            display_output=False,
        )
        stages.append(n6_opt_stage)
        n7_opt_stage = run_stage(
            "author n=7--15 (optimized mode)",
            [sys.executable, "-O", str(ROOT / AUTHOR_SCRIPTS[1]), "--json"],
            display_output=False,
        )
        stages.append(n7_opt_stage)
        n6_opt = parse_json_stage(n6_opt_stage)
        n7_opt = parse_json_stage(n7_opt_stage)
        require(n6_opt == author_n6, "n=6 ordinary and optimized JSON results differ")
        require(n7_opt == author_n7, "n=7--15 ordinary and optimized JSON results differ")
        optimized_match = True

    if not skip_independent:
        expected_sentinels = {
            "independent_n6_fast_audit.py": "ALL 98 N=6 FOUR-VARIABLE CERTIFICATES INDEPENDENTLY MATCHED",
            "independent_n6_univariate_audit.py": "ALL 10 N=6 UNIVARIATE CERTIFICATES INDEPENDENTLY MATCHED",
            "independent_n7_certificate_audit.py": "ALL N7 INDEPENDENT CHECKS PASSED",
        }
        for filename in INDEPENDENT_SCRIPTS:
            stages.append(
                run_stage(
                    f"independent checker: {filename}",
                    [sys.executable, str(ROOT / filename)],
                    expect=expected_sentinels[filename],
                )
            )

    if not skip_tests:
        stages.append(
            run_stage(
                "exact arithmetic unit tests",
                [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
                expect="OK",
            )
        )

    total_elapsed = time.perf_counter() - started
    try:
        maxrss_kb = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    except Exception:
        maxrss_kb = None

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed",
        "environment": env,
        "source_hashes": source_hashes,
        "manuscript_sha256": manuscript_hash,
        "author_results": {
            "n6": author_n6,
            "n7_n15": author_n7,
        },
        "ordinary_and_optimized_results_identical": optimized_match,
        "independent_checkers_run": not skip_independent,
        "unit_tests_run": not skip_tests,
        "stages": [concise_stage_record(stage) for stage in stages],
        "total_elapsed_seconds": round(total_elapsed, 6),
        "children_maxrss_kb": maxrss_kb,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_output(
        env=env,
        source_hashes=source_hashes,
        stages=stages,
        n6=author_n6,
        n7=author_n7,
        total_elapsed=total_elapsed,
        maxrss_kb=maxrss_kb,
    )
    write_sha256sums()
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-independent", action="store_true", help="skip the three independent checkers")
    parser.add_argument("--skip-optimized", action="store_true", help="skip duplicate -O author-verifier runs")
    parser.add_argument("--skip-tests", action="store_true", help="skip exact unit tests")
    parser.add_argument("--json", action="store_true", help="print the final manifest to stdout")
    parser.add_argument("--hashes-only", action="store_true", help="refresh SHA256SUMS without running certificates")
    args = parser.parse_args()

    if args.hashes_only:
        write_sha256sums()
        print(f"Wrote {HASHES_PATH.name}.")
        return

    manifest = verify_full(
        skip_independent=args.skip_independent,
        skip_optimized=args.skip_optimized,
        skip_tests=args.skip_tests,
    )
    if args.json:
        print(json.dumps(manifest, indent=2, sort_keys=True))
    print("\nAll Dittert certificates for dimensions 6 through 15 passed.")


if __name__ == "__main__":
    main()
