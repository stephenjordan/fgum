#!/usr/bin/env python3
"""Regenerate or replay the fifteen FGUM interval certificates."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import mpmath

sys.dont_write_bytecode = True
from generate_tables import generate_tables


HERE = Path(__file__).resolve().parent
PAIRS = tuple((k, D) for k in range(3, 8) for D in range(k + 1, 9))


def run_pair(k, D, output, dps, verify_only):
    """Each CLI invocation gets its own process and interval context."""
    target = output / f"certificate_k{k}_D{D}.json"
    commands = []
    if not verify_only:
        commands.append([
            sys.executable, "-B", str(HERE / "certify_core.py"),
            "--k", str(k), "--D", str(D), "--dps", str(dps),
            "--output", str(target),
        ])
    commands.append([
        sys.executable, "-B", str(HERE / "certify_core.py"),
        "--verify", str(target),
    ])
    for command in commands:
        completed = subprocess.run(command, text=True, capture_output=True)
        if completed.returncode:
            detail = completed.stderr.strip() or completed.stdout.strip()
            raise RuntimeError(f"({k},{D}): {detail}")
        response = json.loads(completed.stdout)
    if response.get("verified") is not True:
        raise RuntimeError(f"({k},{D}): replay did not report verification")
    return {"k": k, "D": D, "file": target.name, "verified": True}


def write_manifest(output, status, dps, results, failures):
    sources = ("certify_core.py", "run_all.py", "generate_tables.py",
               "reference_values.json", "requirements.txt")
    manifest = {
        "schema": "fgum-certificate-run-v1",
        "status": status,
        "provenance": "Fresh reconstruction from the manuscript equations",
        "python": platform.python_version(),
        "mpmath": mpmath.__version__,
        "interval_decimal_precision": dps,
        "expected_parameter_pairs": [list(pair) for pair in PAIRS],
        "verified_parameter_pairs": sorted(results, key=lambda r: (r["k"], r["D"])),
        "failures": failures,
        "source_sha256": {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in sources
        },
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=HERE / "generated",
        help="Result directory (default: reproducibility/generated).",
    )
    parser.add_argument(
        "--jobs", type=int, default=min(4, os.cpu_count() or 1),
        help="Maximum concurrent parameter pairs (default: up to four).",
    )
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument(
        "--verify-only", action="store_true",
        help="Replay existing certificates and check their tables without writing.",
    )
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.dps < 70:
        parser.error("--dps must be at least 70")
    if mpmath.__version__ != "1.3.0":
        parser.error("Install the pinned dependency mpmath==1.3.0")
    output = args.output.resolve()
    if args.verify_only:
        if not output.is_dir():
            parser.error("--verify-only requires an existing result directory")
    else:
        output.mkdir(parents=True, exist_ok=True)
        write_manifest(output, "running", args.dps, [], [])
    results, failures = [], []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        jobs = {
            pool.submit(run_pair, k, D, output, args.dps, args.verify_only): (k, D)
            for k, D in PAIRS
        }
        for future in as_completed(jobs):
            k, D = jobs[future]
            try:
                results.append(future.result())
                print(f"PASS ({k},{D})", flush=True)
            except Exception as exc:
                failures.append({"k": k, "D": D, "error": str(exc)})
                print(f"FAIL ({k},{D}): {exc}", file=sys.stderr, flush=True)
    if not failures:
        try:
            generate_tables(output, check=args.verify_only)
        except Exception as exc:
            failures.append({"error": f"Table verification: {exc}"})
    if not args.verify_only:
        write_manifest(output, "failed" if failures else "passed",
                       args.dps, results, failures)
    if failures:
        print(json.dumps({"status": "failed", "failures": failures}, indent=2),
              file=sys.stderr)
        return 1
    print(f"Verified all {len(PAIRS)} certificates and manuscript table values.")
    print(f"Results: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
