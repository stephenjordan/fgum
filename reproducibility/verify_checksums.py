#!/usr/bin/env python3
"""Check the frozen FGUM package; --write deliberately refreshes its manifest."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MANIFEST = HERE / "SHA256SUMS"


def package_files():
    files = [ROOT / "README.md",
             ROOT / "fgum_data/README.md",
             ROOT / "fgum_data/generate_fgum_table.py"]
    files.extend(path for path in HERE.iterdir()
                 if path.is_file() and path != MANIFEST)
    results = HERE / "results"
    if not results.is_dir():
        raise ValueError("Missing frozen results directory")
    files.extend(path for path in results.rglob("*") if path.is_file())
    return sorted(files, key=lambda path: path.relative_to(ROOT).as_posix())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="Refresh checksums after an intentional package update.")
    args = parser.parse_args()
    files = package_files()
    current = "".join(f"{digest(path)}  {path.relative_to(ROOT).as_posix()}\n"
                      for path in files)
    if args.write:
        MANIFEST.write_text(current, encoding="utf-8")
        print(f"Wrote checksums for {len(files)} files.")
        return 0
    recorded = MANIFEST.read_text(encoding="utf-8")
    if recorded != current:
        print("CHECKSUM FAILURE: a frozen file changed, is missing, or was added.")
        old = dict(line.split("  ", 1)[::-1] for line in recorded.splitlines())
        new = dict(line.split("  ", 1)[::-1] for line in current.splitlines())
        for name in sorted(old.keys() | new.keys()):
            if old.get(name) != new.get(name):
                print(f"  {name}")
        return 1
    print(f"Verified SHA-256 checksums for {len(files)} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
