#!/usr/bin/env python3
"""Make CSV tables from interval certificates and check manuscript rounding.

This script formats already computed results. Use certify_core.py or
run_all.py to verify the underlying mathematical inequalities.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
import io
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def check_rounding(interval, expected, name):
    expected_value = Decimal(expected)
    quantum = Decimal(1).scaleb(expected_value.as_tuple().exponent)
    with localcontext() as ctx:
        ctx.prec = 90
        for side in ("lower", "upper"):
            got = Decimal(interval[side]).quantize(quantum, rounding=ROUND_HALF_EVEN)
            if got != expected_value:
                raise ValueError(f"{name}: {side} rounds to {got}, expected {expected}")


def csv_text(rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def generate_tables(output, check=False):
    references = json.loads((HERE / "reference_values.json").read_text(
        encoding="utf-8"))["pairs"]
    expected_pairs = {(k, D) for k in range(3, 8) for D in range(k + 1, 9)}
    found_pairs = [(row["k"], row["D"]) for row in references]
    if len(found_pairs) != 15 or set(found_pairs) != expected_pairs:
        raise ValueError("Reference data do not contain exactly the fifteen pairs")
    certificate_rows, score_rows = [], []
    for reference in sorted(references, key=lambda r: (r["k"], r["D"])):
        k, D = reference["k"], reference["D"]
        doc = json.loads((output / f"certificate_k{k}_D{D}.json").read_text(
            encoding="utf-8"))
        if (doc["k"], doc["D"]) != (k, D):
            raise ValueError("Certificate parameter pair does not match its filename")
        summary = doc["summary"]
        check_rounding(summary["e_core"], reference["e_core"], f"({k},{D}) e_core")
        check_rounding(summary["fgum_score"], reference["fgum_score"],
                       f"({k},{D}) FGUM score")
        check_rounding(summary["beta"], reference["beta_printed"],
                       f"({k},{D}) beta")
        targets = reference["upper_bounds"]
        for field, actual in (
            ("exterior_U", summary["exterior_U_upper"]),
            ("central_Gamma_second", summary["central_Gamma_second_upper"]),
            ("spectral_radius", summary["spectral_radius"]["upper"]),
        ):
            if Decimal(actual) > Decimal(targets[field]):
                raise ValueError(f"({k},{D}) does not verify manuscript bound {field}")
        certificate_rows.append({
            "k": k, "D": D,
            "e_core_lower": summary["e_core"]["lower"],
            "e_core_upper": summary["e_core"]["upper"],
            "e_core_printed": reference["e_core"],
            "beta_upper": summary["beta"]["upper"],
            "beta_printed": reference["beta_printed"],
            "exterior_U_upper": summary["exterior_U_upper"],
            "central_Gamma_second_upper": summary["central_Gamma_second_upper"],
            "spectral_radius_upper": summary["spectral_radius"]["upper"],
            "manuscript_exterior_U_bound": targets["exterior_U"],
            "manuscript_central_Gamma_second_bound": targets["central_Gamma_second"],
            "manuscript_spectral_radius_bound": targets["spectral_radius"],
            "alpha_max_upper": summary["alpha_max"]["upper"],
            "high_branch_T_lower": summary["branch_T_lower"],
            "branch_cells": summary["cells"]["branch"],
            "central_cells": summary["cells"]["central"],
            "exterior_cells": summary["cells"]["exterior"],
        })
        score_rows.append({
            "k": k, "D": D,
            "e_core_lower": summary["e_core"]["lower"],
            "e_core_upper": summary["e_core"]["upper"],
            "ihat_star_exact": summary["ihat_star_exact"],
            "alpha_min_lower": summary["alpha_min"]["lower"],
            "alpha_min_upper": summary["alpha_min"]["upper"],
            "fgum_score_lower": summary["fgum_score"]["lower"],
            "fgum_score_upper": summary["fgum_score"]["upper"],
            "fgum_score_printed": reference["fgum_score"],
        })
    contents = {
        "certificate_summary.csv": csv_text(certificate_rows),
        "fgum_scores.csv": csv_text(score_rows),
    }
    for name, text in contents.items():
        path = output / name
        if check:
            if path.read_text(encoding="utf-8") != text:
                raise ValueError(f"{name} differs from the verified certificate summaries")
        else:
            path.write_text(text, encoding="utf-8")
    return len(references)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=HERE / "results")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    count = generate_tables(args.input.resolve(), check=args.check)
    print(f"{'Checked' if args.check else 'Generated'} tables for {count} pairs.")


if __name__ == "__main__":
    main()
