"""Fast, independent checks of the optional exact-rank experiments.

The reference rank routine uses a dense matrix and row elimination.  The
reference peeling routine recomputes every row degree after each deletion.
Neither uses the integer-column basis or incremental bookkeeping under test.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
from itertools import product
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest


MODULE_PATH = Path(__file__).with_name("exact_rank_experiments.py")


def dense_rank(columns: list[list[int]], nrows: int) -> int:
    """Independent GF(2) reference using dense row-echelon elimination."""
    rows = [[0 for _ in columns] for _ in range(nrows)]
    for j, support in enumerate(columns):
        for i in support:
            rows[i][j] = 1
    rank = 0
    for j in range(len(columns)):
        pivot = next((i for i in range(rank, nrows) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(rank + 1, nrows):
            if rows[i][j]:
                rows[i] = [a ^ b for a, b in zip(rows[i], rows[rank])]
        rank += 1
        if rank == nrows:
            break
    return rank


def reference_core(
    columns: list[list[int]], nrows: int, reverse: bool = False
) -> tuple[list[int], int]:
    """Strip one weight-one row at a time, rebuilding degrees from scratch."""
    alive = set(range(len(columns)))
    row_order = range(nrows - 1, -1, -1) if reverse else range(nrows)
    while True:
        incident = [[j for j in alive if i in columns[j]] for i in range(nrows)]
        pivot = next((i for i in row_order if len(incident[i]) == 1), None)
        if pivot is None:
            return sorted(alive), sum(bool(row) for row in incident)
        alive.remove(incident[pivot][0])


def load_experiments():
    spec = importlib.util.spec_from_file_location("fgum_exact_rank_under_test", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {MODULE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class ScriptedBlockRandom:
    """Make the channel's Bernoulli decisions observable without sampling noise."""

    def __init__(self, draws=()):
        self.draws = iter(draws)
        self.random_calls = 0
        self.sample_calls = []

    def random(self):
        self.random_calls += 1
        return next(self.draws)

    def sample(self, population, count):
        population = list(population)
        self.sample_calls.append((population, count))
        return population[::-1][:count]


class ExactRankExperimentsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.experiments = load_experiments()

    def assert_core_matches_reference(self, columns, nrows):
        expected, expected_rows = reference_core(columns, nrows)
        reverse_expected, _ = reference_core(columns, nrows, reverse=True)
        self.assertEqual(expected, reverse_expected)
        core = self.experiments.peel_core(columns, nrows)
        self.assertEqual(core.surviving_column_indices, tuple(expected))
        remaining = [columns[j] for j in expected]
        row_support = sorted({i for column in remaining for i in column})
        self.assertEqual(core.surviving_row_indices, tuple(row_support))
        self.assertEqual(core.core_rows, expected_rows)
        self.assertEqual(core.core_columns, len(expected))
        self.assertEqual(core.peeled_columns, len(columns) - len(expected))
        original_rank = dense_rank(columns, nrows)
        core_rank = dense_rank(remaining, nrows)
        self.assertEqual(original_rank - core_rank, core.peeled_columns)
        self.assertEqual(len(columns) - original_rank, len(remaining) - core_rank)
        for row in row_support:
            self.assertGreaterEqual(sum(row in column for column in remaining), 2)

    def test_rank_against_all_binary_three_by_three_matrices(self):
        for entries in product((0, 1), repeat=9):
            columns = [[i for i in range(3) if entries[3 * i + j]] for j in range(3)]
            with self.subTest(entries=entries):
                self.assertEqual(
                    self.experiments.gf2_column_rank(columns, 3), dense_rank(columns, 3)
                )

    def test_rank_empty_rectangular_and_dependent_cases(self):
        fixtures = [
            ([], 0, 0),
            ([[], []], 3, 0),
            ([[0], [0], [1], [0, 1]], 4, 2),
            ([[0, 1], [1, 2], [0, 2]], 3, 2),
            ([[0, 1, 3], [0, 2, 3], [1, 2, 3]], 4, 3),
        ]
        for columns, nrows, expected in fixtures:
            with self.subTest(columns=columns, nrows=nrows):
                self.assertEqual(dense_rank(columns, nrows), expected)
                self.assertEqual(self.experiments.gf2_column_rank(columns, nrows), expected)
        rng = random.Random(619)
        for _ in range(40):
            nrows, ncols = rng.randrange(1, 9), rng.randrange(11)
            columns = [[i for i in range(nrows) if rng.randrange(2)] for _ in range(ncols)]
            self.assertEqual(
                self.experiments.gf2_column_rank(columns, nrows), dense_rank(columns, nrows)
            )

    def test_peeling_cascade_and_order_independence(self):
        # Two initially peelable chains meet a dependent two-column core.
        columns = [[0, 1], [1, 5], [2, 3], [3, 5], [5, 6], [5, 6]]
        self.assertEqual(reference_core(columns, 7)[0], [4, 5])
        self.assert_core_matches_reference(columns, 7)
        permutation = [5, 2, 0, 4, 1, 3]
        relabeled = [[6 - row for row in columns[j]] for j in permutation]
        core = self.experiments.peel_core(relabeled, 7)
        self.assertEqual(sorted(permutation[j] for j in core.surviving_column_indices), [4, 5])
        self.assert_core_matches_reference(relabeled, 7)

    def test_zero_columns_survive_and_nonempty_core_can_have_full_rank(self):
        self.assert_core_matches_reference([[], [0, 1], [1, 2]], 3)
        self.assertEqual(
            self.experiments.peel_core([[], [0, 1], [1, 2]], 3).surviving_column_indices,
            (0,),
        )
        columns = [[0, 1, 3], [0, 2, 3], [1, 2, 3]]
        core = self.experiments.peel_core(columns, 4)
        self.assertEqual(core.core_columns, 3)
        self.assertEqual(self.experiments.gf2_column_rank(columns, 4), 3)

    def test_peeling_preserves_nullity_on_small_varied_matrices(self):
        rng = random.Random(1601)
        for case in range(40):
            nrows, ncols = rng.randrange(1, 8), rng.randrange(10)
            columns = [[i for i in range(nrows) if rng.random() < 0.3] for _ in range(ncols)]
            with self.subTest(case=case):
                self.assert_core_matches_reference(columns, nrows)

    def test_gallager_sampler_layer_degrees_and_first_layer_blocks(self):
        for k, D, blocks in [(3, 4, 4), (4, 5, 3), (7, 8, 2)]:
            with self.subTest(k=k, D=D, blocks=blocks):
                columns = self.experiments.sample_gallager_columns(
                    k, D, blocks, random.Random(402)
                )
                self.assertEqual(len(columns), D * blocks)
                degrees = [0] * (k * blocks)
                for j, support in enumerate(columns):
                    self.assertEqual(len(support), k)
                    self.assertEqual(support[0], j // D)
                    self.assertEqual(sorted(row // blocks for row in support), list(range(k)))
                    for row in support:
                        self.assertGreaterEqual(row, 0)
                        self.assertLess(row, k * blocks)
                        degrees[row] += 1
                self.assertEqual(degrees, [D] * (k * blocks))
                self.assertEqual(
                    columns,
                    self.experiments.sample_gallager_columns(k, D, blocks, random.Random(402)),
                )

    def test_channel_models_have_distinct_documented_semantics(self):
        rng = ScriptedBlockRandom([0.1, 0.39, 0.4, 0.8, 0.2])
        selected = self.experiments.select_failed_blocks(5, 0.4, rng, "bernoulli")
        self.assertEqual(selected, (0, 1, 4))
        self.assertEqual(rng.random_calls, 5)
        self.assertEqual(rng.sample_calls, [])
        fixed_rng = ScriptedBlockRandom()
        selected = self.experiments.select_failed_blocks(5, 0.5, fixed_rng, "fixed-count")
        self.assertEqual(selected, (2, 3, 4))  # Documented half-up rounding: 2.5 -> 3.
        self.assertEqual(fixed_rng.sample_calls, [(list(range(5)), 3)])
        self.assertEqual(fixed_rng.random_calls, 0)

    def test_no_failures_and_all_failures(self):
        for model in ("bernoulli", "fixed-count"):
            with self.subTest(model=model):
                empty = self.experiments.run_trial(3, 4, 1, 0.0, 71, model)
                self.assertEqual(empty.failed_blocks, 0)
                self.assertEqual(empty.actual_failed_fraction, 0)
                self.assertEqual((empty.columns, empty.rank, empty.nullity), (0, 0, 0))
                self.assertTrue(empty.full_column_rank)
                self.assertEqual((empty.core_columns, empty.core_rows, empty.core_nullity), (0, 0, 0))
                full = self.experiments.run_trial(3, 4, 1, 1.0, 71, model)
                # With one row in each layer all four columns are identical.
                self.assertEqual(full.failed_blocks, 1)
                self.assertEqual(full.actual_failed_fraction, 1)
                self.assertEqual((full.columns, full.rank, full.nullity), (4, 1, 3))
                self.assertFalse(full.full_column_rank)
                self.assertEqual((full.core_columns, full.core_rows, full.core_rank), (4, 3, 1))
                self.assertEqual(full.core_nullity, 3)
                self.assertEqual(full.peeled_columns, 0)

    def test_trial_reconstructs_from_recorded_random_streams(self):
        for model in ("bernoulli", "fixed-count"):
            for trial_seed in (5, 191, 6409):
                with self.subTest(model=model, seed=trial_seed):
                    result = self.experiments.run_trial(3, 4, 4, 0.6, trial_seed, model)
                    self.assertNotEqual(result.matrix_seed, result.failure_seed)
                    columns = self.experiments.sample_gallager_columns(
                        3, 4, 4, random.Random(result.matrix_seed)
                    )
                    failed = self.experiments.select_failed_blocks(
                        4, 0.6, random.Random(result.failure_seed), model
                    )
                    restricted = [column for j, column in enumerate(columns) if j // 4 in failed]
                    rank = dense_rank(restricted, 12)
                    surviving, n_core_rows = reference_core(restricted, 12)
                    residual = [restricted[j] for j in surviving]
                    core_rank = dense_rank(residual, 12)
                    self.assertEqual(result.failed_blocks, len(failed))
                    self.assertEqual(result.actual_failed_fraction, len(failed) / 4)
                    self.assertEqual(result.columns, 4 * len(failed))
                    self.assertEqual(result.nonzero_rows, len({i for col in restricted for i in col}))
                    self.assertEqual(result.rank, rank)
                    self.assertEqual(result.nullity, len(restricted) - rank)
                    self.assertEqual(result.full_column_rank, rank == len(restricted))
                    self.assertEqual((result.core_columns, result.core_rows), (len(residual), n_core_rows))
                    self.assertEqual(result.core_rank, core_rank)
                    self.assertEqual(result.core_nullity, len(residual) - core_rank)
                    self.assertEqual(result.nullity, result.core_nullity)
                    self.assertEqual(result.peeled_columns, len(restricted) - len(residual))

    def test_per_trial_seeds_and_grid_order_independence(self):
        seed = self.experiments.derive_trial_seed(17, 3, 4, 4, 0.5, 0, "bernoulli")
        self.assertIsInstance(seed, int)
        self.assertGreaterEqual(seed, 0)
        self.assertLess(seed, 1 << 256)
        self.assertEqual(seed, self.experiments.derive_trial_seed(17, 3, 4, 4, 0.5, 0, "bernoulli"))
        alternatives = [
            (18, 3, 4, 4, 0.5, 0, "bernoulli"),
            (17, 4, 5, 4, 0.5, 0, "bernoulli"),
            (17, 3, 4, 5, 0.5, 0, "bernoulli"),
            (17, 3, 4, 4, 0.6, 0, "bernoulli"),
            (17, 3, 4, 4, 0.5, 1, "bernoulli"),
            (17, 3, 4, 4, 0.5, 0, "fixed-count"),
        ]
        seeds = [seed] + [self.experiments.derive_trial_seed(*args) for args in alternatives]
        self.assertEqual(len(set(seeds)), len(seeds))
        forward = self.experiments.run_grid(3, 4, [2, 4], [0.2, 0.6], 4, 17, "bernoulli")
        backward = self.experiments.run_grid(3, 4, [4, 2], [0.6, 0.2], 4, 17, "bernoulli")
        expanded = self.experiments.run_grid(
            3, 4, [4, 3, 2, 4], [0.6, 0.9, 0.2, 0.6], 4, 17, "bernoulli"
        )
        self.assertEqual(forward, backward)
        self.assertEqual(
            forward,
            [row for row in expanded if row["rate"] in (0.2, 0.6) and row["blocks"] in (2, 4)],
        )
        longer = self.experiments.run_grid(3, 4, [2, 4], [0.2, 0.6], 6, 17, "bernoulli")
        self.assertEqual(forward, [row for row in longer if row["trial_index"] < 4])
        self.assertEqual(len({row["trial_seed"] for row in forward}), len(forward))

    def test_summary_is_computed_from_raw_trials(self):
        raw = []
        for model in ("bernoulli", "fixed-count"):
            raw.extend(self.experiments.run_grid(3, 4, [2, 3], [0.0, 0.5, 1.0], 4, 991, model))
        summaries = self.experiments.summarize_trials(raw)
        self.assertEqual(len(summaries), 12)
        mean_fields = (
            "failed_blocks", "actual_failed_fraction", "columns", "nonzero_rows",
            "rank", "nullity", "core_columns", "core_rows", "core_rank",
            "core_nullity", "peeled_columns",
        )
        for summary in summaries:
            sample = [row for row in raw if row["rate_hex"] == summary["rate_hex"]
                      and row["erasure_model"] == summary["erasure_model"]
                      and row["blocks"] == summary["blocks"]]
            self.assertEqual(summary["trials"], len(sample))
            successes = sum(row["full_column_rank"] for row in sample)
            self.assertEqual(summary["full_column_rank_trials"], successes)
            self.assertEqual(summary["full_column_rank_fraction"], successes / len(sample))
            for field in mean_fields:
                expected = math.fsum(row[field] for row in sample) / len(sample)
                self.assertAlmostEqual(summary["mean_" + field], expected, places=14)

    def cli(self, *args):
        # -S disables site packages: the experiment runner must use only stdlib.
        return subprocess.run(
            [sys.executable, "-S", str(MODULE_PATH), *map(str, args)],
            capture_output=True, text=True, timeout=20,
        )

    def test_cli_writes_raw_summary_and_manifest_without_overwriting(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp) / "run"
            args = ("--k", 3, "--D", 4, "--blocks", 2, "--rates", 0, 0.5, 1,
                    "--trials", 2, "--seed", 17, "--outdir", out)
            completed = self.cli(*args)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            paths = [out / name for name in ("trials.csv", "summary.csv", "run.json")]
            self.assertTrue(all(path.is_file() for path in paths))
            with paths[0].open(newline="") as handle:
                raw = list(csv.DictReader(handle))
            with paths[1].open(newline="") as handle:
                summaries = list(csv.DictReader(handle))
            self.assertEqual(len(raw), 6)
            self.assertEqual(len(summaries), 3)
            expected = self.experiments.run_grid(3, 4, [2], [0.0, 0.5, 1.0], 2, 17, "bernoulli")
            for written, trial in zip(raw, expected):
                self.assertEqual(int(written["trial_seed"]), trial["trial_seed"])
                self.assertEqual(int(written["rank"]), trial["rank"])
                self.assertEqual(int(written["nullity"]), trial["nullity"])
            for summary in summaries:
                sample = [row for row in raw if row["rate_hex"] == summary["rate_hex"]]
                self.assertEqual(int(summary["trials"]), len(sample))
                self.assertAlmostEqual(
                    float(summary["mean_nullity"]),
                    sum(int(row["nullity"]) for row in sample) / len(sample),
                )
            manifest = json.loads(paths[2].read_text())
            self.assertEqual(manifest["parameters"]["total_trials"], 6)
            self.assertEqual(manifest["parameters"]["base_seed"], 17)
            self.assertEqual(manifest["parameters"]["erasure_model"], "bernoulli")
            for path, row_count in ((paths[0], 6), (paths[1], 3)):
                recorded = manifest["output_files"][path.name]
                self.assertEqual(recorded["rows"], row_count)
                self.assertEqual(recorded["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertEqual(
                manifest["script"]["sha256"], hashlib.sha256(MODULE_PATH.read_bytes()).hexdigest()
            )
            self.assertEqual((manifest["reference_certificate"]["k"],
                              manifest["reference_certificate"]["D"]), (3, 4))
            original = {path.name: path.read_bytes() for path in paths}
            repeated = self.cli(*args)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertEqual(original, {path.name: path.read_bytes() for path in paths})

    def test_cli_rejects_invalid_parameters(self):
        invalid = [
            ("--k", "0"), ("--D", "0"), ("--blocks", "0"),
            ("--trials", "0"), ("--rates", "-0.01"),
            ("--rates", "1.01"), ("--rates", "nan"),
            ("--rates", "inf"), ("--erasure-model", "unknown"),
        ]
        with tempfile.TemporaryDirectory() as temp:
            for index, override in enumerate(invalid):
                with self.subTest(arguments=override):
                    completed = self.cli(
                        "--k", 3, "--D", 4, "--blocks", 2, "--rates", 0.5,
                        "--trials", 1, "--outdir", Path(temp) / str(index), *override,
                    )
                    self.assertNotEqual(completed.returncode, 0)
                    self.assertTrue(completed.stderr.strip())
                    self.assertFalse((Path(temp) / str(index)).exists())


if __name__ == "__main__":
    unittest.main()
