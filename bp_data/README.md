# DQI + Belief Propagation

**Reproduces:** DQI+BP column of Table 1 in [arXiv:2604.24633](https://arxiv.org/abs/2604.24633).

BP performance is computed analytically via density evolution — no random instances needed. The BP decoding threshold for each (k, D) is then converted to a DQI+BP satisfaction fraction using the semicircle law.

## Reproduce

```sh
# 1. Run density evolution (or use the pre-computed de_results.txt)
bash de.sh

# 2. Parse thresholds and compute DQI+BP satisfaction fractions
python3 generate_bp_table.py
```

## Pre-computed data

- **`de_results.txt`** — Density evolution output for all (k, D) with 3 ≤ k < D ≤ 8

## Note on build paths

The `de.sh` script references a Bazel build path (`./bazel-bin/src/density`). If you built with `make` from the repository root, either:
- Edit `de.sh` to use `../bin/density` instead, or
- Just use the pre-computed `de_results.txt`

## How it works

1. `./bazel-bin/src/density --k K --D D --bins 2000 --iterations 2000 --precision 0.000001` runs density evolution to find the BP decoding threshold for a (k, D)-regular LDPC code (use `../bin/density` for Make builds).
2. `generate_bp_table.py` parses the thresholds and applies the semicircle law: s(threshold) = 1/2 + sqrt(threshold * (1 - threshold)).
