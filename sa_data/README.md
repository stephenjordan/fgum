# Simulated Annealing

**Reproduces:** SA column of Table 1 in [arXiv:2604.24633](https://arxiv.org/abs/2604.24633).

SA performance is determined empirically: generate random D-regular Max-k-XORSAT instances from Gallager's ensemble, run simulated annealing, and take the best satisfaction fraction across 4 independent trials.

## Reproduce

```sh
# 1. Generate instance files and job commands
python3 generate_sa_jobs.py

# 2. Run the anneals (or use the pre-computed summary.txt)
#    The commands are in sa_jobs.txt — we ran them on a cluster.

# 3. Tabulate results
python3 generate_sa_table.py
```

## Pre-computed data

All results are included so step 2 can be skipped:

- **`instance_k_D_n.tsv`** — 15 random instances, one per (k, D) pair with 3 ≤ k < D ≤ 8
- **`summary.txt`** — Aggregated anneal output (`grep satisfied log*`)
- **`sa_jobs.txt`** — Full command list for reference

## Instance format

Each `.tsv` file encodes one D-regular k-XORSAT instance:

```
n  m                          ← number of variables, number of clauses
sign  var₁  var₂  …  varₖ    ← one line per clause (sign ±0.5, 0-indexed variables)
```

## Executables used

- `bin/random_regular k D n` — generates a random instance
- `bin/full_anneal instance.tsv sweeps` — runs SA with the given number of sweeps

Build with `make` from the repository root.
