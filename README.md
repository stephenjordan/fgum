# Optimization Using Locally-Quantum Decoders

Companion code and data for:

> N. Shutty, S. P. Jordan, M. Shao, S. Hadfield, J. S. Azariah, G. Smith,
> "Optimization Using Locally-Quantum Decoders,"
> [arXiv:2604.24633](https://arxiv.org/abs/2604.24633) (2026).

This repository contains everything needed to reproduce the numerical results in Table 1 of the paper. Each algorithmic approach — simulated annealing, belief propagation via density evolution, Regev+FGUM, and QAOA — has its own data directory and scripts. A self-contained [QAOA verifier](#5-qaoa) is included so that any researcher can independently confirm the QAOA satisfaction fractions from the optimised angle tables.

---

## Table of Contents

| Section | What it reproduces | Language | Runtime |
|---|---|---|---|
| [1. Build](#1-build) | C++ executables for SA and density evolution | C++20 | < 1 min |
| [2. Simulated Annealing](#2-simulated-annealing) | SA rows of Table 1 | C++ + Python | Hours (cluster) |
| [3. DQI + Belief Propagation](#3-dqi--belief-propagation) | DQI+BP rows of Table 1 | C++ + Python | Minutes |
| [4. Regev + FGUM](#4-regev--fgum) | FGUM rows of Table 1 | Python | Seconds |
| [5. QAOA](#5-qaoa) | QAOA rows of Table 1 | Julia | Seconds–hours |

---

## Repository Structure

```
fgum/
├── README.md
├── Makefile                          # Builds all C++ executables
├── src/                              # C++ source code
│   ├── full_anneal.cpp               #   Simulated annealing solver
│   ├── random_regular.cpp            #   Random regular XORSAT instance generator
│   ├── density.cpp                   #   Density evolution for BP threshold
│   ├── density_evolution.{cpp,hpp}   #   Density evolution core
│   ├── distribution.{cpp,hpp}        #   Degree distribution utilities
│   ├── bitmatrix.{cpp,hpp}           #   GF(2) linear algebra
│   ├── xoropt.{cpp,hpp}              #   XOR-SAT optimisation routines
│   ├── utils.{cpp,hpp}               #   Shared utilities
│   ├── test_xortools.cpp             #   Unit tests for bitmatrix/xortools
│   └── argparse.hpp                  #   Argument parser (vendored, MIT)
│
├── sa_data/                          # Simulated annealing data
│   ├── generate_sa_jobs.py           #   Generates instance + anneal commands
│   ├── generate_sa_table.py          #   Parses results → Table 1 SA rows
│   ├── sa_jobs.txt                   #   Generated command list
│   ├── summary.txt                   #   Cluster output (grep satisfied log*)
│   └── instance_*.tsv                #   Random D-regular k-XORSAT instances
│
├── bp_data/                          # Belief propagation via density evolution
│   ├── de.sh                         #   Runs density evolution for all (k,D)
│   ├── de_results.txt                #   Density evolution output
│   └── generate_bp_table.py          #   Parses thresholds → Table 1 DQI+BP rows
│
├── fgum_data/                        # Regev + FGUM calculations
│   └── generate_fgum_table.py        #   Computes FGUM performance → Table 1 rows
│
├── qaoa_data/                        # QAOA optimised angles
│   └── qaoa_angles_for_verifier.csv  #   254 optimised (k,D,p) configurations
│
├── qaoa_verifier/                    # Standalone QAOA verification tool
│   ├── Project.toml                  #   Julia package metadata (zero dependencies)
│   ├── verify.jl                     #   CLI entry point
│   ├── src/QaoaVerifier.jl           #   Verification module (~400 lines)
│   └── test/runtests.jl              #   42 unit tests
│
└── test_data/                        # Bitmatrix/xortools test fixtures
    └── *.txt
```

---

## 1. Build

**Requirements:** A C++20 compiler (GCC 10+, Clang 13+, or Apple Clang 14+) and GNU Make.

```sh
make
```

This builds four executables in `bin/`:

| Executable | Purpose |
|---|---|
| `bin/full_anneal` | Simulated annealing for Max-k-XORSAT instances |
| `bin/random_regular` | Generate random D-regular k-XORSAT instances from Gallager's ensemble |
| `bin/density` | Density evolution for computing BP decoding thresholds |
| `bin/test_xortools` | Unit tests for the GF(2) linear algebra library |

To run the bundled test suite (bitmatrix unit tests + density evolution smoke test):

```sh
make check
```

To clean up:

```sh
make clean
```

---

## 2. Simulated Annealing

**Corresponds to:** SA column of Table 1.

Simulated annealing (SA) performance is determined empirically: generate random D-regular Max-k-XORSAT instances from Gallager's ensemble (with uniformly random variable signs), then run SA and record the best satisfaction fraction across multiple trials.

### 2.1 Generate problem instances and job commands

```sh
cd sa_data
python3 generate_sa_jobs.py
```

This produces `sa_jobs.txt`, which contains commands to:
1. Generate random instances via `random_regular k D n` (stored as `instance_k_D_n.tsv`)
2. Run `full_anneal` on each instance for 1,000,000 sweeps, 4 independent trials per instance

### 2.2 Run the anneals

The commands in `sa_jobs.txt` can be run directly, but the full set takes substantial time. We ran them in parallel on a cluster. Pre-computed results are included:

- **`instance_*.tsv`** — 15 random instances covering all (k,D) pairs with 3 ≤ k < D ≤ 8
- **`summary.txt`** — Aggregated output (`grep satisfied log*`)

### 2.3 Tabulate results

```sh
cd sa_data
python3 generate_sa_table.py
```

For each (k,D), this takes the best of 4 trials and reports the satisfaction ratio.

### Instance format

Each `.tsv` file encodes a random D-regular k-XORSAT instance:
- Line 1: `n m` (number of variables, number of clauses)
- Lines 2+: `sign var₁ var₂ … varₖ` (clause sign ±0.5 and 0-indexed variable indices)

---

## 3. DQI + Belief Propagation

**Corresponds to:** DQI+BP column of Table 1.

The BP performance is computed analytically via density evolution, then converted to DQI+BP performance using the semicircle law.

### 3.1 Run density evolution

```sh
cd bp_data
bash de.sh
```

This runs `bin/density` for all (k,D) pairs with 3 ≤ k < D ≤ 8 at high precision (2000 bins, 2000 iterations). Each run finds the BP decoding threshold — the noise level below which BP decoding succeeds with high probability. Pre-computed results are in `de_results.txt`.

> **Note:** The `de.sh` script references a Bazel build path. To run with the Makefile-built binary, replace `./bazel-bin/src/density` with `../bin/density` in `de.sh`, or simply use the pre-computed results.

### 3.2 Compute DQI+BP performance

```sh
cd bp_data
python3 generate_bp_table.py
```

This parses the density evolution thresholds and applies the semicircle law to obtain the DQI+BP satisfaction fraction for each (k,D).

---

## 4. Regev + FGUM

**Corresponds to:** FGUM column of Table 1.

The FGUM performance is computed entirely from the closed-form equations in the manuscript, with one numerical rootfinding step for e_max.

```sh
cd fgum_data
python3 generate_fgum_table.py
```

This solves for e_max(k,D) by bisection, computes α_min and the expected weight of a uniformly random element of I*_D, and produces the Regev+FGUM satisfaction fraction for all (k,D) with 3 ≤ k < D ≤ 8.

**No build step required** — pure Python with only standard-library dependencies.

---

## 5. QAOA

**Corresponds to:** QAOA column of Table 1.

The QAOA satisfaction fractions are the most computationally demanding results in the paper. For each constraint arity k and variable degree D, the depth-p QAOA circuit's expected performance on D-regular Max-k-XORSAT is computed exactly via a tensor network contraction on the light-cone tree. The angles (γ, β) were optimised using gradient-based methods on GPUs.

### 5.1 Optimised angles

The file `qaoa_data/qaoa_angles_for_verifier.csv` contains 254 optimised configurations covering:

| Parameters | Depth range | Count |
|---|---|---|
| MaxCut: (k=2, D=3) | p = 1 … 16 | 16 |
| XORSAT: all (k,D) with 3 ≤ k < D ≤ 8 | p = 1 … 16 | 238 |

CSV format:
```
k,D,p,ctilde,source,gamma,beta
```
where `gamma` and `beta` are semicolon-separated angle vectors of length p.

The angles were computed using two independent implementations:
- [QOKit](https://github.com/jpmorganchase/QOKit/tree/add-max-k-xor-sat/qokit/max_k_xor_sat) — C++, double-double precision (JP Morgan Chase)
- [QaoaXorsat.jl](https://github.com/johnazariah/qaoa-xorsat) — Julia, Float64 + Double64 (Centre for Quantum Software and Information, UTS)

### 5.2 Verification

A self-contained Julia verifier (`qaoa_verifier/`) independently recomputes every satisfaction fraction from the optimised angles. It implements the Basso/Villalonga branch tensor recurrence with Walsh-Hadamard-accelerated constraint folding and produces bit-exact agreement (Δ < 10⁻¹²) with the claimed values.

**Prerequisites:** [Julia 1.9+](https://julialang.org/downloads/) (no other dependencies).

**Verify all 254 configurations:**

```sh
cd qaoa_verifier
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. verify.jl
```

Expected output (truncated):
```
╔══════════════════════════════════════════════════════════════════╗
║  QAOA Satisfaction Fraction Verifier                            ║
║  Companion to arXiv:2604.24633                                 ║
╚══════════════════════════════════════════════════════════════════╝

Found 254 configurations to verify.
────────────────────────────────────────────────────────────────────
  (k=2, D=3, p= 1) [MaxCut]  c̃=0.692450089730  OK (Δ=1.25e-13)
  (k=2, D=3, p= 2) [MaxCut]  c̃=0.755906458214  OK (Δ=4.12e-13)
  ...
  (k=3, D=4, p= 8) [XORSAT]  c̃=0.854102048055  OK (Δ=3.91e-13)
  ...
✓  All 254 evaluations verified successfully.
```

Low-depth configurations (p ≤ 8) verify in under a second. Higher depths scale as O(p² · 4^p): p = 12 takes ~1 minute, p = 16 takes ~30 minutes.

**Run the unit test suite:**

```sh
cd qaoa_verifier
julia --project=. -e 'using Pkg; Pkg.test()'
```

**Verify a single configuration:**

```sh
cd qaoa_verifier
julia --project=. verify.jl --k 3 --D 4 \
    --gamma "0.550050516856796" \
    --beta "0.287371347852009" \
    --ctilde 0.676056660061
```

**Verify a custom CSV:**

```sh
cd qaoa_verifier
julia --project=. verify.jl path/to/angles.csv
```

The verifier source is a single ~400-line module at `qaoa_verifier/src/QaoaVerifier.jl` with no external dependencies. See [qaoa_verifier/README.md](qaoa_verifier/README.md) for the mathematical details of the evaluation pipeline.

---

## References

- N. Shutty, S. P. Jordan, M. Shao, S. Hadfield, J. S. Azariah, G. Smith, "Optimization Using Locally-Quantum Decoders," [arXiv:2604.24633](https://arxiv.org/abs/2604.24633) (2026).
- M. Basso, E. Farhi, K. Marwaha, B. Villalonga, L. Zhou, "The Quantum Approximate Optimization Algorithm at High Depth for MaxCut on Large-Girth Regular Graphs and the Sherrington-Kirkpatrick Model," [arXiv:2110.14206](https://arxiv.org/abs/2110.14206) (2022).
- E. Farhi, J. Goldstone, S. Gutmann, "A Quantum Approximate Optimization Algorithm," [arXiv:1411.4028](https://arxiv.org/abs/1411.4028) (2014).

## License

MIT. See individual source files for third-party licenses (argparse.hpp).
