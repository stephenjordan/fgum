# Optimization Using Locally-Quantum Decoders

Companion code and data for:

> N. Shutty, A. Mandal, S. Ragavan, Q. Buzet, A. Chailloux, N. C. Rubin,
> A. Khan, S. Boulebnane, R. Shaydulin, J. Azariah, S. P. Jordan,
> "Optimization Using Locally-Quantum Decoders,"
> [arXiv:2604.24633](https://arxiv.org/abs/2604.24633) (2026).

This repository contains the scripts and data for Table 1, together with the interval checks supporting the FGUM thresholds in Section 5 and Appendix A. Each calculation has a directory with its own instructions.

Archived releases are available on Zenodo at [doi:10.5281/zenodo.20173447](https://doi.org/10.5281/zenodo.20173447). This DOI covers all versions and resolves to the latest published release.

## Quick Start

```sh
make          # build C++ tools (SA solver, instance generator, density evolution)
make check    # run unit tests
```

Requires a C++20 compiler and GNU Make. The table-generating Python scripts use only the standard library. The FGUM interval certificate additionally requires `mpmath==1.3.0`, as specified in [`reproducibility/requirements.txt`](reproducibility/requirements.txt). The Julia QAOA verifier requires [Julia 1.9+](https://julialang.org/downloads/). The C++ QAOA verifier requires CMake 3.16+ and a C++17 compiler.

## Table 1 — Directory Guide

| Directory | Table 1 Column | What to run | Details |
|---|---|---|---|
| [`sa_data/`](sa_data/) | Simulated Annealing | `python3 generate_sa_table.py` | [sa_data/README.md](sa_data/README.md) |
| [`bp_data/`](bp_data/) | DQI + BP | `python3 generate_bp_table.py` | [bp_data/README.md](bp_data/README.md) |
| [`fgum_data/`](fgum_data/) | Regev + FGUM | `python3 generate_fgum_table.py` | [fgum_data/README.md](fgum_data/README.md) |
| [`reproducibility/`](reproducibility/) | FGUM threshold certificate (Section 5 and Appendix A) | `python3 run_all.py` | [reproducibility/README.md](reproducibility/README.md) |
| [`qaoa_data/`](qaoa_data/) | QAOA (p ≤ 16) | — (angle table only) | [qaoa_data/README.md](qaoa_data/README.md) |
| [`qaoa_verifier/`](qaoa_verifier/) | QAOA verification (Julia) | `julia --project=. verify.jl` | [qaoa_verifier/README.md](qaoa_verifier/README.md) |
| [`qaoa_verifier_cpp/`](qaoa_verifier_cpp/) | QAOA verification (C++) | `make verify` | [qaoa_verifier_cpp/README.md](qaoa_verifier_cpp/README.md) |

Other directories: `src/` contains the C++ source, `test_data/` contains unit test fixtures.

## Verifying the FGUM Thresholds and Scores

The FGUM score generator uses the core-balance threshold from Section 5 and Appendix A for all fifteen integer pairs `3 <= k < D <= 8`. The separate interval certificate verifies the root enclosures, stability, and finite inequalities used in the appendix.

From the repository root:

```sh
python3 -m pip install -r reproducibility/requirements.txt
python3 reproducibility/verify_checksums.py
python3 reproducibility/run_all.py
```

The driver runs each parameter pair in a separate process and writes fresh results to `reproducibility/generated/`. The committed [`reproducibility/results/`](reproducibility/results/) directory contains the frozen results, including interval bounds and the score table. See [`reproducibility/README.md`](reproducibility/README.md) for the calculation, provenance, and commands for individual pairs.

## Verifying the QAOA Results

The QAOA satisfaction fractions are the most computationally demanding results in the paper. Pre-computed optimised angles for 254 (k, D, p) configurations are in `qaoa_data/`. Two independent verifiers are provided:

### Julia verifier (`qaoa_verifier/`)

Self-contained Julia implementation using Double64 arithmetic (~31 decimal digits). Handles p ≤ 13 on a 64 GB machine.

```sh
cd qaoa_verifier
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. verify.jl              # all 254 configurations
julia --project=. verify.jl --p 10       # only depth p=10
julia --project=. verify.jl --max-p 12   # depths p=1..12
```

See [qaoa_verifier/README.md](qaoa_verifier/README.md) for the mathematical details and API.

### C++ verifier (`qaoa_verifier_cpp/`)

C++ implementation with double-double precision, extracted from [QOKit](https://github.com/jpmorganchase/QOKit). Faster than the Julia verifier and handles all depths through p=16.

```sh
cd qaoa_verifier_cpp
make
./cpp/build/verify qaoa_data/            # all 256 configurations
./cpp/build/verify qaoa_data/ --p 10     # only depth p=10
./cpp/build/verify qaoa_data/ --p 14-16  # depths p=14..16
```

See [qaoa_verifier_cpp/README.md](qaoa_verifier_cpp/README.md) for build instructions and options.

Both verifiers agree to Δ < 10⁻¹² on all tested configurations.

## References

- N. Shutty et al., "Optimization Using Locally-Quantum Decoders," [arXiv:2604.24633](https://arxiv.org/abs/2604.24633) (2026).
- J. Basso, E. Farhi, K. Marwaha, B. Villalonga, L. Zhou, "The QAOA at High Depth for MaxCut…," [arXiv:2110.14206](https://arxiv.org/abs/2110.14206) (2022).
- E. Farhi, J. Goldstone, S. Gutmann, "A Quantum Approximate Optimization Algorithm," [arXiv:1411.4028](https://arxiv.org/abs/1411.4028) (2014).

## License

Except as otherwise noted, the code and data here is released under the MIT license, as described in [LICENSE.txt](LICENSE.txt).

- `qaoa_verifier_cpp/` is licensed under the **Apache License 2.0** (see [qaoa_verifier_cpp/LICENSE](qaoa_verifier_cpp/LICENSE)), extracted from [QOKit](https://github.com/jpmorganchase/QOKit) by JPMorgan Chase.
- `src/argparse.hpp` has its own license — see the file header.
