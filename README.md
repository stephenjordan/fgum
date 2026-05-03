# Optimization Using Locally-Quantum Decoders

Companion code and data for:

> N. Shutty, A. Mandal, S. Ragavan, Q. Buzet, A. Chailloux, N. C. Rubin,
> A. Khan, S. Boulebnane, R. Shaydulin, J. Azariah, S. P. Jordan,
> "Optimization Using Locally-Quantum Decoders,"
> [arXiv:2604.24633](https://arxiv.org/abs/2604.24633) (2026).

Everything needed to reproduce Table 1 of the paper is here. Each column of the table has its own directory with scripts, data, and a README explaining how to run it.

## Quick Start

```sh
make          # build C++ tools (SA solver, instance generator, density evolution)
make check    # run unit tests
```

Requires a C++20 compiler and GNU Make. The Python scripts use only the standard library. The QAOA verifier requires [Julia 1.9+](https://julialang.org/downloads/).

## Table 1 — Directory Guide

| Directory | Table 1 Column | What to run | Details |
|---|---|---|---|
| [`sa_data/`](sa_data/) | Simulated Annealing | `python3 generate_sa_table.py` | [sa_data/README.md](sa_data/README.md) |
| [`bp_data/`](bp_data/) | DQI + BP | `python3 generate_bp_table.py` | [bp_data/README.md](bp_data/README.md) |
| [`fgum_data/`](fgum_data/) | Regev + FGUM | `python3 generate_fgum_table.py` | [fgum_data/README.md](fgum_data/README.md) |
| [`qaoa_data/`](qaoa_data/) | QAOA (p ≤ 16) | — (angle table only) | [qaoa_verifier/README.md](qaoa_verifier/README.md) |
| [`qaoa_verifier/`](qaoa_verifier/) | QAOA verification | `julia --project=. verify.jl` | [qaoa_verifier/README.md](qaoa_verifier/README.md) |

Other directories: `src/` contains the C++ source, `test_data/` contains unit test fixtures.

## Verifying the QAOA Results

The QAOA satisfaction fractions are the most computationally demanding results in the paper. Pre-computed optimised angles for 254 (k, D, p) configurations are in `qaoa_data/`. A self-contained Julia verifier independently recomputes each value:

```sh
cd qaoa_verifier
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. verify.jl
```

All 254 configurations verify to Δ < 10⁻¹² against the claimed values. See [qaoa_verifier/README.md](qaoa_verifier/README.md) for the mathematical details and API.

## References

- N. Shutty et al., "Optimization Using Locally-Quantum Decoders," [arXiv:2604.24633](https://arxiv.org/abs/2604.24633) (2026).
- J. Basso, E. Farhi, K. Marwaha, B. Villalonga, L. Zhou, "The QAOA at High Depth for MaxCut…," [arXiv:2110.14206](https://arxiv.org/abs/2110.14206) (2022).
- E. Farhi, J. Goldstone, S. Gutmann, "A Quantum Approximate Optimization Algorithm," [arXiv:1411.4028](https://arxiv.org/abs/1411.4028) (2014).

## License

MIT. See individual source files for third-party licenses (`src/argparse.hpp`).
