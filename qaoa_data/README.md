# QAOA Optimised Angles

**Reproduces:** QAOA column of Table 1 in [arXiv:2604.24633](https://arxiv.org/abs/2604.24633).

The file `qaoa_angles_for_verifier.csv` contains 254 optimised (k, D, p) configurations — 16 depths for MaxCut (k=2, D=3) and up to 16 depths for each XORSAT (k, D) pair with 3 ≤ k < D ≤ 8.

## CSV format

```
k,D,p,ctilde,source,gamma,beta
```

- `k` — constraint arity (2 for MaxCut, ≥3 for XORSAT)
- `D` — variable degree
- `p` — QAOA depth
- `ctilde` — claimed satisfaction fraction
- `source` — provenance label
- `gamma`, `beta` — semicolon-separated angle vectors of length p

## How the angles were computed

The angles were optimised using gradient-based methods on GPUs, independently by two implementations:

- [QOKit](https://github.com/jpmorganchase/QOKit/tree/add-max-k-xor-sat/qokit/max_k_xor_sat) — Primary implementation, C++, double-double precision (JP Morgan Chase).
- [QaoaXorsat.jl](https://github.com/johnazariah/qaoa-xorsat) — Independent cross-validation of the angle optimisation and development of the verification pipeline (Julia).

## Verification

The satisfaction fractions can be independently verified using [`qaoa_verifier/`](../qaoa_verifier/):

```sh
cd ../qaoa_verifier
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. verify.jl
```

See [qaoa_verifier/README.md](../qaoa_verifier/README.md) for details.
