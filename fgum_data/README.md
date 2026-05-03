# Regev + FGUM

**Reproduces:** Regev+FGUM column of Table 1 in [arXiv:2604.24633](https://arxiv.org/abs/2604.24633).

The FGUM satisfaction fraction is computed from closed-form equations in the manuscript. No C++ build, no instances, no cluster — just one Python script.

## Reproduce

```sh
python3 generate_fgum_table.py
```

This prints `e_max`, `alpha_min`, and the satisfaction fraction for every (k, D) with 3 ≤ k < D ≤ 8.

## How it works

1. Solve the nonlinear equation for `e_max(k, D)` by bisection (Eq. 29 in the paper).
2. Compute `alpha_min = (1 - e_max) * (1 - 2^(1-D))` (Eq. 30).
3. Compute the expected Hamming weight of a uniformly random element of `I*_D` (Eq. 34).
4. Satisfaction fraction = `1 - alpha_min * Î*_D / D` (Eq. 36).

Pure Python, standard library only.
