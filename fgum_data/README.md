# Regev + FGUM

**Reproduces:** Regev+FGUM column of Table 1 in [arXiv:2604.24633](https://arxiv.org/abs/2604.24633).

The FGUM satisfaction fraction is computed from the core-balance threshold
in Section 5 and Appendix A of the manuscript. The table generator uses
only the Python standard library; no C++ build or sampled instances are
needed.

## Reproduce

```sh
python3 generate_fgum_table.py
```

This prints `e_max`, `alpha_min`, and the satisfaction fraction for every
integer pair `(k, D)` with `3 <= k < D <= 8`. The last column gives the
satisfaction fraction rounded to four decimal places, as in Table 1.

## How it works

The support probabilities satisfy the manuscript's core support equations.
The script parameterizes their nonzero solutions by `p`:

```text
phi_D(z) = 1 - (1-z)^(D-1)
pi_D(z)  = Pr[Bin(D,z) >= 2]
y        = phi_D(p)
a        = phi_D(y^(k-1))
e        = p / (a * y^(k-2))
```

The limiting core row and column counts, both normalized by `m/D`, are

```text
rows    = e * pi_D(y^(k-1)) + (k-1) * pi_D(p)
columns = D * p * y
```

The script finds the positive solution of `rows - columns = 0` by
bisection. Its bracket excludes `p=0`, which would describe an empty core.
The positive balance root lies on the stable branch with the largest
support probabilities. Parameterizing by `p` avoids numerical loss when
`y` is close to one. The script also checks numerically that the support
recursion's Jacobian has spectral radius below one at the solution.

For the 15 supported pairs, the manuscript's rank-threshold theorem
identifies this core threshold with `e_max`. The remaining steps are:

1. Compute `alpha_min = (1 - e_max) * (1 - 2^(1-D))`.
2. Compute `ihat_star(D)`, the average Hamming weight of a uniformly random
   element of `I*_D`, using exact rational sums. For even `D`, the set
   includes half of the strings of weight `D/2`.
3. Compute the limiting satisfaction fraction as
   `1 - alpha_min * ihat_star(D) / D`.

As in Section 5, the score is the limiting value approached from below
the erasure threshold. It does not assert vanishing decoding failure at
the threshold itself.

## Scope

This script generates numerical table values using ordinary
floating-point arithmetic. Its stability check is numerical; it does not
perform the directed-interval certificate described in Appendix A.

The functions `e_max(k, D)` and `regev_fgum(k, D)` support only the 15
integer parameter pairs above and raise `ValueError` outside that range.
The command-line output and these function interfaces can also be used
independently of a certificate run.
