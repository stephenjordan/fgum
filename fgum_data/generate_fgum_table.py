#!/usr/bin/env python3
"""Compute the core thresholds and Regev+FGUM scores for 3 <= k < D <= 8.

This is a floating-point table generator, not the directed-interval
certificate for the rank-threshold theorem in Appendix A.
"""

from fractions import Fraction
from math import comb, expm1, fsum, log1p, sqrt


def one_minus_one_minus_e_to_D(e, D):
    if e == 1.0:
        return 1.0
    return -expm1(D * log1p(-e))


def _pi_D(z, D):
    """Probability of at least two successes in Bin(D, z)."""
    # The positive-term sum avoids cancellation when z is small.
    return fsum(comb(D, d) * z**d * (1 - z)**(D - d)
                for d in range(2, D + 1))


def _core_parameters(p, k, D):
    """Parameterize a nonzero support fixed point by p = e*a*y**(k-2)."""
    y = one_minus_one_minus_e_to_D(p, D - 1)
    a = one_minus_one_minus_e_to_D(y**(k - 1), D - 1)
    e = p / (a * y**(k - 2))
    return e, a, y


def _core_balance(p, k, D):
    """Core row count minus column count, both normalized by m/D."""
    e, _, y = _core_parameters(p, k, D)
    return fsum((e * _pi_D(y**(k - 1), D),
                 (k - 1) * _pi_D(p, D),
                 -D * p * y))


def _check_stability(p, e, a, y, k, D):
    """Numerically check stability of the two support recursions."""
    # Their Jacobian is [[0, j12], [j21, j22]]. Its spectral radius is
    # the positive eigenvalue below, since all entries are nonnegative.
    j12 = (D - 1) * (k - 1) * (1 - y**(k - 1))**(D - 2) * y**(k - 2)
    j21 = (D - 1) * (1 - p)**(D - 2) * e * y**(k - 2)
    j22 = (D - 1) * (1 - p)**(D - 2) * e * a * (k - 2) * y**(k - 3)
    spectral_radius = (j22 + sqrt(j22**2 + 4 * j12 * j21)) / 2
    if not (0 < e < 1 and 0 < spectral_radius < 1):
        raise ArithmeticError("The computed core balance is not a stable threshold")


def e_max(k, D):
    """Return the numerical core threshold for an integer pair 3 <= k < D <= 8.

    The manuscript identifies this core threshold with e_max for exactly
    these 15 pairs. No rank-threshold claim is made for other parameters.
    """
    if not (isinstance(k, int) and isinstance(D, int) and 3 <= k < D <= 8):
        raise ValueError("Supported parameters are integers satisfying 3 <= k < D <= 8")

    # At p=1 the balance is k-D < 0. Find a strictly positive lower
    # endpoint by halving, so the trivial zero-core solution p=0 is never
    # an endpoint or a candidate root. The positive balance root is on
    # the stable branch with the largest y.
    lo, hi = 0.5, 1.0
    for _ in range(60):
        if _core_balance(lo, k, D) > 0:
            break
        hi = lo
        lo /= 2
    else:
        raise ArithmeticError("Could not bracket the nonzero core balance")

    for _ in range(100):
        mid = (lo + hi) / 2
        if mid == lo or mid == hi:
            break
        if _core_balance(mid, k, D) > 0:
            lo = mid
        else:
            hi = mid

    # Solving for p, rather than y, avoids ill-conditioning when y is
    # extremely close to one. Appendix A uses the equivalent y parameter.
    p = (lo + hi) / 2
    e, a, y = _core_parameters(p, k, D)
    _check_stability(p, e, a, y, k, D)
    return e


def ihat_star(D):
    # Expected Hamming weight of a uniformly random element of I_D^*.
    numerator = Fraction(0, 1)
    denominator = Fraction(0, 1)
    for w in range(1, D // 2 + 1):
        count = Fraction(comb(D, w), 2 if D % 2 == 0 and w == D // 2 else 1)
        numerator += count * w
        denominator += count
    return float(numerator / denominator)


def regev_fgum(k, D):
    e = e_max(k, D)
    alpha_min = (1 - e) * (1 - 2 ** (1 - D))
    score = 1 - alpha_min * ihat_star(D) / D
    return e, alpha_min, score


if __name__ == "__main__":
    print(" k  D        e_max      alpha_min    Regev+FGUM  rounded")
    print("-- --  -----------  -------------  ----------  -------")
    for k in range(3, 8):
        for D in range(k + 1, 9):
            e, alpha_min, score = regev_fgum(k, D)
            print(f"{k:2d} {D:2d}  {e:11.8f}  {alpha_min:13.10f}  "
                  f"{score:10.8f}  {score:.4f}")
