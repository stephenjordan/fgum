#!/usr/bin/env python3
"""Compute Regev+FGUM satisfaction fractions for 3 <= k < D <= 8."""

from fractions import Fraction
from math import comb, expm1, log1p


def one_minus_one_minus_e_to_D(e, D):
    if e == 1.0:
        return 1.0
    return -expm1(D * log1p(-e))


def f(e, k, D):
    # We solve
    #   e = e/D + ((k-1)/D) * (1 - (1-e)^D).
    # After multiplying by D and dividing by e, the trivial root e=0 is removed.
    if e == 0.0:
        return (k - 1) * D - (D - 1)
    return (k - 1) * one_minus_one_minus_e_to_D(e, D) / e - (D - 1)


def e_max(k, D):
    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if f(mid, k, D) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


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
