# QAOA Satisfaction Fraction Verifier

Self-contained verification tool for the QAOA results in [arXiv:2604.24633](https://arxiv.org/abs/2604.24633).

Given pre-computed QAOA angles (γ, β) and problem parameters (k, D), computes the exact satisfaction fraction c̃(γ, β) for depth-p QAOA on D-regular Max-k-XORSAT.

**Author:** John S Azariah — Centre for Quantum Software and Information, UTS

---

## Quick Start

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'    # one-time setup
julia --project=. verify.jl                             # verify all 254 configurations
julia --project=. -e 'using Pkg; Pkg.test()'            # run 42 unit tests
```

---

## What This Computes

### The QAOA circuit on the light-cone tree

For a D-regular Max-k-XORSAT instance, the depth-p QAOA circuit's expected satisfaction of any single clause depends only on the local structure within the light cone — a tree with branching factor b = (D−1)(k−1) extending p levels from the root clause. The evaluation avoids constructing the exponentially large quantum state and instead uses a **branch tensor recurrence** that contracts the tree from leaves to root.

### The branch tensor recurrence

The algorithm proceeds in three stages:

**Stage 1 — Leaf initialisation.** The branch tensor B⁽⁰⁾ is the all-ones vector over 2^(2p+1) configurations. Each configuration is a (2p+1)-bit string encoding ket and bra indices at each depth level plus a root spin.

**Stage 2 — Iterative contraction** (t = 1, …, p):
1. **Child weights:** w(a) = f(a) · B^(t−1)(a), where f(a) is the mixer weight — a product of cos(β) and i·sin(β) factors from the bit transitions in the branch string.
2. **Constraint fold:** ŵ = WHT(w), then iWHT(κ̂ · ŵ^(k−1)), where κ(a) = cos(Γ·spins(a)/2) and WHT is the Walsh-Hadamard transform on (Z₂^n, ⊕).
3. **Branching power:** B^(t)(a) = [fold(a)]^(D−1).

Per-step max-magnitude normalisation with log-scale accumulation prevents Float64 overflow at high (k, D, p).

**Stage 3 — Root extraction.** Form root message m(a) = (−1)^{a⁰} · f(a) · B^(p)(a), apply a k-fold XOR convolution, contract against root kernel κ_r(a) = i·sin(c_s · Γ·spins(a)/2), and compute:

> c̃ = (1 + c_s · ⟨Z₁···Zₖ⟩) / 2

where c_s = +1 for XORSAT and c_s = −1 for MaxCut.

### Computational cost

| Depth | Space | Time | Wall clock (approx.) |
|---|---|---|---|
| p ≤ 8 | 256 KB | O(4^8 · 8²) | < 1 second |
| p = 12 | 16 MB | O(4^12 · 12²) | ~1 minute |
| p = 16 | 4 GB | O(4^16 · 16²) | ~30 minutes |

---

## API

```julia
using QaoaVerifier

# Max-3-XORSAT on 4-regular, depth 1
c̃ = verify_xorsat(3, 4, [0.550], [0.287])

# MaxCut on 3-regular, depth 1
c̃ = verify_maxcut(3, [0.615], [0.393])

# General interface
c̃ = verify_satisfaction_fraction(k, D, γ, β; clause_sign=±1)
```

---

## Source Overview

The entire verifier is a single file: `src/QaoaVerifier.jl` (~400 lines).

| Section | Lines | What it does |
|---|---|---|
| Types | ~20 | `TreeParams(k,D,p)`, `QAOAAngles(γ,β)` |
| Walsh-Hadamard Transform | ~60 | Cache-oblivious recursive WHT with iterative SIMD kernel |
| Branch Tensor Recurrence | ~200 | Mixer table, constraint kernel, WHT-accelerated fold, normalised iteration |
| Public API | ~30 | `verify_satisfaction_fraction`, `verify_maxcut`, `verify_xorsat` |

---

## References

- Basso, Farhi, Marwaha, Villalonga, Zhou, "The Quantum Approximate Optimization Algorithm at High Depth for MaxCut on Large-Girth Regular Graphs and the Sherrington-Kirkpatrick Model," [arXiv:2110.14206](https://arxiv.org/abs/2110.14206) (2022).
- Farhi, Goldstone, Gutmann, "A Quantum Approximate Optimization Algorithm," [arXiv:1411.4028](https://arxiv.org/abs/1411.4028) (2014).
