# QAOA Objective Verifier

Standalone C++ tool for verifying QAOA objective values (satisfaction fractions)
on D-regular Max-k-XORSAT, using the Basso/Villalonga branch tensor contraction
with double-double precision.

Companion to [arXiv:2604.24633](https://arxiv.org/abs/2604.24633).

Extracted from [QOKit](https://github.com/jpmorganchase/QOKit).

## Requirements

- C++17 compiler (GCC, Clang, or AppleClang)
- CMake 3.16+
- make

## Build & Run

```bash
make
./cpp/build/verify qaoa_data/
```

Or in one step:

```bash
make verify                       # double-double precision (default)
make verify PRECISION=float64     # faster, sufficient for p < 13
```

## Options

```
./cpp/build/verify qaoa_data/ --p 5              # only depth p=5
./cpp/build/verify qaoa_data/ --p 3-8            # depths p=3 through p=8
./cpp/build/verify qaoa_data/ --precision float64
./cpp/build/verify results_k3_D4.json            # single file
```

## Angle Convention

- Phase separator: `exp(-i * gamma * Z^k)`
- Mixer: `Rx(2*beta)`

## Input Format

JSON files (`results_k*_D*.json`):

```json
{
  "k": 3,
  "D": 4,
  "convention": "e^{-ig Z^k}",
  "results": {
    "1": {
      "gammas": [-0.275025258428398],
      "betas": [0.287371347852009],
      "objective": 0.676056660061
    },
    "2": { ... }
  }
}
```

Pass individual files or a directory containing them.
