# fgum

Software to reproduce the results in "Optimization Using Locally-Quantum Decoders" [arXiv:2604.24633].

## Build

This curated archive builds with GNU Make and a C++20 compiler. No Bazel installation is required.

```sh
make
```

This builds the following executables in `bin/`:

- `test_xortools`
- `full_anneal`
- `random_regular`
- `density`

To run the bundled bitmatrix/xortools test suite and a small density-evolution smoke test:

```sh
make check
```

To remove built binaries and generated object/dependency files:

```sh
make clean
```

## Notes on this curated version

The original density-evolution source used Bazel-style include paths such as
`argparse/argparse.hpp`, `belief/density.hpp`, and `irreg/distribution.hpp`.
Those include paths have been converted to the flat `src/` layout used in this
archive so the entire included C++ codebase can be built with the Makefile.
