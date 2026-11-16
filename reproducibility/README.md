# FGUM core thresholds and interval certificates

This directory reproduces the FGUM thresholds and satisfaction fractions for
all fifteen integer pairs `3 <= k < D <= 8`, together with the finite
inequalities used in Section 5 and Appendix A of *Optimization Using
Locally-Quantum Decoders*.

The implementation was reconstructed from the manuscript's equations because
the supplied repository and its available Git history did not contain the
ancillary certificate package. All frozen results here were freshly generated
and verified by this implementation. They are not recovered outputs from an
earlier calculation.

## Requirements

Use Python 3.10 or later and `mpmath==1.3.0`. The committed results were generated
with Python 3.12.14, mpmath 1.3.0, and 80 decimal digits of interval precision.

From the repository root:

```sh
python3 -m pip install -r reproducibility/requirements.txt
```

The interval programs enforce the pinned mpmath version. The ordinary numerical
score generator in `fgum_data/` requires only Python's standard library.

## Verify the frozen results

```sh
python3 reproducibility/verify_checksums.py
python3 reproducibility/run_all.py --verify-only --output reproducibility/results
```

The first command checks the frozen files and the relevant source files against
`SHA256SUMS`. The second recomputes every saved interval inequality, checks that
the interval covers are complete, and verifies both CSV tables against the
certificate outputs and manuscript values. It does not modify the frozen
results.

These commands serve different purposes: a matching checksum establishes file
integrity; replay establishes that the saved certificate inputs pass the
mathematical checks.

## Regenerate all fifteen pairs

```sh
python3 reproducibility/run_all.py
```

From this directory, the equivalent command is `python3 run_all.py`.

Fresh results go into `reproducibility/generated/`, which Git ignores. Each pair
is generated in its own process, then replayed in another process. Up to four
pairs run concurrently; use `--jobs 1` for serial execution or `--jobs N` to set
the limit. The driver exits unsuccessfully if any calculation, replay, or table
comparison fails. Its run manifest records `running`, `failed`, or `passed`
explicitly.

```sh
python3 reproducibility/run_all.py --jobs 1 --output reproducibility/generated
```

To intentionally replace the frozen results after a source update:

```sh
python3 reproducibility/run_all.py --output reproducibility/results
python3 reproducibility/verify_checksums.py --write
python3 reproducibility/verify_checksums.py
```

The checksum-writing option is for maintainers making an intentional update;
it is not part of checking an existing release.

## One pair, scores, and tables

```sh
python3 reproducibility/certify_core.py --k 3 --D 4 \
    --output reproducibility/generated/certificate_k3_D4.json
python3 reproducibility/certify_core.py \
    --verify reproducibility/generated/certificate_k3_D4.json

python3 fgum_data/generate_fgum_table.py
python3 reproducibility/generate_tables.py --input reproducibility/results --check
```

The score generator provides ordinary floating-point values. The CSV score
table additionally records interval enclosures derived from the certified
balance root and exact rational Hamming-weight averages.

## Files and outputs

| File or directory | Contents |
|---|---|
| `certify_core.py` | Root, stability, branch, small-support, exterior, and central interval checks; certificate replay |
| `run_all.py` | Separate-process driver for all fifteen pairs, including replay and table comparisons |
| `generate_tables.py` | CSV generation and checks of the manuscript's displayed precision and conservative bounds |
| `reference_values.json` | Manuscript comparison targets; these are inputs, not computed certificate outputs |
| `MATHEMATICS.md` | Formulas, analytic branch selection, and explanation of interval coverage |
| `results/certificate_kK_DD.json` | One frozen certificate for each parameter pair, including every interval cell and its Chernoff tilts |
| `results/certificate_summary.csv` | Computed interval bounds alongside the conservative bounds displayed in the manuscript |
| `results/fgum_scores.csv` | Threshold, measurement, and score intervals; exact Hamming-weight averages; printed scores |
| `results/run_manifest.json` | Completion status, runtime versions, and source hashes |
| `SHA256SUMS`, `verify_checksums.py` | Frozen package integrity manifest and checker |
| `test_certificate.py` | Replay acceptance and deliberate-corruption rejection checks |

Run the replay regression checks with:

```sh
python3 -m unittest discover -s reproducibility -p 'test_certificate.py'
```

## Interpreting the certificate

Root brackets have radius `1e-45`. The program verifies strict opposite signs
at their endpoints and a negative balance derivative throughout each bracket.
It then certifies complete interval covers, not a finite sampling of the
functions. Numerical root finding and binary64 tilt optimization only propose
inputs; every accepted inequality uses directed interval arithmetic.

All exported decimal interval endpoints are rounded outward from the exact
binary interval endpoints. Cover endpoints are exact rational numbers. Replay
recomputes the bounds from the saved inputs rather than trusting saved results.

The reconstructed calculation verifies all the manuscript's displayed
`U_max` and central-curvature bounds. Its actual interval bounds can be
stronger, because the adaptive subdivision differs. Both the computed bounds
and the retained conservative manuscript bounds are recorded.

Nine spectral-radius entries needed an upward change of `0.000001` to make
their six-decimal values valid upper bounds. The accompanying manuscript and
`reference_values.json` contain the corrected bounds. The threshold and
satisfaction-fraction values are unchanged.

This package verifies the finite numerical inequalities. The analytic
peeling, conditional distribution, concentration, and first-moment arguments
are supplied by the paper; the program is not a formal proof checker for
those arguments.
