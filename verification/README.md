# Verification

Verified September 28, 2026 with MATLAB 26.1.0.3276743 (R2026a) Update 3 on macOS / Apple
silicon, using two MATLAB computational threads per process.

## Executed checks

| Package | Executed checks | Result |
|---|---|---|
| QC-ring elliptic | `testEbNeumannSemantics`, `runAllTests` (T1–T8) | PASS |
| Strict-tree elliptic | `testEbtNeumannSemantics`, `runAllTreeTests` (T1–T9) | PASS |
| Parabolic | `runParabolicRegressions`, `runParabolicTests` (P1–P4) | PASS |
| Diffuse domain | `runDdRegressions`, `runDdTests` (D1/D1b/D2/D3) | PASS |
| New starfish | `testStarfish`, `demoStarfish('full')` | PASS |
| README snippets | QC-ring, tree, and heat examples | PASS |

The suites were also rerun through the public `runVerification` entry point
to verify path isolation. [Suite output](suite-output.txt),
[starfish output](starfish-output.txt), and [README output](readme-examples.txt)
record the executed results. Generated MAT-files and test figures are local
outputs and are excluded from version control.

Selected checks from the fresh run:

- Nonautonomous Crank–Nicolson reaction test: temporal orders approximately
  2.001, 2.000, 2.000.
- Joint space/time BDF2 refinement with AMR: orders approximately 1.83, 1.86.
- Heat-equation mass drift across 39 regrids: approximately 1.09e-10.
- Tree tests: single-parent containment and a feature crossing two parents.
- Starfish distance: 300 deterministic queries checked against exhaustive
  segment distances and `inpolygon`; boundary vertices checked separately.
- Constant transmission solution: maximum error approximately 4.26e-14.

## Full starfish result

| Full starfish run | Measured result |
|---|---:|
| Visible composite cells | 408,562 |
| Uniform grid at the same finest spacing | 1,048,576 |
| Reduction in active cell count | 61.0% |
| Stored AMR cells, including covered coarse cells | 543,384 |
| V-cycles | 10 |
| Relative composite residual | 2.55e-11 |
| Relative L2 difference from uniform 1024² | 1.06e-04 |
| Maximum difference from restricted uniform 1024² | 6.72e-04 |
| Uncovered cell centers in the tested interface band | 0 |

Cell counts measure refinement savings; they are not a wall-time or memory
speedup claim. The comparison is at fixed epsilon against a numerical
reference, with that reference averaged onto the visible AMR cells.

The uniform 512-to-1024 relative L2 difference, after restriction to 512,
was 3.959e-05. The AMR integrated equation
balance defect was 4.866e-15.

[JSON metrics](starfish-full.json) and [residual history](starfish-full-residual.csv)
are produced by the driver. The driver asserts relative residuals below
1e-10, complete sampled coverage in the tested interface band, and an
AMR/uniform relative L2 difference below 1e-3. Timings include setup and
geometry evaluation, depend on hardware/load, and are not a performance
benchmark.

## Reproduction

From the repository root in MATLAB:

```matlab
runVerification
bsamSetup
result = demoStarfish('full');
```

Or with the MATLAB executable on your PATH:

```sh
matlab -batch "maxNumCompThreads(2); runVerification"
matlab -batch "maxNumCompThreads(2); bsamSetup; demoStarfish('full');"
python3 scripts/verify_package.py
```

Each legacy test suite can also run from its own `tests` folder. Avoid
adding all folders recursively: both elliptic suites have `mmsCases.m`.

## Scope and numerical details

The tests exercise the listed cases. The additional geometry/forcing sweeps,
stress tests, and long movie examples included in the copied packages were
not rerun as part of this delivery. Historical `TEST_REPORT.md` files record
other runs and are not fresh verification of every optional example.

Some clustering runs warn that thin clusters are dropped. The starfish
coverage check independently confirms zero missing finest-level cell
centers in `abs(d)<epsilon` for the delivered configuration. It does not
establish coverage for all parameter choices.

The starfish is the two-sided transmission example at a fixed epsilon.
There is no claimed exact sharp-interface solution for its active forcing.
The mathematical documentation also records qualifications to older QIF
matrix-sign arguments; numerical test results do not repair those proofs.

All three gallery PNGs were visually inspected. Their corresponding PDFs
are included for download; the dense 3D surface is rasterized in its PDF to
keep the file small, while hierarchy and validation plots use vector export.
