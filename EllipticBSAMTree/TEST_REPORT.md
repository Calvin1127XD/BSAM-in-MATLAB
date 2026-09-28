# EllipticBSAMTree — Test Report

Strict single-parent tree variant of the mass-conservative
block-structured adaptive FAS multigrid (Feng–Guo–Lowengrub–Wise, JCP
352 (2018) 463–497), tree design after Fortran BSAM 2.0 (SuperBSAM).
All runs: MATLAB R2026a, double precision, headless (`matlab -batch`),
V(3,3) cycles, red–black Gauss–Seidel, QIF ghost interpolation unless
stated. Reproduce with `tests/runAllTreeTests.m` (suite ≈ 37 s on an
Apple Silicon laptop).

**Status: ALL TESTS PASSED (T1–T9).**

## T1–T8 — parity with EllipticBSAM

The tree refactor reproduces the EllipticBSAM suite essentially
verbatim (T1–T7 numbers are identical to the QC-ring implementation;
T8's mesh differs slightly because clustering is per parent — 50 vs 49
patches — with matching quality):

| test | key result |
|---|---|
| T1 uniform, variable D,C, Dirichlet | L∞ rates 1.98, 1.99, 2.00; factors 0.026–0.043 |
| T2 pure-Neumann singular | rates 1.99, 2.00, 2.00 |
| T3 box-in-box 2/3/4-level (QIF+LIF) | all rates 1.98–2.03; mass 7e-17…3e-15; factors 0.045–0.054, V=8 |
| T4 L-shaped patch | rates 1.79 (preasymptotic), 2.04, 2.02, 2.01; 2 patches from m0=32 |
| T5 fully adaptive (ulap) | error ×0.25–0.27 per added level; 4 levels: 1.269e-03 at 25.5k dof |
| T6 h-independence | factors 0.045–0.053 for m0 = 32…256 |
| T7 conservation | corrections on: defect 1.4e-16 / 5.2e-16; off: 2.3e-05 / 2.9e-04 |
| T8 scattered-tag robustness | 50 patches, min width 4, factor 0.052, err 8.37e-05, **mass 1.1e-13** |

## T9 — tree invariants (new)

* `ebtTreeCheck` passes on the fully adaptive T5 hierarchy
  (4 AMR levels, patches [1 3 3 3], max 3 children per patch):
  single-parent containment, box/mbounds consistency, 2-alignment,
  disjoint siblings, children/covered bookkeeping, ring nesting.
* Straddle scenario (a tagged feature centered on the shared edge of two
  level-2 patches): per-parent clustering yields **two level-3 children,
  parents [1 2]** — never a multi-parent child. Solver: V=9,
  rel 3.6e-12, factor 0.033, mass 9.5e-13.

## Regression found and fixed during bring-up

Initial T8 runs showed a mass defect of 9.0e-08 that did **not** shrink
with tighter `RelTol` (a structural inconsistency, not a solver-accuracy
artifact). Cause: coarse patch boundary faces are stored twice (one copy
per abutting patch). When a child's footprint touches its parent's edge,
the fill-down averaged face diffusivity was written only into the
parent's copy, so the abutting sibling computed the flux across the same
physical face with the *analytic* D — breaking the flux telescoping for
variable-coefficient problems (only T8 uses variable D on an adaptive
mesh, which is why the other tests stayed clean; a ring-vs-owner-truth
comparison isolated the coefficient path by showing the ghost/ring
values were exact). Fix: `ebtCoeffFilldown>pushFaceX/pushFaceY` routes
boundary-face averages to the abutting siblings via the owner map (the
EllipticBSAM `rectIsect`-over-all-patches behavior, restricted to the
shared faces). After the fix: T8 mass = 1.1e-13.
