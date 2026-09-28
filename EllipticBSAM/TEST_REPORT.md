# EllipticBSAM — Test Report

Solver: mass-conservative block-structured adaptive FAS multigrid
(Feng–Guo–Lowengrub–Wise, JCP 352 (2018) 463–497), generalized to
`-div(D(x,y) grad u) + C(x,y) u = f` with Dirichlet/Neumann boundaries.
All runs: MATLAB R2025b, double precision, headless (`matlab -batch`),
V(3,3) cycles, red–black Gauss–Seidel, QIF ghost interpolation unless
stated. Reproduce with `tests/unitKernels.m` and `tests/runAllTests.m`
(total suite ≈ 19–45 s on an Apple Silicon laptop). Raw results:
`output/results.mat`; figures in `output/`.
Re-verified on MATLAB R2026a (2026-07-01): all tests pass; numbers below
refreshed from that run. Supplementary stress tests (aspect-ratio
freedom, multi-parent children, anisotropic cells, mixed BCs, `rte`
indicator, 6-level hierarchies, `ebSample`) in `tests/auditExtras.m`.

**Status: ALL TESTS PASSED.**

## 0. Unit tests (`unit_kernels.m`)

| check | result |
|---|---|
| QIF ghost stencils reproduce global quadratics | 8.9e-16 (machine) |
| LIF ghost stencils reproduce global linears | 4.4e-16 (machine) |
| flux-form operator truncation rates (variable D, C) | 1.99, 2.00 |
| restrict∘prolong identity on affine data | 1.1e-16 |
| uniform 64² solve: rel. residual / error | 1.2e-12 / 8.3e-5 in 10 V-cycles |

## 1. T1 — uniform grids, variable D(x,y) and C(x,y), Dirichlet

`u = sin(πx)cos(πy)`, `D = 1 + 0.5 sin(2πx)cos(πy)`, `C = 1 + xy`.

| m | err L∞ | err L2 | V-cycle factor |
|---|---|---|---|
| 32 | 3.290e-04 | 1.278e-04 | 0.026 |
| 64 | 8.311e-05 | 3.243e-05 | 0.035 |
| 128 | 2.086e-05 | 8.185e-06 | 0.042 |
| 256 | 5.225e-06 | 2.057e-06 | 0.043 |

**L∞ rates: 1.98, 1.99, 2.00.**

## 2. T2 — pure-Neumann singular Poisson

`u = cos(2πx)cos(2πy)`, `C = 0`, all-Neumann (kernel = constants;
augmented coarsest solve + mean pinning).

L∞ errors 3.188e-03 → 5.019e-05 over m = 32…256;
**rates 1.99, 2.00, 2.00.**

## 3. T3 — fixed box-in-box hierarchies (paper Tables 2/3 setup)

Periodic Gaussian (paper eq. 5.1), `C = 1`, exact homogeneous Neumann
data on (−0.25, 0.25)²; nested refinement boxes 0.375 / 0.25 / 0.125
wide; 2-, 3-, 4-level meshes.

L∞ error (QIF):

| m0 | 2-level | 3-level | 4-level |
|---|---|---|---|
| 32 | 4.794e-03 | 1.177e-03 | 2.806e-04 |
| 64 | 1.212e-03 | 2.963e-04 | 6.856e-05 |
| 128 | 3.038e-04 | 7.428e-05 | 1.699e-05 |
| 256 | 7.599e-05 | 1.859e-05 | 4.231e-06 |

**All successive rates 1.98–2.03** (LIF is also second order; on the
deepest hierarchies its constants are up to ~2x larger, e.g. 5.44e-04
vs 2.81e-04 at m0=32, 4 levels).
Mass-conservation identity at the solved state: relative defect
**7e-17 … 3e-15 in every run** (24 runs). V-cycle factors 0.045–0.054
(8 cycles from scratch to rel. 1e-10 on every hierarchy).

## 4. T4 — L-shaped refinement patch (paper Table 5 setup)

Three-peak periodic Gaussian (paper eq. 5.2); level-2 region is the
L-shape `[0.375,0.625]×[0.375,0.5] ∪ [0.375,0.5]×[0.5,0.625]`, which the
Berger–Rigoutsos clusterer decomposes into 2 patches sharing an edge —
exercising same-level ghost exchange, the reentrant corner (a coarse
cell receiving **two** flux-balance corrections) and double zombies at
salient corners.

| m0 | err L∞ | factor | level-2 patches | mass defect |
|---|---|---|---|---|
| 16 | 5.670e-03 | 0.040 | 1 | 6.8e-15 |
| 32 | 1.640e-03 | 0.035 | 2 | 6.7e-13 |
| 64 | 3.999e-04 | 0.032 | 2 | 2.3e-11 |
| 128 | 9.839e-05 | 0.037 | 2 | 7.2e-12 |
| 256 | 2.439e-05 | 0.039 | 2 | 1.5e-12 |

(The T4 mass defects scale with the residual tolerance reached,
RelTol = 1e-10 here; tightening the tolerance drives them to 1e-16 as
in T7.)

**Rates: 1.79, 2.04, 2.02, 2.01** (the first interval is preasymptotic:
at m0=16 the min-width constraint yields a single patch).
Figures: `output/t4_lshape_mesh.png`, `output/t4_lshape_solution.png`.

## 5. T5 — fully adaptive refinement (undivided-Laplacian indicator)

Three sharp Gaussian bumps (σ = 0.02) on a 64² root; tagging by the
paper's eq. (4.32) indicator with per-level threshold `1e-2 / 4^(L-1)`,
buffer 2.

| MaxLevels | err L∞ | err L2 | visible dof | factor |
|---|---|---|---|---|
| 1 | 7.607e-02 | 3.324e-03 | 4096 | 0.028 |
| 2 | 1.950e-02 | 7.832e-04 | 5776 | 0.030 |
| 3 | 4.774e-03 | 1.941e-04 | 10399 | 0.038 |
| 4 | 1.269e-03 | 6.195e-05 | 25528 | 0.039 |

**Error drops ×3.8–4.1 per added level** (matching the h² expectation
for a feature fully captured by the next level), at a small fraction of
the 512²-uniform cell count (25.5k vs 262k).
Figures: `output/t5_adaptive_mesh.png`, `output/t5_adaptive_solution.png`.

## 6. T6 — h-independent geometric solver convergence (paper Fig. 8)

Average per-cycle composite-residual reduction factor (solving the
composite system from a zero start), 2-level mesh:

| m0 | 32 | 64 | 128 | 256 |
|---|---|---|---|---|
| factor | 0.045 | 0.050 | 0.052 | 0.053 |

i.e. **~1.3 digits per V(3,3) cycle, independent of h** — matching the
paper's reported factor ≈ 1/20 (their Fig. 8). On 3- and 4-level meshes
the factors are the same 0.05 ± 0.005: independent of depth as well.
History plot: `output/t6_vcycle_history.png`.

Implementation note: physical boundary conditions are *folded exactly*
into the smoother (partial ghosts + diagonal correction, see
`ebBcFold.m` and Lemma 6.2 of `docs/convergence_theory.pdf`); without
the fold, the lagged Dirichlet ghost degrades the cycle from ~10 to ~22
V-cycles for 1e-12.

## 7. T7 — discrete mass conservation and the role of the corrections

Same problem as T3 (m0 = 64). "Corrections off" = plain FAS restriction
of the solution as the coarse interface neighbor (the paper's
non-conservative "QI" variant).

| levels | corrections | mass defect (rel) | err L∞ |
|---|---|---|---|
| 2 | on | 1.4e-16 | 1.212e-03 |
| 2 | off | 2.3e-05 | 1.208e-03 |
| 3 | on | 5.2e-16 | 2.963e-04 |
| 3 | off | 2.9e-04 | 3.402e-04 |

The flux-balance corrections keep the conservation identity at machine
precision **and** protect accuracy: without them the interface defect
grows by an order of magnitude per level pair (2.3e-05 → 2.9e-04) and
the 3-level error is 15% worse in this configuration; on hierarchies
whose refinement boxes hug the solution features more tightly the
accuracy gap has reached 5.6× — the multilevel mechanism behind the
divergence of earlier in-house attempts.

## 8. T8 — robustness against pathological clustering

24 random tag disks (the adversarial scenario that broke previous
implementations: small/irregular clustered patches). 3-level adaptive,
MinBoxWidth 4.

Result: **49 patches, minimum width 4 fine cells, V-cycle factor 0.052,
err 8.38e-05, mass defect 1.1e-13** — no divergence, no residual floor.
Sliver boxes left over after the proper-nesting clip are merged into
neighbors or dropped (warning `eb_tag_cluster:dropthin`) instead of ever
producing under-sized patches. Mesh: `output/t8_scattered_mesh.png`.

## Why this implementation does not reproduce the old failure modes

1. **Conservative coarse–fine coupling everywhere** (QIF/LIF ghosts +
   flux-balanced zombie corrections at *every* level pair, accumulated
   correctly at reentrant corners) — the old "block" backend lacked
   refluxing and the AFAS path lost it during post-sync.
2. **Patch geometry discipline**: patches always aligned to parent
   cells (odd/even index alignment), ≥ 2 parent cells wide (target
   MinBoxWidth 4–8), properly nested with full 8-neighborhood parent
   coverage, 1-conforming. The clusterer cannot emit slivers.
3. **Composite residual** as the convergence control (residual of the
   actual conservative composite system, corrections included), so
   stalls are visible rather than hidden by a fine-level-only norm.
4. **Vectorized kernels + precomputed index transfers**: all hot-path
   operations are array ops; coefficient evaluation is array-in/array-out
   (the old per-cell evaluation was the documented 100× slowdown).
