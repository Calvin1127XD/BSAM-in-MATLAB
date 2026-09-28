# EllipticBSAMTree — strict-tree block-structured adaptive FAS multigrid

A tree-based variant of **EllipticBSAM** (same discretization, same
conservative AFAS multigrid) in which the patch hierarchy is organized as
a **strict single-parent tree**, faithfully following the design of the
Fortran BSAM 2.0 library (`treeops.f90` / `bsamroutines.f90` of
SuperBSAM):

```
-div( D(x,y) grad u ) + C(x,y) u = f(x,y)   on [a,b] x [c,d],
```

with `D > 0`, `C >= 0`, Dirichlet and/or Neumann data per side.

## What "strict tree" changes vs. EllipticBSAM

| | EllipticBSAM (QC-ring) | EllipticBSAMTree |
|---|---|---|
| clustering | Berger–Rigoutsos on the whole level | Berger–Rigoutsos **inside each parent patch** (Fortran `RefineGrid` → `NewSubGrids`) |
| lineage | `parents{p}` may list several coarse patches | `parent(p)` is exactly **one** patch; `mbounds(p,:)` records the footprint in that parent's cell space |
| a feature crossing two coarse patches | one straddling child | **two abutting children**, one per parent, coupled as same-level siblings |
| parent window / QC ring | gathered per cell from all owner patches | **one contiguous slice of the parent's padded array** (ring cells that leave the parent's interior read the parent's ghost cells, synchronized at the parent's level) |
| restriction / FAS load / coefficient fill-down | scattered over all intersecting coarse patches | single sliced write into the parent over `mbounds` |
| flux-balance correction targets | owner-map scatter | unchanged (targets may belong to the parent's siblings — the precomputed analogue of Fortran's `masscorrlayer` strips) |

Two index systems per patch, exactly as in the Fortran code:
`lev.box(p,:)` (global at the patch's level — Fortran `mglobal`; drives
same-level exchange/owner maps) and `lev.mbounds(p,:)` (parent-local —
Fortran `mbounds`; drives all coarse-fine transfer).

Numerics are unchanged: QIF/LIF ghost interpolation, flux-balanced
zombie corrections, red-black Gauss–Seidel with global parity,
conservative 4-point restriction, bilinear prolongation, exact boundary
folding. Convergence and conservation match EllipticBSAM (see
`TEST_REPORT.md`).

One subtlety worth knowing (found the hard way): coarse patch boundary
**faces are stored twice** (in both abutting patches' face arrays). When
a child's footprint touches its parent's edge, the coefficient fill-down
must push the averaged face `D` into the abutting sibling's copy too
(`ebtCoeffFilldown>pushFaceX/Y`), or discrete conservation breaks at
variable-coefficient interfaces.

Proper nesting is enforced *pre-emptively* by clipping tags to cells
whose full 8-neighborhood is covered by the current finest level
(`erode8`), which guarantees the Fortran 2-level-conformity invariant
without the `levellandscape`/defect-restart machinery. Tag *buffering*
is global (crosses patch boundaries) exactly like Fortran
`BufferAndList`/`InflateEdgeTags`; only the clustering is per-parent.

## Layout

```
src/        solver source (ebt*.m), entry point ebtSolve.m,
            tree validator ebtTreeCheck.m
examples/   demoTreeBasicUniform.m, demoTreePoissonAdaptive.m
tests/      unitTreeKernels.m, smokeTreeAmr.m, runAllTreeTests.m (T1-T9)
```

## Quick start

```matlab
addpath EllipticBSAMTree/src
uex  = @(x,y) sin(2*pi*x) .* sin(2*pi*y);
f    = @(x,y) (1 + 8*pi^2) * uex(x,y);              % -Lap u + u = f
prob = ebtProblem([0 1 0 1], 1, 1, f, {'dirichlet', 0}, uex);
sol  = ebtSolve(prob, 'BaseCells', [64 64], 'MaxLevels', 4, ...
                'TagThreshold', 1e-3);
sol.err                       % composite error vs exact
sol.mass                      % discrete conservation check
ebtTreeCheck(sol.H)           % validate the tree invariants
ebtPlot(sol.H, 'mesh', 'mesh.png');
```

Options are identical to `ebOptions` (see `ebtOptions`).

## Verification

`tests/runAllTreeTests.m`: T1–T8 are the EllipticBSAM suite (same
manufactured problems; the tree refactor reproduces the same errors,
rates 1.98–2.03, V-cycle factors 0.026–0.054, conservation to 1e-13);
T9 validates the tree invariants (`ebtTreeCheck`) and asserts that a
feature crossing two parents produces per-parent children with exactly
one parent each.

Full write-up: `EllipticBSAMTree_Documentation.pdf` in this folder (the
strict tree, its two index systems, per-parent clustering, and the
single-parent transfer machinery).

The documentation now includes compatibility, error-norm and proof-scope
clarifications. Only all-Neumann `C==0` problems have a constant gauge;
`sol.H.compatibility.sourceShift` records the constant removed from their
sampled source after each refinement. Positive reaction retains absolute
error. The multidimensional QIF matrix has positive off-diagonal entries,
so the companion theory's blanket matrix-sign argument needs revision.
See `tests/testEbtNeumannSemantics.m` for the compatibility regressions.
The optional reinforcement-learning experiments from the research checkout
are outside this solver collection.
