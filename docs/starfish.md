# Starfish benchmark

## Reference

The geometry, box, epsilon, coefficients, and forcing follow
[`DiffuseDomainHelmholtz`, commit `a9d8449`](https://github.com/Calvin1127XD/DiffuseDomainHelmholtz/tree/a9d8449005f4ab030e64cf4c346c7225998fa01e/DiffuseDomainMultigrid/2DCC).
The relevant files are `generateBoundary.m` and `elliptic2DMGCCDiff.m`.
The active code specifies the five-armed starfish and the forcing in this
repository's README. Some comments still mention an ellipse. In particular,
the ellipse reference field in that script is not an exact solution of its
active starfish problem; it is not used to report an error here.

## Geometry and phase field

`ddStarfishGeometry(3000)` constructs the same boundary polygon, including
the repeated endpoint. Its signed distance is positive inside. A tree of
segment bounding boxes computes the nearest-segment distance directly;
the original repository instead propagates near-boundary distances with
an eikonal/Fast Marching method. Thus the geometry and PDE parameters match,
while the numerical distance construction differs.

The polygon is an approximation of the smooth polar curve. At points where
the distance is differentiable, its gradient has magnitude one, which gives
`delta = sech(d/epsilon)^2/(2*epsilon)`. Normals are not unique on the
polygon's medial axis. The benchmark fixes this polygon in all comparisons.
`testStarfish` compares 300 deterministic random queries with an exhaustive
segment search and MATLAB `inpolygon`, checks boundary vertices, and solves
an independent exact constant transmission problem.

## Refinement and comparisons

The default adaptive hierarchy has a 64-by-64 base grid and five levels.
At level L it tags `abs(d) < max(3*epsilon,3*h_L)`, then buffers the tags by
two cells. A separate check confirms that all 1024-grid cell centers in
`abs(d)<epsilon` lie in a finest-level patch. This is a sampled coverage
check of that band, not a proof for arbitrary geometries and settings.

At fixed epsilon, the full driver solves uniform 512 and 1024 grids with
the same distance and coefficients. For the main AMR comparison, the
uniform 1024 cell averages are restricted to each visible AMR cell, and
the difference uses cell-area weights. Covered coarse cells are excluded.
The comparison panel displays that same difference on each visible AMR
cell footprint; its block structure reflects the local grid spacing.

These comparisons measure agreement between numerical discretizations of
the fixed diffuse problem. They do not measure error against a known sharp
interface solution or establish an epsilon-convergence rate. The residual
and the integrated balance defect report different algebraic diagnostics.

## Reproduce

```matlab
bsamSetup
demoStarfish('full')
```

Outputs: three PNG/PDF pairs in `figures/`, numerical metrics in
`verification/starfish-full.json`, and the V-cycle history in
`verification/starfish-full-residual.csv`. The quick mode writes separate
filenames and does not replace the full gallery.

The requested relative residual is `1e-10` for both adaptive and uniform
solves. An initial trial with a `1e-11` uniform target reached approximately
`1.6e-11` on the 1024 grid after 100 cycles; the final driver uses the
shared `1e-10` target and asserts that every solve reaches it.
