# EBSAM (EllipticBSAM) — mass-conservative block-structured adaptive FAS multigrid

A clean, vectorized MATLAB implementation of the cell-centered,
mass-conservative adaptive FAS multigrid method of

> W. Feng, Z. Guo, J. S. Lowengrub, S. M. Wise,
> *A mass-conservative adaptive FAS multigrid solver for cell-centered
> finite difference methods on block-structured, locally-cartesian grids*,
> J. Comput. Phys. **352** (2018) 463–497,

generalized to the variable-coefficient problem

```
-div( D(x,y) grad u ) + C(x,y) u = f(x,y)   on [a,b] x [c,d],
```

with `D > 0`, `C >= 0`, and Dirichlet and/or Neumann data on each side of
the box.

## Method summary

* **Discretization** — cell-centered finite differences in flux form
  (paper eqs. 4.4–4.5); `D` evaluated at face midpoints; quadratic ghost
  elimination for Dirichlet sides, natural fluxes for Neumann sides.
* **Grid hierarchy** — one ladder of levels with refinement ratio 2:
  uniform full-domain grids below the root (classic multigrid coarsening,
  direct solve at the coarsest), the root grid, and adaptive levels above
  it made of rectangular patches with an explicit parent/child tree
  (`children`/`parents` per patch), 1-conforming interfaces and proper
  nesting enforced at regridding.
* **Coarse–fine coupling** — fine-grid ghost cells by **QIF** quadratic
  interpolation along the diagonal Λ stencil with weights
  `(8/15, 2/3, −1/5)` (paper eq. 3.12 / Fig. 4a), tangential end-cell
  weights `(5/32, 15/16, −3/32)`; **LIF** linear variant available.
* **Mass conservation** — the coarse-level force is corrected at every
  cell adjacent to a coarse–fine interface by the flux-balance ("zombie")
  correction `C = D_c (Z − R u)/h_c²` (paper eqs. 4.10/4.20), which makes
  the coarse equations see the fine fluxes exactly.  The discrete
  conservation identity holds to roundoff (see `ebMassCheck`).
* **Solver** — adaptive FAS (AFAS) V-cycles (paper alg. a0–a9 with
  recursion), red–black Gauss–Seidel smoothing with *global* cell parity,
  conservative 4-point restriction, bilinear correction prolongation,
  exact folding of the physical boundary conditions into the smoother
  (boundary cells perform exact Gauss–Seidel of the ghost-eliminated
  equations).
* **Regridding** — undivided-Laplacian indicator (paper eq. 4.32), the
  relative truncation error `Indicator='rte'`
  (τ = L₂ₕ(Ru) − R(Lₕu), paper §4.4), or a user tag function; buffered
  tags, Berger–Rigoutsos signature clustering, automatic enforcement of
  even alignment, minimum patch width, proper nesting; sliver clusters
  are merged or dropped rather than ever producing fragile patches (the
  failure mode of earlier attempts).

## Layout

```
src/        solver source (eb_*.m), entry point ebSolve.m
examples/   demoBasicUniform.m, demoPoissonAdaptive.m,
            demoVisualWalkthrough.m  <- animated step-by-step tour of
                                          the solver (slideshow in
                                          desktop MATLAB; also writes
                                          output/walkthrough/*.png,
                                          walkthrough.gif, walkthrough.mp4)
tests/      unitKernels.m, smokeAmr.m, runAllTests.m, mmsCases.m
docs/       convergence_theory.tex/pdf  (discrete theory + proofs)
output/     figures and results.mat written by tests (created on run)
TEST_REPORT.md
```

## Quick start

```matlab
addpath EllipticBSAM/src
uex = @(x,y) sin(pi*x).*cos(pi*y);                 % manufactured solution
D = @(x,y) 1 + 0.5*sin(2*pi*x).*cos(pi*y);
C = @(x,y) 1 + x.*y;
f = ...;                                           % see examples/
prob = ebProblem([0 1 0 1], D, C, f, {'dirichlet', uex}, uex);
sol  = ebSolve(prob, 'BaseCells', [64 64], 'MaxLevels', 3, ...
                'Indicator', 'ulap', 'TagThreshold', @(L) 1e-2/4^(L-1));
sol.err          % composite error norms vs exact
sol.mass         % discrete conservation check
ebPlot(sol.H, 'mesh', 'mesh.png');
```

Boundary conditions per side:
`bc.west = {'dirichlet', g}` or `{'neumann', g}` with `g(x,y)` the
boundary value resp. the outward normal derivative; a single
`{type, g}` applies to all sides.  The pure-Neumann, `C == 0` singular
case is handled (augmented coarsest solve + mean pinning).

## Verification

Run `tests/runAllTests.m` (see `TEST_REPORT.md` for current numbers):
second-order convergence on uniform, fixed multi-level, L-shaped, and
fully adaptive hierarchies; h-independent V-cycle factors ~0.01–0.05;
conservation to ~1e-16; robustness under adversarial scattered tagging.

The [implementation audit appendix](docs/implementation_audit.tex), also
included in `BSAM_Code_Documentation.pdf`, explains the corrected Neumann
gauge and source-projection conventions and an explicit multidimensional
QIF counterexample to the companion theory's matrix-sign claim. The
theory's dependent proof claims remain review items; numerical convergence
results apply to the tested families. `tests/testEbNeumannSemantics.m`
checks these diagnostic and compatibility regressions.

For singular Neumann problems, inspect `sol.H.compatibility.sourceShift`:
the solver subtracts this constant from the originally sampled source to
enforce discrete compatibility after every refinement. Error means are
removed only when `C==0` and all boundaries are Neumann, for every source
quadrature; positive reaction uses absolute error. Historical tables may
have used a different mean-alignment convention.
