# DiffuseDomainBSAM — diffuse-domain method on the adaptive solver

Solves elliptic boundary-value problems posed on an **irregular domain**
`Omega_1` (circle / ellipse / star-shaped flower) by the diffuse domain
method, using the general solver in `../EllipticBSAM` with the particular
diffuse-domain coefficients `D(x,y)`, `C(x,y)`, `f(x,y)`, and
block-structured adaptive refinement concentrated on the diffuse
interfacial layer.

## Formulation

The boundary `Gamma = boundary(Omega_1)` is described by a level set
`psi` (`psi < 0` inside) and smoothed by the phase field

```
phi(x) = ( 1 - tanh( 3 psi(x) / eps ) ) / 2,
```

as in Feng et al. (2018), sec. 6.4.  Two reformulations on the enclosing
box (homogeneous Neumann outer boundary) are provided by `ddProblem`:

* **Neumann data** `D du/dn = g` on Gamma (paper eq. 6.18; Li, Lowengrub,
  Rätz & Voigt 2009):

  ```
  beta phi u - div( phi D grad u ) = phi f + |grad phi| g
  ```
  i.e. the general solver is called with `Ddd = phi D`, `Cdd = phi beta`,
  `fdd = phi f + |grad phi| g`.

* **Dirichlet data** `u = g` on Gamma (penalty form):

  ```
  -div( phi D grad u ) + phi C u + (mu/eps^3)(1-phi)(u - g) = phi f
  ```
  i.e. `Cdd = phi C + (mu/eps^3)(1-phi)`, `fdd = phi f + (mu/eps^3)(1-phi) g`.

A floor `phi >= PhiFloor` (default `1e-8`) keeps the operator uniformly
elliptic where `tanh` saturates in double precision; see
`docs/ddm_analysis.pdf` for the perturbation bound.

**Refinement near the internal boundary** is geometric and deterministic:
cells with `|psi| < TagWidth * eps` (default `TagWidth = 2.5`, covering
`|phi - 1/2| < 0.4999...`) are tagged on every level, so the patch
hierarchy tracks Gamma and the hierarchy is intended to follow Gamma. Verify the actual band
coverage and local `h/eps`: center-sampled tags, dropped thin patches,
and `MaxLevels` can leave parts of a narrow layer underresolved.

## Layout

```
src/        ddGeometry.m, ddProblem.m, ddSolve.m, ddErrorInside.m
examples/   demoStarNeumann.m, demoDiscDirichlet.m,
            demoPentagonDirichlet.m
tests/      runDdTests.m        (paper Table-6 replication + ablation)
            runDdGeometries.m   (circle/ellipse/pentagon/polar sweep)
            runDdForced.m       (prescribed forcing, Cauchy self-
                                   convergence, normal-derivative jump)
docs/       ddm_analysis.tex/pdf  (asymptotic eps-convergence analysis)
output/     figures + results_*.mat (created on run)
TEST_REPORT.md
```

Geometries (`ddGeometry`): `'circle'`, `'ellipse'`, the paper's
`'star'`, generic `'polar'` curves r = R(θ) > 0 (`'flower3'`,
`'peanut'`), and exact signed-distance regular `'polygon'`s /
`'pentagon'` (corners: see TEST_REPORT section G).

## Quick start

```matlab
run DiffuseDomainBSAM/src/ddAddpaths.m
geom = ddGeometry('circle', 1);
u = @(x,y) (x.^2+y.^2)/4;
data = struct('f', @(x,y) u(x,y)-1, ...
              'g', @(x,y) 0.5*ones(size(x)), 'C', 1, 'exact', u);
sol = ddSolve(geom, 0.2, [-2 2 -2 2], 'neumann', data, ...
              'BaseCells', [32 32], 'MaxLevels', 4);
sol.errInside                 % visible-cell interior error vs sharp u
sol.mass                      % discrete integrated equation defect
sol.stats.stages{end}.finalRel % composite algebraic residual
```

## Verification

`tests/runDdTests.m` replicates the eps -> 0 study of the paper's
sec. 6.4 (star domain, schedule eps = 0.8 ... 0.1 with eps/h fixed):
measured convergence on a specified coupled `(eps,h)` schedule, small
integrated equation defects, and savings in cell count. The September
2026 run found the adaptive star-domain L2 rates near 1.5 on the finest
intervals; lower cell count did not imply lower wall time on this machine.  Disabling the
flux-balance corrections (the "QI" scheme) visibly degrades both
conservation and accuracy, as in the paper's Table 6.

## Mathematical interpretation and accuracy checks

`ddGeometry('ellipse',...)` uses a dimensionless elliptical level set,
and the star/polar geometries use `r-R(theta)`; these are not exact signed
distances. Near a smooth interface the physical layer width is about
`eps/(3*abs(grad psi))`. In contrast, the circle and polygon geometries
use signed distance (the polygon normal is not unique at corners or on
its medial axis). The `DdmRefinementPath` experiments use a different
profile, `tanh(psi/eps)`, so their numerical epsilon is not directly
interchangeable with this package's.

`g` in the internal Neumann formulation means **flux** `D*du/dn`.
`OuterBC` is passed to `ebProblem`, whose Neumann datum means **outward
normal derivative** `du/dn`. Keeping these conventions separate is
necessary when `D` is not one.

For an extension `u_e`, choosing `f = C*u_e-div(D*grad(u_e))` and
`g = D*grad(u_e)·grad(psi)/abs(grad(psi))` gives a pointwise identity for
the unfloored diffuse PDE. It is an exact boundary-value solution only
if `u_e` also satisfies the outer boundary condition. `PhiFloor` modifies
the diffusion/reaction coefficients but not the phase-weighted source,
so floor and outer-boundary errors remain possible.

Separate four errors in a convergence study: diffuse-model error,
geometry/extension error, spatial discretization error, and algebraic
solve error. Refine `h` at fixed `eps` before interpreting an epsilon
sweep; a fit over a few coupled runs is not a theorem. `ddErrorInside`
reports the unnormalized L2 sum weighted by visible-cell area and the
maximum cell-center error for `psi<0`; it does not integrate cut cells
exactly. A small `sol.mass.rel` measures integrated balance, not pointwise
accuracy or conservation of the integral of `u` for every PDE.

The corrected [analysis](docs/ddm_analysis.pdf) distinguishes continuous
coercivity, formal asymptotics and observed rates, and states the
conditions omitted by the older unconditional QIF/mesh-coverage claims.

```sh
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('DiffuseDomainBSAM/tests'); runDdRegressions; runDdTests; runDdGeometries; runDdForced"
```
