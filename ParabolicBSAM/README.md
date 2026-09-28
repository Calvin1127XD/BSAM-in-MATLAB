# ParabolicBSAM (EBSAM-P) — time-dependent block-structured adaptive multigrid

Solves the 2D parabolic problem

```
du/dt - div( D(x,y[,t]) grad u ) + C(x,y[,t]) u = f(x,y,t)   on [a,b] x [c,d]
```

by implicit time stepping in which **every step is one mass-conservative
AFAS solve** of the EBSAM elliptic kernel (`../EllipticBSAM`) — the
approach of Feng, Guo, Lowengrub & Wise (JCP 352, 2018, §6), where the
same machinery time-steps Cahn–Hilliard-type systems.

## Method

* **Schemes** — backward Euler, **BDF2** (default; BE bootstrap on the
  first step), Crank–Nicolson (θ = ½ in its
  doubled form, so the implicit operator is always
  `L_phys + (α/Δt) I` with `α = 1, 3/2, 2`).  The shift enters as
  `C_eff = C + α/Δt`. For the supported `D > 0`, `C >= 0`, `Δt > 0`
  problem it removes the constant Neumann nullspace; a convergence rate
  still depends on the mesh, coefficients and smoother.
* **Warm starts** — uⁿ initializes the solve for u^{n+1}; measured cost
  is **3–7 V-cycles per step** (3.5 for a smoothing heat equation, ~7
  for a fast-moving bump) vs ~10 for a cold solve at the same
  tolerance.
* **Moving refinement** — every `RegridEvery` steps the hierarchy is
  rebuilt from the current solution (any of the indicators: undivided
  Laplacian, RTE, or a user tag function) and the state is transferred
  **conservatively**: direct copies where the resolution is unchanged,
  conservative-linear prolongation (children average exactly to the
  parent) where refined, and — via the maintained fill-down state —
  conservative averages where derefined.  The composite integral of u is
  preserved exactly across regrids; combined with the flux-balanced
  interfaces, the heat equation has a discrete mass-balance identity. Algebraic residuals and floating-point
  errors can accumulate over many steps.
* **Time-dependent data** — `D`, `C`, `f`, and boundary data may be
  `@(x,y,t)`; arity is detected and only what changes is re-evaluated
  (coefficients, folded diagonals, coarsest factorization, BC vectors).
* **Composite consistency** — the history term uⁿ/Δt on every level uses
  the post-fill-down state (covered coarse = restriction of fine), so
  the FAS fixed point is the conservative composite discretization at
  each step.

## Layout

```
src/        ebp*.m  (entry point ebpSolve.m; thin layer over eb*)
examples/   demoHeatMovie.m  (rotating source; writes heat_movie.gif/.mp4)
tests/      runParabolicTests.m
output/     figures, movies, results_parabolic.mat (created on run)
TEST_REPORT.md
```

## Quick start

```matlab
run ParabolicBSAM/src/ebpAddpaths.m
pprob = ebpProblem([0 1 0 1], D, C, f, {'dirichlet', g}, u0, exact);
sol = ebpSolve(pprob, 'Scheme', 'bdf2', 'Dt', 1e-3, 'TFinal', 0.5, ...
                'BaseCells', [64 64], 'MaxLevels', 3, ...
                'Indicator', 'ulap', 'TagThreshold', @(L) 1e-3/4^(L-1), ...
                'RegridEvery', 5, 'TagBuffer', 3);
sol.mass      % composite integral of u per step (conservation history)
sol.cycles    % V-cycles per step (warm-started)
sol.err       % error vs exact at TFinal (if provided)
```

**Buffer rule for moving features**: between regrids a feature moves
`speed × RegridEvery × Δt`; choose `TagBuffer` (in cells of each level)
plus the indicator's own footprint to cover that distance, or regrid
more often.

## Verification

`tests/runParabolicTests.m`: temporal orders on a traveling-Gaussian
space–time MMS (BE ≈ 1, BDF2/CN ≈ 2) with time-dependent Dirichlet data
and regridding every 2 steps; heat-equation mass conservation across
~40 regrids including derefinement; mesh tracking of a rotating source;
warm-start efficiency and large-Δt robustness.  See `TEST_REPORT.md`.

## Discrete equations and time-dependent operators

Write `L(t)u = -div(D(t) grad u) + C(t)u`, with the physical boundary
conditions included. On a fixed composite hierarchy, the implemented
steps are

```text
BE:    (L(t[n+1]) + 1/dt) u[n+1] = f(t[n+1]) + u[n]/dt
BDF2:  (L(t[n+1]) + 3/(2dt)) u[n+1]
                                      = f(t[n+1]) + (2u[n] - u[n-1]/2)/dt
CN:    (L(t[n+1]) + 2/dt) u[n+1]
                                      = 2f(t[n+1/2]) + 2u[n]/dt - L(t[n])u[n]
```

CN uses midpoint forcing and endpoint operators, a second-order variant
of the trapezoidal rule for smooth data. Both old and new diffusion
operators include the coarse–fine flux corrections when `MassCorrection`
is enabled. `ebpRhs` evaluates the old operator before `ebpSetTime`
refreshes coefficients and boundary data. Regridding transfers both
BDF2 solution histories conservatively, but it can still introduce local
spatial interpolation error.

The actual constant step is `TFinal/max(1,round(TFinal/Dt))`; it can be
slightly larger or smaller than the requested `Dt`. `Dt` and `TFinal`
must be finite and positive. Inspect `sol.t`, `sol.relres`, and
`sol.cycles`: reaching the V-cycle cap does not prove that the algebraic
solve reached its requested tolerance.

For zero source, zero reaction, and homogeneous Neumann boundaries,
composite flux cancellation gives BE/CN mass preservation and the BDF2
recurrence `3M[n+1] - 4M[n] + M[n-1] = 0`, up to the volume-integrated
algebraic residual. The BE bootstrap makes the BDF2 mass history constant
when those residuals vanish. Dirichlet data, nonzero flux, source, or
reaction change this balance.

## MATLAB R2026a regression audit

From the repository root:

```sh
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('ParabolicBSAM/tests'); runParabolicRegressions; runParabolicTests"
```

The September 10, 2026 run passed the full suite and regressions for
nonautonomous CN order, time-dependent diffusion/boundary data, and CN
composite mass balance. See [the portfolio verification](../verification/README.md)
for a fresh run of the copied code, commands, measured results and coverage. Recorded timing claims above
are examples, not hardware-independent performance guarantees.
