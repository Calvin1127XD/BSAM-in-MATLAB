# ParabolicBSAM — test report

All runs: MATLAB R2025b headless (`-batch`), macOS. Suite:
`tests/run_parabolic_tests.m` (27 s total). Saved results:
`output/results_parabolic.mat`.

The solver advances `du/dt − ∇·(D∇u) + Cu = f` by implicit stepping in
which each step is one warm-started AFAS solve of the EBSAM elliptic
kernel with diagonal shift `C_eff = C + α/Δt` (α = 1 BE, 3/2 BDF2,
2 CN-doubled).

## P1a — temporal orders (Δt-Cauchy on a fixed uniform 128² mesh)

Traveling-Gaussian space–time MMS, σ = 0.06, center moving at
(0.4, 0.3), T = 0.4, time-dependent Dirichlet data. Comparing final
fields between consecutive Δt on the *same* mesh cancels the spatial
error exactly, isolating the temporal order:

| scheme | Cauchy diffs (Δt = 0.04, /2, /4, /8) | rates |
|---|---|---|
| BE   | 1.23e-3, 6.19e-4, 3.10e-4 | **0.99, 1.00** |
| BDF2 | 3.57e-4, 9.02e-5, 2.26e-5 | **1.98, 2.00** |
| CN   | 1.04e-2, 2.64e-3, 6.54e-4 | **1.97, 2.01** |

(CN's larger constant is its weak damping of the bump's stiff modes;
the order is clean.)

## P1b — joint (h, Δt) refinement with the moving adaptive mesh

Same MMS; BDF2; root 32²/64²/128² with `MaxLevels` 2, Δt halved with h,
**regridding every 2 steps** (15 regrids at the finest), conservative
state transfer at every regrid:

| root | finest h | Δt | ‖e‖∞ at T | rate |
|---|---|---|---|---|
| 32²  | 1/64  | 0.05   | 8.19e-3 | — |
| 64²  | 1/128 | 0.025  | 2.30e-3 | 1.83 |
| 128² | 1/256 | 0.0125 | 6.33e-4 | 1.86 |

Combined 2nd order (O(Δt²) + O(h²)) survives the moving mesh and the
repeated transfers. Two lessons surfaced while building this test:

* the `ulap` indicator scales like h²·|Δu|, so a refinement sweep must
  scale `TagThreshold` by 4⁻ʲ — otherwise the finest run silently stops
  refining and the apparent rate collapses to ~1;
* the BDF2 history u^{n−1} must be transferred at 2nd order
  (conservative-linear) — piecewise-constant prolongation injects
  O(h·|∇u|) noise exactly where fresh cells appear (the feature's
  leading edge).

## P2 — heat-equation mass conservation across regrids

`du/dt = ∇·(D∇u)`, D = 1 + 0.5 sin(2πx)cos(2πy), homogeneous Neumann,
f = 0; two Gaussian bumps over a 0.5 background. BDF2, Δt = 2e-3,
200 steps, root 64², `MaxLevels` 3, regrid every 5 steps (39 regrids).

* mass drift over the full run: **max|m(tₙ) − m(t₀)| = 1.09e-10**
  (relative ~2e-10; tracks the solver tolerance 1e-9, not the step or
  regrid count);
* the hierarchy **derefines** as the solution smooths — 3 levels → 1
  level by the end — so the conservative *coarsening* branch of the
  transfer (fill-down averages) is exercised, not just refinement;
* 3.5 V-cycles/step mean.

This is the discrete counterpart of the conservation identity: flux
form + flux-balanced coarse–fine interfaces (zombie corrections) make
each step conservative; the regrid transfer (same-level copies,
conservative-linear prolongation, fill-down averages) preserves the
composite integral exactly, so the only drift is algebraic.

## P3 — moving-source tracking

Rotating Gaussian source (radius 0.22, period 1, σ = 0.03), BDF2,
Δt = 5e-3, root 64², `MaxLevels` 3, regrid every 4 steps (49 regrids).

* final source position (0.720, 0.500) is covered by the finest level;
* composite size stays ≈ **5.7k cells vs 65.5k** for the uniform 256²
  equivalent (11.5×) while the nest orbits with the source;
* figures: `output/p3_tracking_mesh.png`, `p3_tracking_solution.png`;
  movie of the same physics: `output/heat_movie.gif/.mp4`
  (`examples/demo_heat_movie.m`).

## P4 — warm starts and robustness

* mean V-cycles/step (P2): **3.51** (cold elliptic solves at the same
  tolerance need ~10); the moving-bump problems run at 5.8–7 cycles
  per step — the faster the solution changes, the more cycles a step
  needs, but never close to a cold start;
* backward Euler with Δt = 0.1 (diffusive CFL ≈ 1e-5 at this h, i.e.
  Δt/CFL ≈ 10⁴): bounded, monotone, max|u| = 0.513 — the implicit
  solver is unconditionally stable and the multigrid does not care
  about the size of α/Δt.

## Verdict

| test | check | result |
|---|---|---|
| P1a | BE/BDF2/CN temporal orders 1/2/2 | **pass** (0.99 / 2.00 / 2.01) |
| P1b | joint order ≥ 1.5 on moving mesh | **pass** (1.86) |
| P2  | relative mass drift < 1e-9 | **pass** (2e-10) |
| P3  | finest level covers source at T | **pass** |
| P4  | ≤ 4.5 cycles/step; large-Δt bounded | **pass** (3.51; 0.513) |
