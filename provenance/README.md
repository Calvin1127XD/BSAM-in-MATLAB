# Source provenance

## Solver snapshot

The four packages were copied from the working tree of
[`Calvin1127XD/BSAMElliptic`](https://github.com/Calvin1127XD/BSAMElliptic)
on September 28, 2026. That working tree included local corrections beyond
its base commit. [source-manifest.json](source-manifest.json) records the
base commit and SHA-256 of each retained file at the time of copying.
The original research checkout was left in place.

All retained MATLAB solver, example, and test files match that snapshot.
The new starfish adapter, geometry, driver, plotting code, regression test,
setup function, and suite runner are additions. Two package READMEs were
adjusted to link to this repository's verification and to remove a link to
the omitted optional reinforcement-learning experiments.

The manifest also lists files excluded after the initial copy: unrelated
editor drafts, optional reinforcement-learning code, a historical generated
convergence report, and a dissertation-only script with external paths.
Caches, compiled binaries, generated test output, and other solver families
were excluded during copying. Existing PDF documentation and historical
test reports are included unchanged. The old audit reports can refer to the
larger research checkout; freshly executed results are in [verification](../verification/).

## Starfish reference

Reference repository:
[`Calvin1127XD/DiffuseDomainHelmholtz`](https://github.com/Calvin1127XD/DiffuseDomainHelmholtz).

Pinned commit: `a9d8449005f4ab030e64cf4c346c7225998fa01e`.

| Reference file | Reused specification |
|---|---|
| `DiffuseDomainMultigrid/2DCC/generateBoundary.m` | Polar curve and 3000-point boundary |
| `DiffuseDomainMultigrid/2DCC/elliptic2DMGCCDiff.m` | Box, phase profile, coefficients, forcing, outer Neumann condition |
| `README.md` | Documentation structure and linked PNG/PDF figures |
| `LICENSE` | Apache 2.0 license text |

The solver here is BSAM. Its distance evaluator computes distances directly
to boundary segments; it does not copy the reference's Fast Marching solver.
The reference's leftover ellipse error field is not used. See
[starfish.md](../docs/starfish.md) for the exact comparison scope.

## Check the package

```sh
python3 scripts/verify_package.py
```

This checks source hashes for retained MATLAB files, local Markdown links,
figure sizes, and the recorded benchmark's acceptance conditions. It does
not rerun the PDE solvers; use `runVerification` and `demoStarfish('full')`
in MATLAB for numerical reproduction.
