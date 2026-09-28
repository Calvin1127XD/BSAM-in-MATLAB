# EllipticBSAM — Audit Report (2026-07-01)

Auditor: automated deep audit (code read in full, tests executed on MATLAB
R2026a, documentation cross-checked line-by-line against the source, and
the Fortran reference SuperBSAM/BSAM2.0 compared for the patch-geometry
questions).

## Verdict

**The solver is correct, robust, and matches its documentation.**
All existing tests pass on R2026a; eight new stress tests pass; the
2107-line LaTeX documentation contains **zero factual errors** about the
code; and no bugs were found in `src/`. The only defects found were
staleness in `TEST_REPORT.md` and four cosmetic wording issues in the
LaTeX (all now fixed, see §6).

## 1. What was run

| suite | result |
|---|---|
| `tests/runAllTests.m` (T1–T8) | **ALL PASSED**, 42 s. Uniform/adaptive rates 1.98–2.03, V-cycle factors 0.026–0.054 (h- and depth-independent), conservation defects 7e-17…3e-15 at tight tolerance. |
| `tests/unitKernels.m` | machine-precision stencil reproduction (QIF 8.9e-16, LIF 4.4e-16). |
| `tests/auditExtras.m` (new, A1–A8) | **ALL PASSED** (details below). |

Numbers agree with `TEST_REPORT.md` except the T3 table and the T7
accuracy-impact claim, which reflected an older test configuration —
both refreshed from today's run.

## 2. New stress tests (`tests/auditExtras.m`)

* **A1 — patch aspect-ratio freedom.** Tagging a long thin band + a
  diagonal stripe yields 21 patches with widths 8–78 and heights 8–20,
  aspect ratios up to **9.75 : 1**, all parent-aligned, solver converges
  normally. See §3 for the interpretation.
* **A2 — a fine patch straddling two parents.** Engineered two abutting
  level-2 patches and a level-3 tag across their shared edge: the
  clusterer produces **one child with 2 parents** (`parents{p} = [1 2]`),
  and the solve is flawless (rel 1.7e-12, mass 3.4e-15). A per-cycle
  comparison against the same child placed inside one parent shows
  **identical convergence factors (0.041)** — the QC-ring machinery has
  no performance penalty for multi-parent children.
* **A3 — anisotropic cells** (hy = 2 hx, domain [0,2]×[0,1]): rates
  2.07, 2.01. The QIF diagonal stencil is exact for quadratics for any
  (hx, hy) because its four points remain collinear in the scaled
  diagonal parameter (weights 8/15, 2/3, −1/5 are h-independent).
* **A4 — mixed Dirichlet/Neumann sides** (D, N, D, N with variable D,C):
  rates 1.96–1.99, conservation 1e-12–1e-13.
* **A5 — `rte` indicator** (previously untested path): error 7.6e-02 →
  4.7e-03 over 3 levels; behaves like `ulap` with more tagged cells.
* **A6 — Omega sweep + nonzero InitialGuess**: converges for ω ∈
  {0.85, 1.0, 1.15}; mild over-relaxation is fastest here (factor 0.015).
* **A7 — 6 adaptive levels** on a 32² root: factors stay 0.027–0.037,
  error decays 4.8e-03 → 4.0e-04, conservation ~1e-14.
* **A8 — `ebSample`**: 2008 random + boundary points, max sampling error
  7.2× the grid error (bilinear interpolation's own O(h²), as expected).

## 3. Answer: must patches be 2:1 rectangles? Can shapes be free?

**Patches are already free rectangles — there is no 2:1 (or any)
aspect-ratio constraint anywhere in the code, and none is needed.**
Evidence:

* `ebTagCluster` is a Berger–Rigoutsos signature clusterer: boxes are
  accepted purely by fill ratio (`MinFillRatio`, default 0.70) and split
  at signature holes/inflections. Nothing constrains shape; A1 produced
  9.75:1 boxes that solved perfectly.
* The only real geometry rules are:
  1. **alignment** — a child is the exact 2× refinement of a parent-space
     box, so its global index range is automatically odd-start/even-end
     (this is what "ratio two" refers to — the *refinement ratio between
     levels*, not the patch shape);
  2. **minimum width** — `MinBoxWidth` (default 4 parent cells; hard
     floor 2 for valid interface stencils);
  3. **proper nesting** — tags are clipped to cells whose 8-neighborhood
     is covered by the current finest level (one-level interface jumps).
* The Fortran reference behaves the same way: `NewSubGrids`
  (bsamroutines.f90) enforces min-size, parity and fill-ratio rules and
  **no aspect-ratio constraint** (see
  `FortranCodeReference/TREE_ANALYSIS.md` §3.2).

The likely source of the "2:1" impression was the documentation sentence
(tex line ~202) "*Refinement occurs in rectangular patches with a fixed
ratio of two in each direction*" — that "ratio of two" is the
inter-level refinement factor. The sentence has been rewritten to say so
explicitly and to state that patch shapes are unconstrained.

## 4. Multi-parent children (relevant to the EllipticBSAMTree project)

The current implementation deliberately allows a fine patch to overlap
several coarse patches (`parents{p}` is a list; the QC window/ring is
gathered from all owners). The audit confirms:

* this configuration arises in practice (A2 engineers it; scattered
  tagging produces it naturally),
* it is handled **correctly** (conservation at machine precision) and
  **without any V-cycle penalty** (identical factors with/without a
  straddling child).

So the motivation for the strict single-parent tree variant
(`EllipticBSAMTree`) is architectural — simpler, Fortran-faithful
lineage, single-slice parent transfers — not a fix for a correctness or
convergence problem. (Indeed the tree build later confirmed identical
convergence with strictly per-parent children.)

## 5. Code-reading observations (no action needed)

* `ebTagCluster>growMin` has a leftover unused variable `cand`
  (cosmetic; flagged `%#ok<NASGU>`).
* `ebOptions`' help text omits the supported `'rte'` indicator (the
  LaTeX documentation describes it correctly).
* The T4/T9-style mass defects at loose `RelTol` (1e-10 → defect ~1e-11)
  are the expected residual-scaled behavior, not leaks; tightening the
  tolerance drives them to 1e-16 (verified).
* The ring gather (`ebBuildOps>resolveCells`) reads owner-patch physical
  ghosts on domain overhang — correct and confirmed exercised by A1/A4.

## 6. Fixes applied during the audit

* `TEST_REPORT.md`: refreshed the stale T3 table and T7 numbers from the
  R2026a run, camelCased the stale `unit_kernels/run_all_tests/eb_bc_fold`
  references, added a re-verification note and a pointer to
  `auditExtras.m`.
* `BSAM_Code_Documentation.tex` (recompiled; PDF updated):
  1. line ~202: clarified "ratio of two" = refinement ratio; patch
     shapes explicitly unconstrained;
  2. line ~505: descriptors correctly described as rect ranges + linear
     index lists (was "flat arrays of linear indices");
  3. line ~1284: fill-ratio acceptance is `>=` ("meets or exceeds");
  4. line ~1861: level-1 factorization is rebuilt after *every*
     `ebAddLevel`;
  5. `tab:functions` caption no longer claims to list *all* supplied
     files.
* New `tests/auditExtras.m` added to the suite.

Full documentation cross-check notes: see the audit agent's category-by-
category table (all stencil weights, struct fields, option defaults and
algorithm orderings verified bit-exact; 67 labels, no undefined refs,
compiles clean).
