%auditExtras Supplementary audit tests for EllipticBSAM (2026-07-01 audit).
%   Exercises paths not covered by runAllTests:
%     A1  patch aspect-ratio freedom (are boxes forced to be 2:1? -> no)
%     A2  a fine patch straddling TWO parent patches (multi-parent lineage)
%     A3  anisotropic cells (hx ~= hy) + rectangular domain convergence
%     A4  mixed Dirichlet/Neumann sides convergence
%     A5  'rte' (relative truncation error) indicator path
%     A6  under/over-relaxation and nonzero InitialGuess smoke
%     A7  deep hierarchy (6 adaptive levels) error decay
%     A8  ebSample point sampling accuracy
%   Prints PASS/FAIL per test and errors out at the end on any failure.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
PASS = true;
say = @(varargin) fprintf(varargin{:});

% ======================================================================
say('\n=== A1: patch aspect-ratio freedom ===============================\n');
% Tag a long thin horizontal band and a diagonal stripe. If the clusterer
% imposed any fixed shape (e.g. 2:1), the resulting boxes could not track
% these features. We record every level-2 box's width:height ratio.
uex = @(X, Y) sin(pi*X) .* cos(pi*Y);
f   = @(X, Y) (2*pi^2 + 1) * uex(X, Y);
prob = ebProblem([0 1 0 1], 1, 1, f, {'dirichlet', uex}, uex);
bandfn = @(X, Y, L) (abs(Y - 0.70) < 0.035 & X > 0.06 & X < 0.94) | ...
                    (abs(Y - X) < 0.05);
opts = ebOptions('BaseCells', [64 64], 'MaxLevels', 2, 'Indicator', 'fun', ...
  'TagFunction', bandfn, 'TagBuffer', 0, 'Verbose', 0, 'RelTol', 1e-9);
sol = ebSolve(prob, opts);
bx = sol.levels(2).box;                    % level-2 boxes, fine index space
w = bx(:,2) - bx(:,1) + 1;
h = bx(:,4) - bx(:,3) + 1;
ar = max(w./h, h./w);
say('  %d level-2 patches; widths %d..%d, heights %d..%d\n', ...
  size(bx,1), min(w), max(w), min(h), max(h));
say('  aspect ratios: min %.2f, median %.2f, max %.2f\n', ...
  min(ar), median(ar), max(ar));
% even alignment to parent cells must hold for every box
alignOk = all(mod(bx(:,1),2)==1) && all(mod(bx(:,2),2)==0) && ...
          all(mod(bx(:,3),2)==1) && all(mod(bx(:,4),2)==0);
say('  parent-cell alignment (ilo odd, ihi even, ...): %d\n', alignOk);
if ~alignOk, say('  ** A1 FAIL (alignment)\n'); PASS = false; end
if max(ar) < 4                              % the band should be much longer
  say('  ** A1 FAIL (expected skinny boxes; clusterer is shape-constrained?)\n');
  PASS = false;
else
  say('  => boxes are free rectangles (no 2:1 or square constraint)\n');
end
if sol.stats.stages{end}.finalRel > 1e-8
  say('  ** A1 FAIL (solver on skinny patches)\n'); PASS = false;
end

% ======================================================================
say('\n=== A2: fine patch straddling two parent patches =================\n');
% Engineer two adjacent level-2 patches sharing the edge x = 0.5: tag a
% tall box left of x = 0.5 and a shorter box right of it. With a high
% MinFillRatio the clusterer must split at the signature inflection
% (x = 0.5), giving exactly two side-by-side patches. A level-3 tag box
% centered ON that edge then yields one child intersecting BOTH: the
% current (QC-ring) scheme permits multi-parent children by design.
twoBoxes = @(X, Y) (X >= 0.25 & X < 0.5 & Y >= 0.25 & Y < 0.75) | ...
                   (X >= 0.5  & X < 0.75 & Y >= 0.40 & Y < 0.60);
straddle = @(X, Y) abs(X - 0.5) < 0.04 & abs(Y - 0.5) < 0.03;
tagfn = @(X, Y, L) (L == 1) .* twoBoxes(X, Y) + (L == 2) .* straddle(X, Y);
cs = mmsCases('gauss52');
% RelTol 1e-11 (not 1e-10): the 1e-11 mass-defect assert below is only
% consistent with the residual level actually reached -- the defect scales
% with the converged residual (2.3e-11 at RelTol 1e-10, 9.5e-13 at 1e-11,
% 3.4e-15 at the roundoff floor), matching the tree suite's T9 twin which
% also solves this scenario at RelTol 1e-11.
opts = ebOptions('BaseCells', [64 64], 'MaxLevels', 3, 'Indicator', 'fun', ...
  'TagFunction', tagfn, 'TagBuffer', 0, 'MinFillRatio', 0.95, ...
  'Verbose', 0, 'RelTol', 1e-11, 'MaxVCycles', 40);
sol = ebSolve(cs.prob, opts);
H = sol.H;
kTop = H.nlev;                             % the level-3 AMR level
maxParents = 0;
for p = 1:H.lev{kTop}.np
  maxParents = max(maxParents, numel(H.lev{kTop}.parents{p}));
end
say('  level-2 patches: %d,  level-3 patches: %d,  max #parents of a level-3 patch: %d\n', ...
  H.lev{kTop-1}.np, H.lev{kTop}.np, maxParents);
say('  solver: V=%d rel = %.2e, factor = %.3f, mass = %.1e\n', ...
  sol.stats.stages{end}.vcycles, sol.stats.stages{end}.finalRel, ...
  sol.stats.stages{end}.avgFactor, sol.mass.rel);
if maxParents < 2
  say('  ** A2 FAIL (no straddling child; expected a multi-parent patch)\n');
  PASS = false;
end
if sol.stats.stages{end}.finalRel > 5e-10 || sol.mass.rel > 1e-11 ...
    || sol.stats.stages{end}.avgFactor > 0.2
  say('  ** A2 FAIL (straddling-child solve degraded)\n'); PASS = false;
end

% ======================================================================
say('\n=== A3: anisotropic cells + rectangular domain ====================\n');
% Domain [0,2]x[0,1] with BaseCells [m, m/2]/... chosen so hy = 2*hx (true
% cell anisotropy). QIF's diagonal stencil is exact for quadratics for ANY
% (hx, hy) (the four stencil points stay collinear in the scaled diagonal
% parameter), so second order must survive.
uex = @(X, Y) sin(pi*X) .* cos(pi*Y);
D = @(X, Y) 1 + 0.5*X;
C = @(X, Y) 1 + 0*X;
% f = -(D_x u_x + D (u_xx + u_yy)) + C u  with D_x = 1/2:
f = @(X, Y) -0.5*pi*cos(pi*X).*cos(pi*Y) ...
          + 2*pi^2*(1 + 0.5*X).*sin(pi*X).*cos(pi*Y) ...
          + uex(X, Y);
prob = ebProblem([0 2 0 1], D, C, f, {'dirichlet', uex}, uex);
centerfn = @(X, Y, L) (X - 1).^2 + (Y - 0.5).^2 < 0.2^2;
errs = [];
for base = {[64 16], [128 32], [256 64]}
  opts = ebOptions('BaseCells', base{1}, 'MaxLevels', 2, 'Indicator', 'fun', ...
    'TagFunction', centerfn, 'Verbose', 0, 'RelTol', 1e-10);
  sol = ebSolve(prob, opts);
  errs(end+1) = sol.err.linf; %#ok<*SAGROW>
  say('  base=[%3d %2d] (hx=%.4f hy=%.4f): err_inf=%.4e factor=%.3f mass=%.1e\n', ...
    base{1}(1), base{1}(2), 2/base{1}(1), 1/base{1}(2), ...
    sol.err.linf, sol.stats.stages{end}.avgFactor, sol.mass.rel);
end
r3 = log2(errs(1:end-1) ./ errs(2:end));
say('  L_inf rates (hy = 2*hx anisotropy): %s\n', mat2str(r3, 3));
if any(r3 < 1.8), say('  ** A3 FAIL\n'); PASS = false; end

% ======================================================================
say('\n=== A4: mixed Dirichlet/Neumann sides =============================\n');
uex = @(X, Y) exp(X) .* sin(pi*Y) + 0.5*X.^2;
D = @(X, Y) 1 + 0.5*X.*Y;
C = @(X, Y) 1 + X + Y;
% f = -div(D grad u) + C u, with ux = e^x sin(pi y) + x, uy = pi e^x cos(pi y)
%   Dx = 0.5 Y, Dy = 0.5 X, uxx = e^x sin(pi y) + 1, uyy = -pi^2 e^x sin(pi y)
f = @(X, Y) -( 0.5*Y .* (exp(X).*sin(pi*Y) + X) ...
             + (1 + 0.5*X.*Y) .* (exp(X).*sin(pi*Y) + 1) ...
             + 0.5*X .* (pi*exp(X).*cos(pi*Y)) ...
             + (1 + 0.5*X.*Y) .* (-pi^2*exp(X).*sin(pi*Y)) ) ...
           + C(X, Y) .* uex(X, Y);
bc.west  = {'dirichlet', @(X, Y) uex(X, Y)};
bc.east  = {'neumann',   @(X, Y) exp(X).*sin(pi*Y) + X};   % du/dn = +du/dx
bc.south = {'dirichlet', @(X, Y) uex(X, Y)};
bc.north = {'neumann',   @(X, Y) pi*exp(X).*cos(pi*Y)};    % du/dn = +du/dy
prob = ebProblem([0 1 0 1], D, C, f, bc, uex);
cornerfn = @(X, Y, L) (X - 0.75).^2 + (Y - 0.75).^2 < 0.15^2;
errs = [];
for m = [32 64 128 256]
  opts = ebOptions('BaseCells', [m m], 'MaxLevels', 2, 'Indicator', 'fun', ...
    'TagFunction', cornerfn, 'Verbose', 0, 'RelTol', 1e-10);
  sol = ebSolve(prob, opts);
  errs(end+1) = sol.err.linf;
  say('  m=%4d: err_inf=%.4e factor=%.3f mass=%.1e\n', m, sol.err.linf, ...
    sol.stats.stages{end}.avgFactor, sol.mass.rel);
end
r4 = log2(errs(1:end-1) ./ errs(2:end));
say('  L_inf rates: %s\n', mat2str(r4, 3));
if any(r4 < 1.8), say('  ** A4 FAIL\n'); PASS = false; end

% ======================================================================
say('\n=== A5: rte indicator path ========================================\n');
cs = mmsCases('threebumps');
errsR = [];
for nlv = 1:3
  opts = ebOptions('BaseCells', [64 64], 'MaxLevels', nlv, 'Indicator', 'rte', ...
    'TagThreshold', @(L) 2e-3/4^(L-1), 'TagBuffer', 2, 'Verbose', 0, ...
    'RelTol', 1e-9, 'MaxVCycles', 40);
  sol = ebSolve(cs.prob, opts);
  errsR(end+1) = sol.err.linf;
  say('  MaxLevels=%d: err_inf=%.4e dof=%d factor=%.3f\n', nlv, ...
    sol.err.linf, sol.err.dof, sol.stats.stages{end}.avgFactor);
end
if errsR(3) > 0.35 * errsR(1)
  say('  ** A5 FAIL (rte tagging did not reduce error)\n'); PASS = false;
end

% ======================================================================
say('\n=== A6: Omega and InitialGuess smoke ==============================\n');
cs = mmsCases('vardc');
for om = [0.85 1.0 1.15]
  opts = ebOptions('BaseCells', [64 64], 'Omega', om, 'Verbose', 0, ...
    'RelTol', 1e-10, 'InitialGuess', @(X, Y) cos(3*X) .* Y);
  sol = ebSolve(cs.prob, opts);
  say('  Omega=%.2f: V=%d factor=%.3f err=%.3e\n', om, ...
    sol.stats.stages{end}.vcycles, sol.stats.stages{end}.avgFactor, sol.err.linf);
  if sol.stats.stages{end}.finalRel > 1e-9
    say('  ** A6 FAIL (Omega=%.2f did not converge)\n', om); PASS = false;
  end
end

% ======================================================================
say('\n=== A7: deep hierarchy (6 adaptive levels) ========================\n');
cs = mmsCases('threebumps');
errs7 = [];
for nlv = [4 5 6]
  opts = ebOptions('BaseCells', [32 32], 'MaxLevels', nlv, ...
    'Indicator', 'ulap', 'TagThreshold', @(L) 4e-2/4^(L-1), 'TagBuffer', 2, ...
    'Verbose', 0, 'RelTol', 1e-9, 'MaxVCycles', 60);
  sol = ebSolve(cs.prob, opts);
  errs7(end+1) = sol.err.linf;
  say('  MaxLevels=%d (got %d): err_inf=%.4e dof=%d factor=%.3f mass=%.1e\n', ...
    nlv, numel(sol.levels), sol.err.linf, sol.err.dof, ...
    sol.stats.stages{end}.avgFactor, sol.mass.rel);
  if sol.stats.stages{end}.avgFactor > 0.2
    say('  ** A7 FAIL (V-cycle degraded at depth)\n'); PASS = false;
  end
end
if errs7(end) > 0.15 * errs7(1)
  say('  ** A7 FAIL (deep refinement not reducing error)\n'); PASS = false;
end

% ======================================================================
say('\n=== A8: ebSample accuracy =========================================\n');
cs = mmsCases('vardc');
opts = ebOptions('BaseCells', [64 64], 'MaxLevels', 3, 'Indicator', 'ulap', ...
  'TagThreshold', @(L) 1e-2/4^(L-1), 'Verbose', 0, 'RelTol', 1e-10);
sol = ebSolve(cs.prob, opts);
rng(7);
xq = rand(2000, 1); yq = rand(2000, 1);
% include exact domain corners and edge midpoints
xq = [xq; 0; 1; 0; 1; 0.5; 0.5; 0; 1];
yq = [yq; 0; 0; 1; 1; 0; 1; 0.5; 0.5];
v = ebSample(sol.H, xq, yq);
ve = cs.prob.exact(xq, yq);
sampErr = max(abs(v - ve));
say('  max |sample - exact| over %d points = %.3e (grid err_inf = %.3e)\n', ...
  numel(xq), sampErr, sol.err.linf);
% bilinear sampling adds its own O(h^2); allow a modest multiple
if sampErr > 50 * sol.err.linf
  say('  ** A8 FAIL (sampling inconsistent with grid error)\n'); PASS = false;
end

% ======================================================================
say('\n');
if PASS
  say('ALL AUDIT-EXTRA TESTS PASSED\n');
else
  error('auditExtras: at least one test failed');
end
