%testNeumannCompat Multi-level SINGULAR pure-Neumann convergence -- the
%   regression test for the composite-quadrature compatibility stall
%   (STRESS_REPORT.md section 4, fixed by ebCompatProject /
%   ebtCompatProject).
%
%   Before the fix, box-in-box hierarchies on the all-Neumann C = 0
%   problem with the default midpoint source quadrature (SourceQuad = 1)
%   stalled at an O(h^2) incompatibility floor (factor -> 1.000); only
%   SourceQuad >= 3 cured it.  With the source projection the default
%   quadrature must converge at the usual ~0.05 factor and second order,
%   for BOTH solvers.
%
%   MMS: u = cos(2 pi x) cos(2 pi y) on [0,1]^2, D = 1, C = 0,
%   homogeneous Neumann on all sides (du/dn = 0 exactly), f = 8 pi^2 u.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
addpath(fullfile(here, '..', '..', 'EllipticBSAMTree', 'src'));
PASS = true;
say = @(varargin) fprintf(varargin{:});

uex = @(X, Y) cos(2*pi*X) .* cos(2*pi*Y);
f   = @(X, Y) 8*pi^2 * uex(X, Y);
tagfn = @(X, Y, L) ...
  (L == 1) .* (max(abs(X - 0.5), abs(Y - 0.5)) < 0.30) + ...
  (L == 2) .* (max(abs(X - 0.5), abs(Y - 0.5)) < 0.15);

for solver = 1:2
  if solver == 1
    name = 'ring';
    prob = ebProblem([0 1 0 1], 1, 0, f, {'neumann', 0}, uex);
    run1 = @(m, L) ebSolve(prob, 'BaseCells', [m m], 'MaxLevels', L, ...
      'Indicator', 'fun', 'TagFunction', tagfn, 'TagBuffer', 0, ...
      'Verbose', 0, 'RelTol', 1e-10, 'MaxVCycles', 40, 'SourceQuad', 1);
  else
    name = 'tree';
    prob = ebtProblem([0 1 0 1], 1, 0, f, {'neumann', 0}, uex);
    run1 = @(m, L) ebtSolve(prob, 'BaseCells', [m m], 'MaxLevels', L, ...
      'Indicator', 'fun', 'TagFunction', tagfn, 'TagBuffer', 0, ...
      'Verbose', 0, 'RelTol', 1e-10, 'MaxVCycles', 40, 'SourceQuad', 1);
  end
  say('\n=== %s: 3-level singular pure-Neumann, SourceQuad = 1 ===\n', name);
  errs = [];
  for m = [32 64 128]
    sol = run1(m, 3);
    st = sol.stats.stages{end};
    errs(end+1) = sol.err.linf; %#ok<SAGROW>
    say('  m=%4d: err=%.4e  V=%2d  factor=%.3f  rel=%.1e  mass=%.1e\n', ...
      m, sol.err.linf, st.vcycles, st.avgFactor, st.finalRel, sol.mass.rel);
    if st.avgFactor > 0.15 || st.finalRel > 1e-9
      say('  ** FAIL (%s stalled: the compatibility projection is not working)\n', name);
      PASS = false;
    end
  end
  r = log2(errs(1:end-1) ./ errs(2:end));
  say('  rates: %s\n', mat2str(r, 3));
  if any(r < 1.8), say('  ** FAIL (%s rates)\n', name); PASS = false; end
end

if PASS
  say('\nNEUMANN COMPATIBILITY TESTS PASSED\n');
else
  error('test_neumann_compat: failed');
end
