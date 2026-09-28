%runAllTests EllipticBSAM verification suite.  Writes results to
%   ../output/results.mat, figures to ../output/, and a log to stdout.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
outdir = fullfile(here, '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end
R = struct();
t0all = tic;
PASS = true;

logf = @(varargin) fprintf(varargin{:});

% ======================================================================
logf('\n=== T1: uniform MMS, variable D(x,y), C(x,y), Dirichlet =========\n');
cs = mmsCases('vardc');
ms = [32 64 128 256];
T1 = [];
for i = 1:numel(ms)
  opts = ebOptions('BaseCells', ms(i)*[1 1], 'Verbose', 0, 'RelTol', 1e-11);
  t = tic; sol = ebSolve(cs.prob, opts); tt = toc(t);
  st = sol.stats.stages{end};
  T1(i,:) = [ms(i), sol.err.linf, sol.err.l2, st.avgFactor, st.vcycles, tt]; %#ok<*SAGROW>
  logf('  m=%4d  err_inf=%.4e  err_l2=%.4e  factor=%.3f  V=%d  t=%.2fs\n', ...
    ms(i), sol.err.linf, sol.err.l2, st.avgFactor, st.vcycles, tt);
end
r1 = log2(T1(1:end-1,2)./T1(2:end,2));
logf('  L_inf rates: %s\n', mat2str(r1.', 3));
if any(r1 < 1.85), logf('  ** T1 FAIL\n'); PASS = false; end
R.T1 = T1; R.T1rates = r1;

% ======================================================================
logf('\n=== T2: uniform pure-Neumann singular Poisson ====================\n');
cs = mmsCases('pureneumann');
T2 = [];
for i = 1:numel(ms)
  opts = ebOptions('BaseCells', ms(i)*[1 1], 'Verbose', 0, 'RelTol', 1e-11);
  sol = ebSolve(cs.prob, opts);
  st = sol.stats.stages{end};
  T2(i,:) = [ms(i), sol.err.linf, sol.err.l2, st.avgFactor, st.vcycles];
  logf('  m=%4d  err_inf=%.4e  err_l2=%.4e  factor=%.3f  V=%d\n', ...
    ms(i), sol.err.linf, sol.err.l2, st.avgFactor, st.vcycles);
end
r2 = log2(T2(1:end-1,2)./T2(2:end,2));
logf('  L_inf rates: %s\n', mat2str(r2.', 3));
if any(r2 < 1.85), logf('  ** T2 FAIL\n'); PASS = false; end
R.T2 = T2; R.T2rates = r2;

% ======================================================================
logf('\n=== T3: fixed box-in-box hierarchies (paper Tables 2/3) ==========\n');
cs = mmsCases('gauss51');
% Concentric square tag regions: level L refines a fixed centered box, so
% the hierarchy reproduces the paper's box-in-box meshes exactly.
tagfn = @(X, Y, L) ...
  (L == 1) .* (max(abs(X), abs(Y)) < 0.1875) + ...
  (L == 2) .* (max(abs(X), abs(Y)) < 0.125) + ...
  (L == 3) .* (max(abs(X), abs(Y)) < 0.0625);
m0s = [32 64 128 256];
for meth = ["qif", "lif"]
  T3 = nan(numel(m0s), 7);
  for i = 1:numel(m0s)
    for nlv = 2:4
      opts = ebOptions('BaseCells', m0s(i)*[1 1], 'MaxLevels', nlv, ...
        'Indicator', 'fun', 'TagFunction', tagfn, 'TagBuffer', 0, ...
        'GhostInterp', char(meth), 'Verbose', 0, 'RelTol', 1e-10, ...
        'MaxVCycles', 40);
      t = tic; sol = ebSolve(cs.prob, opts); tt = toc(t);
      st = sol.stats.stages{end};
      T3(i, 1) = m0s(i);
      T3(i, nlv) = sol.err.linf;          % cols 2,3,4 = 2,3,4-level err
      T3(i, nlv+3) = st.avgFactor;        % cols 5,6,7 = factors
      if nlv == 3, mass3(i) = sol.mass.rel; end
      logf('  [%s] m0=%4d L=%d: err=%.4e factor=%.3f mass=%.1e V=%d t=%.1fs\n', ...
        meth, m0s(i), nlv, sol.err.linf, st.avgFactor, sol.mass.rel, ...
        st.vcycles, tt);
    end
  end
  rates = log2(T3(1:end-1, 2:4) ./ T3(2:end, 2:4));
  logf('  [%s] rates (2/3/4-level): %s\n', meth, mat2str(rates, 3));
  if any(rates(:) < 1.8), logf('  ** T3 %s FAIL\n', meth); PASS = false; end
  R.(sprintf('T3_%s', meth)) = T3;
  R.(sprintf('T3_%s_rates', meth)) = rates;
end
R.T3_mass3 = mass3;

% ======================================================================
logf('\n=== T4: L-shaped refinement patch (paper Table 5) ================\n');
cs = mmsCases('gauss52');
Lshape = @(X, Y) (X >= 0.375 & X <= 0.625 & Y >= 0.375 & Y <= 0.5) | ...
                 (X >= 0.375 & X <= 0.5   & Y >= 0.5   & Y <= 0.625);
m0s4 = [16 32 64 128 256];
T4 = [];
for i = 1:numel(m0s4)
  opts = ebOptions('BaseCells', m0s4(i)*[1 1], 'MaxLevels', 2, ...
    'Indicator', 'fun', 'TagFunction', @(X,Y,L) Lshape(X,Y), ...
    'TagBuffer', 0, 'MinFillRatio', 0.9, 'Verbose', 0, 'RelTol', 1e-10);
  sol = ebSolve(cs.prob, opts);
  st = sol.stats.stages{end};
  np2 = sol.levels(2).np;
  T4(i,:) = [m0s4(i), sol.err.linf, st.avgFactor, np2, sol.mass.rel];
  logf('  m0=%4d  err_inf=%.4e  factor=%.3f  level-2 patches=%d  mass=%.1e\n', ...
    m0s4(i), sol.err.linf, st.avgFactor, np2, sol.mass.rel);
  if i == 3
    ebPlot(sol.H, 'mesh', fullfile(outdir, 't4_lshape_mesh.png'), ...
      'L-shaped two-level mesh (m_0=64)');
    ebPlot(sol.H, 'solution', fullfile(outdir, 't4_lshape_solution.png'), ...
      'Composite solution, L-shaped refinement');
  end
end
r4 = log2(T4(1:end-1,2)./T4(2:end,2));
logf('  L_inf rates: %s\n', mat2str(r4.', 3));
if any(r4(2:end) < 1.8), logf('  ** T4 FAIL\n'); PASS = false; end
% at m0=16 the min-width constraint correctly yields a single box; from
% m0=32 on, the L must be decomposed into >= 2 patches (reentrant corner)
if any(T4(2:end,4) < 2), logf('  ** T4 FAIL (expected >= 2 patches)\n'); PASS = false; end
R.T4 = T4; R.T4rates = r4;

% ======================================================================
logf('\n=== T5: fully adaptive (undivided-Laplacian tagging) =============\n');
cs = mmsCases('threebumps');
T5 = [];
for nlv = 1:4
  % undivided Laplacian scales like h^2 |lap u|: scale the threshold by
  % 1/4 per level so the same physical feature stays tagged when refined
  opts = ebOptions('BaseCells', [64 64], 'MaxLevels', nlv, ...
    'Indicator', 'ulap', 'TagThreshold', @(L) 1e-2/4^(L-1), 'TagBuffer', 2, ...
    'Verbose', 0, 'RelTol', 1e-10, 'MaxVCycles', 40);
  t = tic; sol = ebSolve(cs.prob, opts); tt = toc(t);
  st = sol.stats.stages{end};
  nA = numel(sol.levels);
  T5(nlv,:) = [nlv, nA, sol.err.linf, sol.err.l2, st.avgFactor, ...
               sol.err.dof, tt];
  logf('  MaxLevels=%d (got %d): err_inf=%.4e  err_l2=%.4e  factor=%.3f  dof=%d  t=%.1fs\n', ...
    nlv, nA, sol.err.linf, sol.err.l2, st.avgFactor, sol.err.dof, tt);
  if nlv == 4
    ebPlot(sol.H, 'mesh', fullfile(outdir, 't5_adaptive_mesh.png'), ...
      'Adaptive 4-level mesh, three-bump problem');
    ebPlot(sol.H, 'solution', fullfile(outdir, 't5_adaptive_solution.png'), ...
      'Composite solution, three bumps');
  end
end
dec = T5(2:end,3) ./ T5(1:end-1,3);
logf('  err_inf decrease factors per added level: %s\n', mat2str(dec.', 3));
if any(dec > 0.6), logf('  ** T5 FAIL (insufficient error reduction)\n'); PASS = false; end
R.T5 = T5;

% ======================================================================
logf('\n=== T6: h-independence of the V-cycle (2-level, paper Fig 8) =====\n');
% from T3 qif: factors col 5 (2-level)
f2 = R.T3_qif(:, 5);
logf('  m0 = %s\n  factors = %s\n', mat2str(m0s), mat2str(f2.', 3));
if max(f2) > 0.2
  logf('  ** T6 FAIL (factor too large)\n'); PASS = false;
end
R.T6 = [m0s(:), f2];

% residual history figure (m0 = 256, 3-level, fresh run, verbose capture)
cs6 = mmsCases('gauss51');
hists = {};
for m0 = [64 128 256]
  opts = ebOptions('BaseCells', m0*[1 1], 'MaxLevels', 3, ...
    'Indicator', 'fun', 'TagFunction', tagfn, 'TagBuffer', 0, ...
    'Verbose', 0, 'RelTol', 1e-12, 'MaxVCycles', 14);
  sol = ebSolve(cs6.prob, opts);
  hists{end+1} = sol.stats.stages{end}.resHist;
end
fig = figure('Visible', 'off');
co = lines(numel(hists)); hold on;
for i = 1:numel(hists)
  semilogy(0:numel(hists{i})-1, hists{i}/hists{i}(1), '-o', 'Color', co(i,:));
end
set(gca, 'YScale', 'log'); grid on;
legend({'m_0=64', 'm_0=128', 'm_0=256'});
xlabel('V-cycle'); ylabel('|r|_\infty / |r^0|_\infty');
title('h-independent geometric convergence (3-level mesh)');
exportgraphics(fig, fullfile(outdir, 't6_vcycle_history.png'), 'Resolution', 150);
close(fig);
R.T6hists = hists;

% ======================================================================
logf('\n=== T7: discrete mass conservation (QIF vs no corrections) =======\n');
cs = mmsCases('gauss51');
mass = nan(2, 3);
for nlv = 2:3
  for mc = [true false]
    opts = ebOptions('BaseCells', [64 64], 'MaxLevels', nlv, ...
      'Indicator', 'fun', 'TagFunction', tagfn, 'TagBuffer', 0, ...
      'MassCorrection', mc, 'Verbose', 0, 'RelTol', 1e-11, 'MaxVCycles', 60);
    sol = ebSolve(cs.prob, opts);
    mass(2-mc, nlv-1) = sol.mass.rel;
    logf('  levels=%d  corrections=%d : mass defect (rel) = %.3e   err=%.3e\n', ...
      nlv, mc, sol.mass.rel, sol.err.linf);
  end
end
if mass(1,1) > 1e-12 || mass(1,2) > 1e-12
  logf('  ** T7 FAIL (conservative scheme defect too large)\n'); PASS = false;
end
if mass(2,1) < 1e3*mass(1,1)
  logf('  note: uncorrected defect not dramatically larger (check)\n');
end
R.T7 = mass;

% ======================================================================
logf('\n=== T8: pathological scattered tagging robustness =================\n');
cs = mmsCases('vardc');
% 24 random tag disks of radius 0.02: produces many small fragmented
% clusters, stressing the clustering/nesting/min-width machinery.
rng(42);
pts = rand(24, 2);
scatfn = @(X, Y, L) reshape(any( ...
  (X(:) - pts(:,1).').^2 + (Y(:) - pts(:,2).').^2 < 0.0004, 2), size(X));
opts = ebOptions('BaseCells', [64 64], 'MaxLevels', 3, ...
  'Indicator', 'fun', 'TagFunction', scatfn, 'TagBuffer', 1, ...
  'MinBoxWidth', 4, 'Verbose', 0, 'RelTol', 1e-10, 'MaxVCycles', 60);
sol = ebSolve(cs.prob, opts);
minw = inf; npat = 0;
for L = 2:numel(sol.levels)
  bx = sol.levels(L).box;
  npat = npat + size(bx, 1);
  minw = min([minw; bx(:,2)-bx(:,1)+1; bx(:,4)-bx(:,3)+1]);
end
st = sol.stats.stages{end};
logf('  patches=%d  min width=%d (fine cells)  factor=%.3f  err=%.4e  mass=%.1e\n', ...
  npat, minw, st.avgFactor, sol.err.linf, sol.mass.rel);
ebPlot(sol.H, 'mesh', fullfile(outdir, 't8_scattered_mesh.png'), ...
  'Scattered-tag robustness mesh');
if st.finalRel > 1e-10 || st.avgFactor > 0.3
  logf('  ** T8 FAIL (solver degraded on fragmented patches)\n'); PASS = false;
end
% hard floor for valid interface stencils is 2 parent = 4 fine cells;
% MinBoxWidth=4 is a target the nesting clip may locally undercut
if minw < 4, logf('  ** T8 FAIL (min patch width < 4 fine cells)\n'); PASS = false; end
R.T8 = struct('npat', npat, 'minw', minw, 'factor', st.avgFactor, ...
              'err', sol.err.linf);

% ======================================================================
R.elapsed = toc(t0all);
save(fullfile(outdir, 'results.mat'), 'R');
logf('\nTotal suite time: %.1f s\n', R.elapsed);
if PASS
  logf('ALL FOLDER-1 TESTS PASSED\n');
else
  error('run_all_tests: at least one test failed');
end
