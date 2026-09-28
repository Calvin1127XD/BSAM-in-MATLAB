%RUNPARABOLICTESTS ParabolicBSAM verification suite.
%   P1 space-time MMS, traveling Gaussian: temporal orders (BE 1,
%      BDF2 2, CN 2) under (h, Dt) -> (h, Dt)/2, with the moving mesh
%      regridding every 2 steps and time-dependent Dirichlet data.
%   P2 heat equation (Neumann, f = 0, C = 0): discrete mass conservation
%      over 200 steps including ~40 conservative regrids (with
%      derefinement as the solution smooths out).
%   P3 moving source: the hierarchy tracks the source path.
%   P4 efficiency/stability: warm-start V-cycles per step; large-Dt BE.
clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'src', 'ebpAddpaths.m'));
outdir = fullfile(here, '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end
warning('off', 'eb_tag_cluster:dropthin');
R = struct();
PASS = true;
t0all = tic;

% ======================================================================
fprintf('\n=== P1: traveling-bump MMS, temporal orders =======================\n');
sg = 0.06; vx = 0.4; vy = 0.3;
xc = @(t) 0.30 + vx*t;  yc = @(t) 0.35 + vy*t;
uex = @(x,y,t) exp(-((x-xc(t)).^2 + (y-yc(t)).^2)/(2*sg^2));
ut  = @(x,y,t) uex(x,y,t) .* ((x-xc(t))*vx + (y-yc(t))*vy)/sg^2;
lap = @(x,y,t) uex(x,y,t) .* (((x-xc(t)).^2 + (y-yc(t)).^2)/sg^4 - 2/sg^2);
ff  = @(x,y,t) ut(x,y,t) - lap(x,y,t) + uex(x,y,t);   % du/dt - lap u + u = f
pprob = ebpProblem([0 1 0 1], 1, 1, ff, {'dirichlet', uex}, ...
                    @(x,y) uex(x,y,0), uex);
T = 0.4;

% ---- P1a: temporal orders, isolated by Dt-Cauchy differences on a FIXED
% uniform 128^2 mesh (the spatial error cancels exactly in differences of
% solutions on the same grid)
fprintf('  -- P1a: Dt-Cauchy on fixed uniform 128^2 --\n');
schemes = {'be', 'bdf2', 'cn'};
dts = 0.04 ./ 2.^(0:3);
P1a = nan(3, numel(dts)-1);
for si = 1:3
  U = cell(1, numel(dts));
  for j = 1:numel(dts)
    sol = ebpSolve(pprob, 'Scheme', schemes{si}, 'Dt', dts(j), ...
      'TFinal', T, 'RegridEvery', 0, 'StepVerbose', 0, ...
      'BaseCells', [128 128], 'MaxLevels', 1, 'RelTol', 1e-10);
    Hf = sol.H;
    U{j} = Hf.lev{Hf.nlev}.Q{1}(2:end-1, 2:end-1);
  end
  for j = 1:numel(dts)-1
    P1a(si, j) = max(max(abs(U{j} - U{j+1})));
  end
  rate = log2(P1a(si, 1:end-1) ./ P1a(si, 2:end));
  fprintf('  %-4s  Cauchy diffs %s  rates %s\n', schemes{si}, ...
    mat2str(P1a(si, :), 3), mat2str(rate, 3));
end
r1a = log2(P1a(:, 1:end-1) ./ P1a(:, 2:end));
if r1a(1, end) < 0.9,  fprintf('  ** P1a FAIL (BE order)\n'); PASS = false; end
if r1a(2, end) < 1.8,  fprintf('  ** P1a FAIL (BDF2 order)\n'); PASS = false; end
if r1a(3, end) < 1.8,  fprintf('  ** P1a FAIL (CN order)\n'); PASS = false; end
R.P1a = P1a; R.P1arates = r1a;

% ---- P1b: joint (h, Dt) refinement WITH the moving adaptive mesh and
% regridding every 2 steps; BDF2 should show the combined 2nd order
fprintf('  -- P1b: joint refinement with AMR + regrids (BDF2) --\n');
P1b = nan(1, 3);
for j = 0:2
  t = tic;
  % ulap ~ h^2|lap u|: scale the threshold with the root h^2 so the SAME
  % physical region is tagged at every resolution of the sweep
  thrj = 5e-2 / 4^j;
  sol = ebpSolve(pprob, 'Scheme', 'bdf2', ...
    'Dt', 0.05/2^j, 'TFinal', T, 'RegridEvery', 2, 'StepVerbose', 0, ...
    'BaseCells', 32*2^j*[1 1], 'MaxLevels', 2, ...
    'Indicator', 'ulap', 'TagThreshold', @(L) thrj/4^(L-1), ...
    'TagBuffer', 3, 'RelTol', 1e-9);
  P1b(j+1) = sol.err.linf;
  fprintf('  bdf2  h=1/%-3d Dt=%-7.4g  err_inf=%.4e  <V>=%.1f  regrids=%d  t=%.1fs\n', ...
    64*2^j, 0.05/2^j, sol.err.linf, mean(sol.cycles), sol.nregrids, toc(t));
end
r1b = log2(P1b(1:end-1) ./ P1b(2:end));
fprintf('  joint rates: %s\n', mat2str(r1b, 3));
if r1b(end) < 1.5
  fprintf('  ** P1b FAIL (joint order with moving mesh)\n'); PASS = false;
end
R.P1b = P1b; R.P1brates = r1b;

% ======================================================================
fprintf('\n=== P2: heat equation mass conservation across regrids ============\n');
D2 = @(x,y) 1 + 0.5*sin(2*pi*x).*cos(2*pi*y);
u02 = @(x,y) 0.5 + exp(-((x-0.35).^2 + (y-0.4).^2)/(2*0.03^2)) ...
           + 0.7*exp(-((x-0.7).^2 + (y-0.65).^2)/(2*0.04^2));
pprob2 = ebpProblem([0 1 0 1], D2, 0, 0, {'neumann', 0}, u02);
sol2 = ebpSolve(pprob2, 'Scheme', 'bdf2', 'Dt', 2e-3, 'TFinal', 0.4, ...
  'RegridEvery', 5, 'StepVerbose', 50, ...
  'BaseCells', [64 64], 'MaxLevels', 3, ...
  'Indicator', 'ulap', 'TagThreshold', @(L) 1e-3/4^(L-1), ...
  'TagBuffer', 2, 'RelTol', 1e-9);
drift = max(abs(sol2.mass - sol2.mass(1)));
fprintf('  steps=%d  regrids=%d  <V/step>=%.1f  mass drift=%.3e (mass=%.4f)\n', ...
  numel(sol2.cycles), sol2.nregrids, mean(sol2.cycles), drift, sol2.mass(1));
fprintf('  final levels: %d (started with 3 - derefinement as u smooths)\n', ...
  numel(sol2.levels));
if drift > 1e-9 * abs(sol2.mass(1))
  fprintf('  ** P2 FAIL (mass drift)\n'); PASS = false;
end
R.P2 = struct('drift', drift, 'mass0', sol2.mass(1), ...
  'meanV', mean(sol2.cycles), 'nregrids', sol2.nregrids);

% ======================================================================
fprintf('\n=== P3: moving heat source - mesh tracking ========================\n');
om = 2*pi; Rc = 0.22; ctr = [0.5 0.5]; ss = 0.03;
cx = @(t) ctr(1) + Rc*cos(om*t);  cy = @(t) ctr(2) + Rc*sin(om*t);
src = @(x,y,t) 25*exp(-((x-cx(t)).^2 + (y-cy(t)).^2)/(2*ss^2));
pprob3 = ebpProblem([0 1 0 1], 1, 0, src, {'dirichlet', 0}, @(x,y) 0*x);
sol3 = ebpSolve(pprob3, 'Scheme', 'bdf2', 'Dt', 5e-3, 'TFinal', 1.0, ...
  'RegridEvery', 4, 'StepVerbose', 0, ...
  'BaseCells', [64 64], 'MaxLevels', 3, ...
  'Indicator', 'ulap', 'TagThreshold', @(L) 2e-3/4^(L-1), ...
  'TagBuffer', 3, 'RelTol', 1e-8);
% the finest level must cover the final source position
H3 = sol3.H;
lev = H3.lev{H3.nlev};
dom = [0 1 0 1];
covered = false;
for p = 1:lev.np
  b = lev.box(p, :);
  x0 = dom(1) + (b(1)-1)*lev.h(1); x1 = dom(1) + b(2)*lev.h(1);
  y0 = dom(3) + (b(3)-1)*lev.h(2); y1 = dom(3) + b(4)*lev.h(2);
  if cx(1.0) >= x0 && cx(1.0) <= x1 && cy(1.0) >= y0 && cy(1.0) <= y1
    covered = true;
  end
end
nc = 0;
for L = 1:numel(sol3.levels), nc = nc + sol3.levels(L).cells; end
fprintf('  source at t=1: (%.3f, %.3f); finest level covers it: %d\n', ...
  cx(1.0), cy(1.0), covered);
fprintf('  cells=%d (uniform 256^2 would be %d), regrids=%d, <V>=%.1f\n', ...
  nc, 256^2, sol3.nregrids, mean(sol3.cycles));
if ~covered, fprintf('  ** P3 FAIL (mesh lost the source)\n'); PASS = false; end
ebPlot(H3, 'mesh', fullfile(outdir, 'p3_tracking_mesh.png'), ...
  'Mesh tracking the rotating source (t = 1)');
ebPlot(H3, 'solution', fullfile(outdir, 'p3_tracking_solution.png'), ...
  'Temperature at t = 1');
R.P3 = struct('covered', covered, 'cells', nc);

% ======================================================================
fprintf('\n=== P4: warm-start efficiency and large-Dt stability ==============\n');
fprintf('  P2 mean V-cycles/step: %.2f\n', mean(sol2.cycles));
if mean(sol2.cycles) > 4.5
  fprintf('  ** P4 FAIL (too many cycles per step)\n'); PASS = false;
end
solS = ebpSolve(pprob2, 'Scheme', 'be', 'Dt', 0.1, 'TFinal', 0.5, ...
  'RegridEvery', 0, 'StepVerbose', 0, 'BaseCells', [64 64], ...
  'MaxLevels', 1, 'RelTol', 1e-9);
umax = -inf;
Hf = solS.H;
for p = 1:Hf.lev{Hf.nlev}.np
  umax = max(umax, max(abs(Hf.lev{Hf.nlev}.Q{p}(:))));
end
fprintf('  BE with Dt=0.1 (5 steps): max|u| = %.3f (bounded, no ringing)\n', umax);
if ~isfinite(umax) || umax > 10
  fprintf('  ** P4 FAIL (large-Dt instability)\n'); PASS = false;
end
R.P4 = struct('meanV', mean(sol2.cycles), 'umaxBigDt', umax);

% ======================================================================
R.elapsed = toc(t0all);
save(fullfile(outdir, 'results_parabolic.mat'), 'R');
fprintf('\nTotal parabolic suite time: %.1f s\n', R.elapsed);
if PASS
  fprintf('ALL PARABOLIC TESTS PASSED\n');
else
  error('run_parabolic_tests: at least one test failed');
end
