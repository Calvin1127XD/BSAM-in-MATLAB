function makeConvergenceStudy()
%makeConvergenceStudy  Numerical-ground-truth convergence study + figures.
%
%   Reproduces the discrete-truth benchmark for the standalone report
%   (EllipticBSAM/convergence_study/).  One uniform 1024^2 solve provides a
%   numerical ground truth u_fine; the adaptive composite solver is then
%   measured AGAINST u_fine (not against the analytic MMS solution) along
%   two orthogonal refinement axes:
%
%     (A) BREADTH  -- fix the depth (MaxLevels = 4), lower the tag threshold
%                     so the finest-cell COVERAGE grows; watch the error fall.
%     (B) DEPTH    -- fix the threshold (0.0003), add adaptive LEVELS so the
%                     tagged regions refine deeper; watch the error fall.
%
%   Error is the composite discrete-truth norm: every visible composite cell
%   is compared to the uniform-fine solution conservatively block-restricted
%   to that cell's resolution.  L-inf = max|e|;  L2 = sqrt( sum h_k^2 e^2 ).
%
%   Writes figs/{surface,mesh,coverage,levels}.png|pdf and prints every
%   number used in the report.  Rerun:  >> makeConvergenceStudy

here   = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
figdir = fullfile(here, '..', 'convergence_study', 'figs');
if ~exist(figdir, 'dir'), mkdir(figdir); end

% -------------------- parameters ---------------------------------------
rootFine = 1024;             % uniform ground-truth resolution
rootAMR  = 128;              % adaptive root resolution
maxLev   = 4;                % root + 3 refinements -> 1024 at finest
buf      = 3;
relTol   = 1e-12;            % solve everything below any discretization scale
thrFix   = 0.0003;           % fixed threshold for the DEPTH study (B)
thrSweep = [0.02 0.005 0.0012 0.0003 0.0001 0.00003];  % BREADTH study (A)

% -------------------- manufactured solution (homogeneous Dirichlet) ----
%   u = 0.1 sin(pi x) sin(2pi y) + 4 interior Gaussian bumps; sigmas chosen
%   so every bump is << discretization error on the boundary, hence u=0 on
%   the boundary to machine-relevant accuracy.
xi = [0.50 0.50 0.30 0.80]; yi = [0.50 0.55 0.30 0.80];
si = [0.08 0.04 0.04 0.02]; Ai = [1.00 5.00 -1.00 0.50];
g  = @(x,y,a,b,s) exp(-((x-a).^2 + (y-b).^2)/(2*s^2));
lg = @(x,y,a,b,s) g(x,y,a,b,s).*(((x-a).^2 + (y-b).^2)/s^4 - 2/s^2);
ss = @(x,y) 0.1*sin(pi*x).*sin(2*pi*y);
lss= @(x,y) -0.5*pi^2*sin(pi*x).*sin(2*pi*y);
uex= @(x,y) ss(x,y) + Ai(1)*g(x,y,xi(1),yi(1),si(1)) + Ai(2)*g(x,y,xi(2),yi(2),si(2)) ...
   + Ai(3)*g(x,y,xi(3),yi(3),si(3)) + Ai(4)*g(x,y,xi(4),yi(4),si(4));
lap= @(x,y) lss(x,y) + Ai(1)*lg(x,y,xi(1),yi(1),si(1)) + Ai(2)*lg(x,y,xi(2),yi(2),si(2)) ...
   + Ai(3)*lg(x,y,xi(3),yi(3),si(3)) + Ai(4)*lg(x,y,xi(4),yi(4),si(4));
f  = @(x,y) uex(x,y) - lap(x,y);                 % -div(grad u) + u = f, D=C=1
prob = ebProblem([0 1 0 1], 1, 1, f, {'dirichlet', 0}, uex);

% -------------------- numerical ground truth ---------------------------
fprintf('Ground truth: uniform %d^2 solve ...\n', rootFine);
t = tic;
sol1  = ebSolve(prob, 'BaseCells', rootFine*[1 1], 'Verbose', 0, ...
  'RelTol', relTol, 'MaxVCycles', 80);
ufine = sol1.H.lev{sol1.H.rootIdx}.Q{1}(2:end-1, 2:end-1);
fprintf('  done (%.1fs); finalRel=%.1e; vs-analytic Linf=%.3e\n', ...
  toc(t), sol1.stats.stages{end}.finalRel, sol1.err.linf);

% ============================ STUDY A: BREADTH =========================
fprintf('\n===== STUDY A: fixed depth (%d levels), sweep threshold =====\n', maxLev);
fprintf(' thr      cover%%   glinf      gl2       fineLinf   vsAnalytic  cells\n');
A = struct('thr',{},'cover',{},'glinf',{},'gl2',{},'fineLinf',{},'analytic',{},'cells',{});
Hrep = [];   % thr=0.0003 hierarchy (surface + high-coverage mesh)
Hlo  = [];   % thr=0.02   hierarchy (localized low-coverage mesh)
for it = 1:numel(thrSweep)
  th = thrSweep(it);
  s = ebSolve(prob, 'BaseCells', rootAMR*[1 1], 'MaxLevels', maxLev, ...
    'Indicator', 'ulap', 'TagThreshold', @(L) th/4^(L-1), 'TagBuffer', buf, ...
    'Verbose', 0, 'RelTol', relTol, 'MaxVCycles', 80);
  [gli, gl2, ncell] = globalDiscErr(s.H, ufine, rootFine);
  [fli, ~, ncov]    = finestAlignedErr(s.H, ufine, rootFine);
  A(it) = struct('thr',th,'cover',100*ncov/rootFine^2,'glinf',gli,'gl2',gl2, ...
    'fineLinf',fli,'analytic',s.err.linf,'cells',ncell);
  fprintf('  %8.5f %5.1f  %.3e  %.3e  %.3e  %.3e  %d\n', ...
    th, A(it).cover, gli, gl2, fli, s.err.linf, ncell);
  if abs(th - thrFix) < 1e-12, Hrep = s.H; end     % save the thr=0.0003 run
  if abs(th - 0.02)  < 1e-12, Hlo  = s.H; end     % save the thr=0.02 run
end

% ============================ STUDY B: DEPTH ===========================
fprintf('\n===== STUDY B: fixed threshold %.4f, sweep #levels =====\n', thrFix);
fprintf(' levels  finest   glinf      gl2       vsAnalytic  cells\n');
levList = 1:maxLev;
B = struct('lev',{},'finest',{},'glinf',{},'gl2',{},'analytic',{},'cells',{});
for it = 1:numel(levList)
  ml = levList(it);
  s = ebSolve(prob, 'BaseCells', rootAMR*[1 1], 'MaxLevels', ml, ...
    'Indicator', 'ulap', 'TagThreshold', @(L) thrFix/4^(L-1), 'TagBuffer', buf, ...
    'Verbose', 0, 'RelTol', relTol, 'MaxVCycles', 80);
  [gli, gl2, ncell] = globalDiscErr(s.H, ufine, rootFine);
  B(it) = struct('lev',ml,'finest',s.H.lev{s.H.nlev}.n(1),'glinf',gli, ...
    'gl2',gl2,'analytic',s.err.linf,'cells',ncell);
  fprintf('  %4d    %5d^2  %.3e  %.3e  %.3e  %d\n', ...
    ml, B(it).finest, gli, gl2, s.err.linf, ncell);
end
fprintf('  (uniform %d^2 ground truth used %d cells)\n', rootFine, rootFine^2);

% ============================ FIGURES ==================================
set(groot, 'defaultAxesFontSize', 13, 'defaultLineLineWidth', 1.8);

% --- Fig 1: 3D surface of the composite solution -----------------------
if isempty(Hrep), Hrep = sol1.H; end
fig = ebPlot(Hrep, 'solution', '', ...
  sprintf('composite solution  (root %d^2 + %d levels, thr=%.4f)', ...
  rootAMR, maxLev-1, thrFix));
exportgraphics(fig, fullfile(figdir, 'surface.png'), 'Resolution', 300);
close(fig);

% --- Fig 2: adaptive patch layout, low vs high coverage ----------------
if isempty(Hlo), Hlo = Hrep; end
covLo = A([A.thr] == 0.02).cover;
fig = ebPlot(Hlo, 'mesh', '', sprintf('thr = 0.02  (%.1f%% finest coverage)', covLo));
exportgraphics(fig, fullfile(figdir, 'mesh_lo.png'), 'Resolution', 300);
close(fig);
covHi = A([A.thr] == thrFix).cover;
fig = ebPlot(Hrep, 'mesh', '', sprintf('thr = %.4f  (%.1f%% finest coverage)', thrFix, covHi));
exportgraphics(fig, fullfile(figdir, 'mesh_hi.png'), 'Resolution', 300);
close(fig);

% --- Fig 3: BREADTH -- error vs finest-cell coverage -------------------
cov = [A.cover]; gli = [A.glinf]; gl2 = [A.gl2]; fli = [A.fineLinf];
fig = figure('Visible', 'off', 'Position', [100 100 760 560]);
ax = axes(fig); hold(ax, 'on');
semilogy(ax, cov, gli, '-o', 'Color', [0.20 0.30 0.75], 'MarkerFaceColor', [0.20 0.30 0.75]);
semilogy(ax, cov, gl2, '-s', 'Color', [0.85 0.45 0.10], 'MarkerFaceColor', [0.85 0.45 0.10]);
semilogy(ax, cov, fli, '--^', 'Color', [0.45 0.45 0.45], 'MarkerFaceColor', [0.45 0.45 0.45]);
set(ax, 'YScale', 'log'); grid(ax, 'on');
xlabel(ax, 'finest-cell coverage of the $1024^2$ domain (\%)', 'Interpreter', 'latex');
ylabel(ax, 'error vs numerical ground truth');
legend(ax, {'global $\|e\|_\infty$', 'global $\|e\|_2$', ...
  'refined-region $\|e\|_\infty$'}, 'Interpreter', 'latex', 'Location', 'southwest');
title(ax, 'Study A: breadth -- lower threshold \rightarrow more coverage', 'Interpreter', 'tex');
for it = 1:numel(A)
  text(ax, cov(it), gli(it)*1.7, sprintf('%.0e', A(it).thr), ...
    'FontSize', 9, 'HorizontalAlignment', 'center', 'Color', [0.20 0.30 0.75]);
end
exportgraphics(fig, fullfile(figdir, 'coverage.pdf'), 'ContentType', 'vector');
close(fig);

% --- Fig 4: DEPTH -- error vs number of adaptive levels ----------------
lv = [B.lev]; bgli = [B.glinf]; bgl2 = [B.gl2]; fin = [B.finest];
fig = figure('Visible', 'off', 'Position', [100 100 760 560]);
ax = axes(fig); hold(ax, 'on');
semilogy(ax, lv, bgli, '-o', 'Color', [0.20 0.30 0.75], 'MarkerFaceColor', [0.20 0.30 0.75]);
semilogy(ax, lv, bgl2, '-s', 'Color', [0.85 0.45 0.10], 'MarkerFaceColor', [0.85 0.45 0.10]);
set(ax, 'YScale', 'log'); grid(ax, 'on');
xlabel(ax, 'number of adaptive levels  (root + refinements)');
ylabel(ax, 'error vs numerical ground truth');
xticks(ax, lv);
xticklabels(ax, arrayfun(@(a,b) sprintf('%d  (%d^2)', a, b), lv, fin, 'UniformOutput', false));
legend(ax, {'global $\|e\|_\infty$', 'global $\|e\|_2$'}, ...
  'Interpreter', 'latex', 'Location', 'northeast');
title(ax, sprintf('Study B: depth -- fixed threshold %.4f', thrFix), 'Interpreter', 'tex');
exportgraphics(fig, fullfile(figdir, 'levels.pdf'), 'ContentType', 'vector');
close(fig);

fprintf('\nFigures written to %s\n', figdir);
fprintf('DONE.\n');
end

% =======================================================================
function [linf, l2, ncell] = globalDiscErr(H, ufine, rootFine) %#ok<INUSD>
%globalDiscErr  Composite error of H vs the block-restricted fine truth,
%   over EVERY visible (uncovered) composite cell.
linf = 0; l2sq = 0; ncell = 0;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  uFL = restrictTo(ufine, lev.n(1));
  da  = lev.h(1) * lev.h(2);
  for p = 1:lev.np
    vis = ~lev.covered{p};
    b   = lev.box(p, :);
    e   = lev.Q{p}(2:end-1, 2:end-1) - uFL(b(1):b(2), b(3):b(4));
    e   = e(vis);
    if isempty(e), continue; end
    linf  = max(linf, max(abs(e)));
    l2sq  = l2sq + sum(e.^2) * da;
    ncell = ncell + numel(e);
  end
end
l2 = sqrt(l2sq);
end

% =======================================================================
function [linf, l2, ncmp] = finestAlignedErr(H, ufine, rootFine)
%finestAlignedErr  Error on the finest level only, compared directly to
%   the aligned subset of the (same-resolution) fine truth.
nlev = H.nlev;
if H.lev{nlev}.n(1) ~= rootFine
  linf = NaN; l2 = NaN; ncmp = 0; return;
end
linf = 0; l2sq = 0; ncmp = 0; da = (1/rootFine)^2;
for p = 1:H.lev{nlev}.np
  b = H.lev{nlev}.box(p, :);
  d = abs(H.lev{nlev}.Q{p}(2:end-1, 2:end-1) - ufine(b(1):b(2), b(3):b(4)));
  linf = max(linf, max(d(:)));
  l2sq = l2sq + sum(d(:).^2) * da;
  ncmp = ncmp + numel(d);
end
l2 = sqrt(l2sq);
end

% =======================================================================
function uc = restrictTo(uf, n)
%restrictTo  Conservative block-average of uf onto an n x n grid.
r = size(uf, 1) / n;
if r == 1, uc = uf; return; end
uc = squeeze(mean(mean(reshape(uf, r, n, r, n), 1), 3));
end
