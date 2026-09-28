function expAmrVsUniform2()
%expAmrVsUniform2  Discrete-truth AMR-vs-uniform experiment, variant B.
%   Same protocol as expAmrVsUniform (uniform 1024^2 "ground truth"
%   vs a 128^2 root + 3 adaptive levels composite), but with a manufactured
%   solution that has a GENTLE, low-amplitude smooth background and SHARP,
%   well-localized bumps -- i.e. real scale separation, the regime where
%   block-structured AMR is supposed to win.
%
%       u(x,y) = 0.1 sin(pi x) sin(2pi y)
%                + sum_i A_i exp(-((x-a_i)^2+(y-b_i)^2)/(2 sigma_i^2))
%
%   bumps (center, sigma, height):
%       (0.50,0.50)  0.08   1.0
%       (0.50,0.55)  0.04   5.0
%       (0.30,0.30)  0.04  -1.0
%       (0.80,0.80)  0.02   0.5
%   (sigmas are the user's originals divided by 5 so every bump decays to
%   << discretization error at the boundary; widest bump's boundary value
%   ~3.4e-9).  The sine background vanishes on all four sides, so the
%   problem is genuinely homogeneous Dirichlet, u = 0 on the boundary.
%
%   Metrics reported (identical to variant A):
%     1. ||u_comp - u_fine||  at the aligned finest-cell subset (the pure
%        composite-vs-uniform discretization DIFFERENCE; no analytic floor).
%     2. per-level |e|_{L,inf}/|e|_{L,2} vs the restricted fine truth.
%     3. a threshold sweep: as more of the smooth background is refined to
%        the finest level, watch the discrete-truth error collapse.
%
%   Rerun: just  >> expAmrVsUniform2   (edit the params block to vary).

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

% -------------------- parameters (edit to vary) ------------------------
rootFine = 1024;             % Case 1 uniform resolution (the ground truth)
rootAMR  = 128;              % Case 2 root resolution
maxLev   = 4;                % 4 AMR levels = root + 3 refinements -> 1024
thr      = @(L) 0.02 / 4^(L-1);   % undivided-Laplacian per-level threshold
buf      = 3;
relTol   = 1e-12;            % solve both well below any discretization scale

% -------------------- problem: homogeneous Dirichlet -------------------
%   gentle sine background (amplitude 0.1) + four interior Gaussian bumps.
xi = [0.50 0.50 0.30 0.80]; yi = [0.50 0.55 0.30 0.80];
si = [0.08 0.04 0.04 0.02]; Ai = [1.00 5.00 -1.00 0.50];   % sigmas already /5
g  = @(x,y,a,b,s) exp(-((x-a).^2 + (y-b).^2)/(2*s^2));
lg = @(x,y,a,b,s) g(x,y,a,b,s).*(((x-a).^2 + (y-b).^2)/s^4 - 2/s^2);
ss = @(x,y) 0.1*sin(pi*x).*sin(2*pi*y);
lss= @(x,y) -0.5*pi^2*sin(pi*x).*sin(2*pi*y);              % laplacian of ss
uex= @(x,y) ss(x,y) + Ai(1)*g(x,y,xi(1),yi(1),si(1)) + Ai(2)*g(x,y,xi(2),yi(2),si(2)) ...
   + Ai(3)*g(x,y,xi(3),yi(3),si(3)) + Ai(4)*g(x,y,xi(4),yi(4),si(4));
lap= @(x,y) lss(x,y) + Ai(1)*lg(x,y,xi(1),yi(1),si(1)) + Ai(2)*lg(x,y,xi(2),yi(2),si(2)) ...
   + Ai(3)*lg(x,y,xi(3),yi(3),si(3)) + Ai(4)*lg(x,y,xi(4),yi(4),si(4));
f  = @(x,y) uex(x,y) - lap(x,y);                 % -div(grad u) + u = f, D=C=1
prob = ebProblem([0 1 0 1], 1, 1, f, {'dirichlet', 0}, uex);

% -------------------- Case 1: uniform fine "ground truth" --------------
fprintf('Case 1: uniform %d^2 ground-truth solve ...\n', rootFine);
t = tic;
sol1 = ebSolve(prob, 'BaseCells', rootFine*[1 1], 'Verbose', 0, ...
  'RelTol', relTol, 'MaxVCycles', 80);
ufine = sol1.H.lev{sol1.H.rootIdx}.Q{1}(2:end-1, 2:end-1);   % rootFine x rootFine
fprintf('  done (%.1fs); finalRel = %.1e; vs analytic Linf = %.3e\n', ...
  toc(t), sol1.stats.stages{end}.finalRel, sol1.err.linf);

% -------------------- Case 2: adaptive composite -----------------------
fprintf('Case 2: %d^2 root + %d adaptive levels ...\n', rootAMR, maxLev-1);
t = tic;
sol2 = ebSolve(prob, 'BaseCells', rootAMR*[1 1], 'MaxLevels', maxLev, ...
  'Indicator', 'ulap', 'TagThreshold', thr, 'TagBuffer', buf, ...
  'Verbose', 0, 'RelTol', relTol, 'MaxVCycles', 80);
H2 = sol2.H; nlev = H2.nlev;
fprintf('  done (%.1fs); finalRel = %.1e; vs analytic Linf = %.3e\n', ...
  toc(t), sol2.stats.stages{end}.finalRel, sol2.err.linf);

% sanity: the finest level must be exactly the ground-truth resolution
if H2.lev{nlev}.n(1) ~= rootFine
  error('exp_amr_vs_uniform2:res', ...
    'finest level is %d^2, not %d^2 -- adjust threshold/levels so the bumps refine fully.', ...
    H2.lev{nlev}.n(1), rootFine);
end

% -------------------- direct comparison at the aligned subset ----------
linf = 0; l2sq = 0; ncmp = 0; da = (1/rootFine)^2;
for p = 1:H2.lev{nlev}.np
  b = H2.lev{nlev}.box(p, :);
  uc = H2.lev{nlev}.Q{p}(2:end-1, 2:end-1);
  uf = ufine(b(1):b(2), b(3):b(4));            % aligned subset of Case 1
  d  = abs(uc - uf);
  linf = max(linf, max(d(:)));
  l2sq = l2sq + sum(d(:).^2) * da;
  ncmp = ncmp + numel(d);
end

fprintf('\n================= DISCRETE-TRUTH COMPARISON =================\n');
fprintf('finest-level cells compared: %d  (%.1f%% of the %d^2 grid)\n', ...
  ncmp, 100*ncmp/rootFine^2, rootFine);
fprintf('||u_comp - u_fine||_inf  = %.3e   <-- the headline number\n', linf);
fprintf('||u_comp - u_fine||_2    = %.3e\n', sqrt(l2sq));
fprintf('-------------------------------------------------------------\n');
fprintf('for reference, vs the ANALYTIC solution:\n');
fprintf('   Case 1 (uniform 1024^2): Linf = %.3e\n', sol1.err.linf);
fprintf('   Case 2 (adaptive)      : Linf = %.3e\n', sol2.err.linf);
fprintf('mass defect: case1 %.1e, case2 %.1e\n', sol1.mass.rel, sol2.mass.rel);

% -------------------- per-level error breakdown ------------------------
fprintf('\n=========== PER-LEVEL ERROR vs fine ground truth ===========\n');
fprintf('  level L | mesh    | #cells  | |e|_{L,inf} | |e|_{L,2}\n');
r0 = H2.rootIdx;
for L = 0:(nlev - r0)
  k = r0 + L;
  nL = H2.lev{k}.n(1);
  uFL = restrictTo(ufine, nL);               % fine truth averaged to level L
  li = 0; l2s = 0; nc = 0; daL = (1/nL)^2;
  for p = 1:H2.lev{k}.np
    b = H2.lev{k}.box(p, :);
    d = abs(H2.lev{k}.Q{p}(2:end-1, 2:end-1) - uFL(b(1):b(2), b(3):b(4)));
    li = max(li, max(d(:)));
    l2s = l2s + sum(d(:).^2) * daL;
    nc = nc + numel(d);
  end
  fprintf('   L = %d  | %4d^2  | %7d | %.3e  | %.3e\n', L, nL, nc, li, sqrt(l2s));
end

% -------------------- diagnosis sweep: coverage vs floor ---------------
fprintf('\n========== SWEEP: refine more of the background ============\n');
fprintf(' threshold | finest-cell coverage | ||u_comp-u_fine||_inf | vs analytic\n');
for th = [0.02 0.005 0.0012 0.0003 0.0001 0.00003]
  s = ebSolve(prob, 'BaseCells', rootAMR*[1 1], 'MaxLevels', maxLev, ...
    'Indicator', 'ulap', 'TagThreshold', @(L) th/4^(L-1), 'TagBuffer', buf, ...
    'Verbose', 0, 'RelTol', relTol, 'MaxVCycles', 80);
  [li, ~, nc] = compareToFine(s.H, ufine, rootFine);
  fprintf('  %8.5f |      %5.1f%%          |      %.3e        |  %.3e\n', ...
    th, 100*nc/rootFine^2, li, s.err.linf);
end
end

% ------------------------------------------------------------------------
function [linf, l2, ncmp] = compareToFine(H2, ufine, rootFine)
nlev = H2.nlev;
if H2.lev{nlev}.n(1) ~= rootFine
  linf = NaN; l2 = NaN; ncmp = 0; return;   % did not refine fully
end
linf = 0; l2sq = 0; ncmp = 0; da = (1/rootFine)^2;
for p = 1:H2.lev{nlev}.np
  b = H2.lev{nlev}.box(p, :);
  d = abs(H2.lev{nlev}.Q{p}(2:end-1, 2:end-1) - ufine(b(1):b(2), b(3):b(4)));
  linf = max(linf, max(d(:)));
  l2sq = l2sq + sum(d(:).^2) * da;
  ncmp = ncmp + numel(d);
end
l2 = sqrt(l2sq);
end

% ------------------------------------------------------------------------
function uc = restrictTo(uf, n)
%restrictTo  Conservative (block-average) restriction of the fine field uf
%   onto an n x n grid; n must divide size(uf,1).  Identity when n == size.
r = size(uf, 1) / n;
if r == 1, uc = uf; return; end
uc = squeeze(mean(mean(reshape(uf, r, n, r, n), 1), 3));
end
