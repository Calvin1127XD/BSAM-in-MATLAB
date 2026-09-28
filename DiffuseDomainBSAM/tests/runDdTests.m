%runDdTests DiffuseDomainBSAM verification suite.
%   D1: paper sec. 6.4 star-shaped domain, Neumann BC via |grad phi| g:
%       eps -> 0 convergence on uniform vs adaptive (QIF) vs adaptive
%       without mass corrections; mass conservation; mesh localization.
%   D2: Dirichlet penalty formulation on a disc (exact solution).
%   D3: efficiency: DOFs / wall time, adaptive vs uniform.
%   Writes ../output/results_dd.mat and figures.
clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'src', 'ddAddpaths.m'));
outdir = fullfile(here, '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end
R = struct();
PASS = true;
t0all = tic;
warning('off', 'eb_tag_cluster:dropthin');   % benign band-edge slivers

% ======================================================================
% D1: star domain (paper eq. 6.19), sharp problem  u - lap u = f  in
% Omega_1, grad u . n = g on Gamma; exact u = (x^2+y^2)/4, so
% f = r^2/4 - 1 and g = (x psix + y psiy) / (2 |grad psi|).
fprintf('\n=== D1: star-shaped domain, Neumann DDM (paper Table 6) ==========\n');
geom = ddGeometry('star', 0.5);
box = [-2.5 2.5 -2.5 2.5];
uex = @(x,y) (x.^2 + y.^2)/4;
fex = @(x,y) (x.^2 + y.^2)/4 - 1;
gfun = @(X,Y) gfunStar(geom, X, Y);   % g = grad u . n (closure over geom)
data = struct('f', fex, 'g', gfun, 'C', 1, 'D', 1, 'exact', uex);

epss = [0.8 0.4 0.2 0.1 0.05];
Nuni = [128 256 512 1024 2048];
D1 = nan(numel(epss), 12);
for i = 1:numel(epss)
  ep = epss(i); N = Nuni(i);

  % ---- uniform reference
  t = tic;
  solU = ddSolve(geom, ep, box, 'neumann', data, ...
    'BaseCells', N*[1 1], 'MaxLevels', 1, 'Verbose', 0, ...
    'RelTol', 1e-10, 'MaxVCycles', 30);
  tU = toc(t);

  % ---- adaptive: root N/8, 4 levels (effective N), QIF + corrections
  t = tic;
  solA = ddSolve(geom, ep, box, 'neumann', data, ...
    'BaseCells', (N/8)*[1 1], 'MaxLevels', 4, 'TagWidth', 2.5, ...
    'TagBuffer', 2, 'Verbose', 0, 'RelTol', 1e-10, ...
    'MaxVCycles', 30);
  tA = toc(t);

  % ---- adaptive without the flux-balance corrections ("QI")
  solQ = ddSolve(geom, ep, box, 'neumann', data, ...
    'BaseCells', (N/8)*[1 1], 'MaxLevels', 4, 'TagWidth', 2.5, ...
    'TagBuffer', 2, 'MassCorrection', false, ...
    'Verbose', 0, 'RelTol', 1e-10, 'MaxVCycles', 30);

  stA = solA.stats.stages{end};
  D1(i,:) = [ep, ...
    solU.errInside.l2, solU.errInside.linf, ...
    solA.errInside.l2, solA.errInside.linf, ...
    solQ.errInside.l2, solQ.errInside.linf, ...
    solU.mass.rel, solA.mass.rel, solQ.mass.rel, ...
    stA.avgFactor, tA / max(tU, eps)];
  fprintf(['  eps=%4.2f  uni(L2)=%.3e  amr(L2)=%.3e  noCorr(L2)=%.3e | ' ...
           'mass: %0.0e / %0.0e / %0.0e | amr factor %.3f | t_amr/t_uni=%.2f\n'], ...
    ep, D1(i,2), D1(i,4), D1(i,6), D1(i,8), D1(i,9), D1(i,10), ...
    D1(i,11), D1(i,12));

  R.D1cells(i,:) = [solU.err.dof, solA.err.dof];
  R.D1times(i,:) = [tU, tA];

  if i == 4
    ebPlot(solA.H, 'mesh', fullfile(outdir, 'd1_star_mesh.png'), ...
      sprintf('Star domain: adaptive 4-level mesh, eps=%.2f', ep));
    ebPlot(solA.H, 'solution', fullfile(outdir, 'd1_star_solution.png'), ...
      sprintf('Diffuse-domain solution u_eps, eps=%.2f', ep));
  end

  % ---- mesh localization: every refined patch must intersect the band
  Hh = solA.H; okloc = true;
  for k = Hh.rootIdx+1:Hh.nlev
    lev = Hh.lev{k};
    for p = 1:lev.np
      b = lev.box(p, :);
      xc = box(1) + ((b(1):b(2)) - 0.5) * lev.h(1);
      yc = box(3) + ((b(3):b(4)) - 0.5) * lev.h(2);
      [XC, YC] = ndgrid(xc, yc);
      % every refined patch must touch the tag band |psi| < 2.5 eps.
      % Slack: tags live on PARENT cells (spacing 2h) and boxes carry
      % buffer (2) + min-width growth (<=4) + merge slack (~2) cells;
      % convert cells to psi-units with the local |grad psi| (up to ~7
      % for the star's lobes).
      gm = max(geom.gpsimag(XC(:), YC(:)));
      if min(abs(geom.psi(XC(:), YC(:)))) > 2.5*ep + 8 * (2*max(lev.h)) * gm
        okloc = false;
      end
    end
  end
  if ~okloc
    fprintf('  ** D1 FAIL: refinement not localized at the interface\n');
    PASS = false;
  end
end
rU = log2(D1(1:end-1,2) ./ D1(2:end,2));
rA = log2(D1(1:end-1,4) ./ D1(2:end,4));
rQ = log2(D1(1:end-1,6) ./ D1(2:end,6));
fprintf('  L2 rates in eps:  uniform %s | adaptive %s | no-corr %s\n', ...
  mat2str(rU.',2), mat2str(rA.',2), mat2str(rQ.',2));
R.D1 = D1; R.D1rates = [rU rA rQ];
% the star level set is NOT a distance function (|grad psi| varies along
% Gamma), so the tanh profile is asymmetric in the normal coordinate and
% the DDM model error sits between O(eps) and O(eps^2); see docs.  The
% clean O(eps^2) assertion is made on the circle (D1b).  Here we assert:
% (a) the adaptive solution matches the uniform one (AMR adds no error),
% (b) rates stay clearly superlinear, (c) conservation to roundoff.
if any(rA < 1.2), fprintf('  ** D1 FAIL (adaptive eps-rate < 1.2)\n'); PASS = false; end
amrclose = abs(D1(:,4) - D1(:,2)) ./ D1(:,2);
fprintf('  |amr - uniform| / uniform (L2): %s\n', mat2str(amrclose.', 2));
if any(amrclose > 0.10)
  fprintf('  ** D1 FAIL (adaptive deviates from uniform by > 10%%)\n'); PASS = false;
end
if any(D1(:,9) > 1e-12), fprintf('  ** D1 FAIL (conservative mass defect)\n'); PASS = false; end
if max(D1(:,10)) < 1e4 * max(D1(:,9))
  fprintf('  note: no-correction mass defect not >> conservative one\n');
end

% ======================================================================
% D1b: circle (true signed distance): isolate the DDM MODEL error from
% the layer discretization error by Richardson extrapolation in h at
% each eps (E(h) = M(eps) + c h^2 within the asymptotic regime), then
% rate the extrapolated model error M(eps): expect O(eps^2)
% (Lervag & Lowengrub).
%
% NOTE on extensions: with the analytic extension g_ext = grad u . n_psi
% the identity |grad phi| g_ext = -grad phi . grad u holds pointwise, so
% u = r^2/4 solves the CONTINUOUS diffuse-domain problem exactly and the
% model error is identically zero (see docs/ddm_analysis.pdf, exact-
% reproduction lemma) -- D1 above therefore measures pure discretization
% error.  To expose a genuine O(eps^2) model error we use here the
% canonical CONSTANT-NORMAL extension g_ext = 1/2 (the value of
% grad u . n on Gamma).
fprintf('\n=== D1b: circle, Neumann DDM: model error via extrapolation ======\n');
geoc = ddGeometry('circle', 1.0);
boxc = [-2 2 -2 2];
gfunc = @(X,Y) 0.5 * ones(size(X));           % constant-normal extension
datac = struct('f', fex, 'g', gfunc, 'C', 1, 'D', 1, 'exact', uex);
epsc = [0.4 0.2 0.1];
Nc   = [128 256 512];
D1b = nan(numel(epsc), 5);
for i = 1:numel(epsc)
  E = nan(1, 2);
  for j = 1:2
    N = Nc(i) * 2^(j-1);
    solC = ddSolve(geoc, epsc(i), boxc, 'neumann', datac, ...
      'BaseCells', (N/8)*[1 1], 'MaxLevels', 4, 'TagWidth', 2.5, ...
      'TagBuffer', 2, 'Verbose', 0, ...
      'RelTol', 1e-10, 'MaxVCycles', 30);
    E(j) = solC.errInside.l2;
  end
  M = max((4*E(2) - E(1)) / 3, 1e-14);     % h -> 0 limit of the error
  D1b(i,:) = [epsc(i), E(1), E(2), M, log2(E(1)/E(2))];
  fprintf('  eps=%4.2f  E(eps/h=12.8)=%.3e  E(25.6)=%.3e  h-rate=%.2f  model M=%.3e\n', ...
    epsc(i), E(1), E(2), D1b(i,5), M);
end
r1b = log2(D1b(1:end-1,4) ./ D1b(2:end,4));
fprintf('  model-error rates in eps: %s\n', mat2str(r1b.', 2));
R.D1b = D1b; R.D1brates = r1b;
if any(r1b < 1.5), fprintf('  ** D1b FAIL (model eps-rate < 1.5)\n'); PASS = false; end

% ======================================================================
fprintf('\n=== D2: disc, Dirichlet penalty DDM ==============================\n');
% -lap u = 4 in unit disc, u = 0 on Gamma; exact u = 1 - r^2.
geo2 = ddGeometry('circle', 1.0);
box2 = [-1.5 1.5 -1.5 1.5];
uex2 = @(x,y) 1 - x.^2 - y.^2;
data2 = struct('f', 4, 'g', 0, 'C', 0, 'D', 1, 'exact', uex2);
eps2 = [0.4 0.2 0.1 0.05];
N2   = [96 192 384 768];
D2 = nan(numel(eps2), 5);
for i = 1:numel(eps2)
  ep = eps2(i);
  sol = ddSolve(geo2, ep, box2, 'dirichlet', data2, ...
    'BaseCells', (N2(i)/8)*[1 1], 'MaxLevels', 4, 'TagWidth', 2.5, ...
    'TagBuffer', 2, 'Mu', 1, 'Verbose', 0, 'RelTol', 1e-10, 'MaxVCycles', 40);
  st = sol.stats.stages{end};
  D2(i,:) = [ep, sol.errInside.l2, sol.errInside.linf, st.avgFactor, ...
             sol.err.dof];
  fprintf('  eps=%5.3f  L2=%.3e  Linf=%.3e  factor=%.3f  dof=%d\n', ...
    ep, D2(i,2), D2(i,3), D2(i,4), D2(i,5));
  if i == 3
    ebPlot(sol.H, 'mesh', fullfile(outdir, 'd2_disc_mesh.png'), ...
      sprintf('Disc (Dirichlet penalty): adaptive mesh, eps=%.2f', ep));
  end
end
r2 = log2(D2(1:end-1,2) ./ D2(2:end,2));
fprintf('  L2 rates in eps: %s\n', mat2str(r2.', 2));
R.D2 = D2; R.D2rates = r2;
if any(r2 < 0.8), fprintf('  ** D2 FAIL (penalty eps-rate < 0.8)\n'); PASS = false; end

% ======================================================================
fprintf('\n=== D3: efficiency, adaptive vs uniform ===========================\n');
fprintf('  eps    cells(uni)  cells(amr)  ratio   t_uni    t_amr\n');
for i = 1:numel(epss)
  fprintf('  %4.2f   %9d  %9d   %4.1fx   %6.2fs  %6.2fs\n', epss(i), ...
    R.D1cells(i,1), R.D1cells(i,2), R.D1cells(i,1)/R.D1cells(i,2), ...
    R.D1times(i,1), R.D1times(i,2));
end

% V-cycle factor across eps (degenerate-coefficient h-independence)
fprintf('  adaptive V-cycle factors across eps: %s\n', mat2str(D1(:,11).', 2));
if max(D1(:,11)) > 0.35
  fprintf('  ** D3 FAIL (V-cycle factor degraded on DDM coefficients)\n');
  PASS = false;
end

% ======================================================================
R.elapsed = toc(t0all);
save(fullfile(outdir, 'results_dd.mat'), 'R');
fprintf('\nTotal DD suite time: %.1f s\n', R.elapsed);
if PASS
  fprintf('ALL FOLDER-2 TESTS PASSED\n');
else
  error('run_dd_tests: at least one test failed');
end

% ---------------------------------------------------------------------
function g = gfunStar(geom, X, Y)
% boundary data g = grad u . n for u = r^2/4, n = grad psi / |grad psi|
gg = geom.gpsi(X, Y);
g = (X.*gg{1} + Y.*gg{2}) ./ (2 * max(geom.gpsimag(X, Y), 1e-12));
end
