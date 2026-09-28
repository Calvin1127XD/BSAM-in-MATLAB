%runDdGeometries Diffuse-domain geometry sweep: circle, ellipse, regular
%   pentagon (exact signed distance; 5 corners with gradient jumps along
%   the bisectors), and non-self-intersecting polar curves (3-petal
%   flower, peanut, and the paper's 5-lobe star).
%
%   G1: Neumann form on every geometry (manufactured u = r^2/4 with the
%       analytic extension g = grad u . n_psi, so by the exact-reproduction
%       lemma the model error vanishes and the runs measure layer
%       discretization + solver robustness across geometries).
%       Asserts: conservation to roundoff, geometric V-cycles, refinement
%       localized at Gamma, error decreasing with eps (eps/h fixed).
%   G2: Dirichlet penalty on the PENTAGON (u = r^2/4, g = u on Gamma):
%       a genuine model-error test on a domain with corners; reports the
%       eps-convergence rates of the corner-limited DDM asymptotics.
%
%   Writes ../output/results_geometries.mat and per-geometry mesh figures.
clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'src', 'ddAddpaths.m'));
outdir = fullfile(here, '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end
warning('off', 'eb_tag_cluster:dropthin');
R = struct();
PASS = true;
t0all = tic;

box = [-2.5 2.5 -2.5 2.5];
uex = @(x,y) (x.^2 + y.^2)/4;
fex = @(x,y) (x.^2 + y.^2)/4 - 1;     % f = u - lap u (beta = 1)

% {name, geometry, minimum eps-decay rate for the Neumann run}
% Smooth boundaries: expect ~1.4-1.6 (layer discretization, eps/h fixed).
% Polygon: the level-set normal jumps across the corner bisectors, so the
% extended flux data g = grad u . n_psi is DISCONTINUOUS along five lines
% inside the band; the corner neighborhoods converge at reduced order in
% L2 (L_inf still decays at ~1.7).  This is intrinsic to flux boundary
% data on corners; the Dirichlet penalty form (G2), which never uses the
% normal, recovers rates 1.5-1.8 on the same pentagon.
geoms = { ...
  'circle',   ddGeometry('circle', 1.0),        0.8; ...
  'ellipse',  ddGeometry('ellipse', 1.4, 0.8),  0.8; ...
  'pentagon', ddGeometry('pentagon', 1.2),      0.4; ...
  'flower3',  ddGeometry('flower3'),            0.8; ...
  'peanut',   ddGeometry('peanut'),             0.8; ...
  'star5',    ddGeometry('star', 0.5),          0.8};

epss = [0.2 0.1];
Ns   = [256 512];

fprintf('\n=== G1: Neumann DDM across geometries =============================\n');
G1 = nan(size(geoms,1), numel(epss), 6);
for gidx = 1:size(geoms, 1)
  gname = geoms{gidx, 1};
  geom = geoms{gidx, 2};
  gfun = @(X,Y) neumannG(geom, X, Y);
  data = struct('f', fex, 'g', gfun, 'C', 1, 'D', 1, 'exact', uex);
  for i = 1:numel(epss)
    ep = epss(i); N = Ns(i);
    sol = ddSolve(geom, ep, box, 'neumann', data, ...
      'BaseCells', (N/8)*[1 1], 'MaxLevels', 4, 'TagWidth', 2.5, ...
      'TagBuffer', 2, 'Verbose', 0, 'RelTol', 1e-10, 'MaxVCycles', 30);
    st = sol.stats.stages{end};
    loc = localizationOk(sol.H, geom, ep, box);
    G1(gidx, i, :) = [sol.errInside.l2, sol.errInside.linf, ...
                      sol.mass.rel, st.avgFactor, sol.err.dof, loc];
    fprintf('  %-9s eps=%4.2f: L2=%.3e Linf=%.3e mass=%.0e factor=%.3f dof=%6d loc=%d\n', ...
      gname, ep, sol.errInside.l2, sol.errInside.linf, sol.mass.rel, ...
      st.avgFactor, sol.err.dof, loc);
    if sol.mass.rel > 1e-11
      fprintf('  ** G1 FAIL (%s: conservation)\n', gname); PASS = false;
    end
    if st.avgFactor > 0.35
      fprintf('  ** G1 FAIL (%s: V-cycle factor)\n', gname); PASS = false;
    end
    if ~loc
      fprintf('  ** G1 FAIL (%s: refinement not localized)\n', gname); PASS = false;
    end
    if i == numel(epss)
      ebPlot(sol.H, 'mesh', ...
        fullfile(outdir, sprintf('g1_%s_mesh.png', gname)), ...
        sprintf('%s: adaptive 4-level mesh, eps=%.2f', geom.name, ep));
    end
  end
  rate = log2(G1(gidx, 1, 1) / G1(gidx, 2, 1));
  ratei = log2(G1(gidx, 1, 2) / G1(gidx, 2, 2));
  fprintf('  %-9s decay rates (eps 0.2 -> 0.1, eps/h fixed): L2 %.2f, Linf %.2f\n', ...
    gname, rate, ratei);
  if rate < geoms{gidx, 3}
    fprintf('  ** G1 FAIL (%s: error not decreasing with eps)\n', gname);
    PASS = false;
  end
end
R.G1 = G1;
R.G1geoms = geoms(:, 1);

% ======================================================================
fprintf('\n=== G2: pentagon, Dirichlet penalty (corners + model error) =======\n');
geomP = ddGeometry('pentagon', 1.2);
dataP = struct('f', fex, 'g', uex, 'C', 1, 'D', 1, 'exact', uex);
epsP = [0.4 0.2 0.1 0.05];
NP   = [128 256 512 1024];
G2 = nan(numel(epsP), 5);
for i = 1:numel(epsP)
  sol = ddSolve(geomP, epsP(i), box, 'dirichlet', dataP, ...
    'BaseCells', (NP(i)/8)*[1 1], 'MaxLevels', 4, 'TagWidth', 2.5, ...
    'TagBuffer', 2, 'Mu', 1, 'Verbose', 0, 'RelTol', 1e-10, 'MaxVCycles', 40);
  st = sol.stats.stages{end};
  G2(i,:) = [epsP(i), sol.errInside.l2, sol.errInside.linf, ...
             st.avgFactor, sol.err.dof];
  fprintf('  eps=%5.3f  L2=%.3e  Linf=%.3e  factor=%.3f  dof=%d\n', ...
    epsP(i), G2(i,2), G2(i,3), G2(i,4), G2(i,5));
  if i == 3
    ebPlot(sol.H, 'mesh', fullfile(outdir, 'g2_pentagon_mesh.png'), ...
      sprintf('Pentagon (Dirichlet penalty): adaptive mesh, eps=%.2f', epsP(i)));
    ebPlot(sol.H, 'solution', fullfile(outdir, 'g2_pentagon_solution.png'), ...
      sprintf('Pentagon: diffuse-domain solution, eps=%.2f', epsP(i)));
  end
end
r2 = log2(G2(1:end-1,2) ./ G2(2:end,2));
fprintf('  pentagon Dirichlet L2 rates in eps: %s\n', mat2str(r2.', 2));
R.G2 = G2; R.G2rates = r2;
if any(diff(G2(:,2)) >= 0)
  fprintf('  ** G2 FAIL (pentagon error not monotone in eps)\n'); PASS = false;
end

% ======================================================================
R.elapsed = toc(t0all);
save(fullfile(outdir, 'results_geometries.mat'), 'R');
fprintf('\nTotal geometry-sweep time: %.1f s\n', R.elapsed);
if PASS
  fprintf('ALL GEOMETRY TESTS PASSED\n');
else
  error('run_dd_geometries: at least one test failed');
end

% ----------------------------------------------------------------------
function g = neumannG(geom, X, Y)
% boundary data g = grad u . n for u = r^2/4, n = grad psi / |grad psi|
gg = geom.gpsi(X, Y);
g = (X.*gg{1} + Y.*gg{2}) ./ (2 * max(geom.gpsimag(X, Y), 1e-12));
end

% ----------------------------------------------------------------------
function ok = localizationOk(H, geom, ep, box)
% every refined patch must touch the tag band |psi| < 2.5 eps, with
% cell-quantization slack converted to psi units by the local |grad psi|
ok = true;
for k = H.rootIdx+1:H.nlev
  lev = H.lev{k};
  for p = 1:lev.np
    b = lev.box(p, :);
    xc = box(1) + ((b(1):b(2)) - 0.5) * lev.h(1);
    yc = box(3) + ((b(3):b(4)) - 0.5) * lev.h(2);
    [XC, YC] = ndgrid(xc, yc);
    gm = max(geom.gpsimag(XC(:), YC(:)));
    if min(abs(geom.psi(XC(:), YC(:)))) > 2.5*ep + 8 * (2*max(lev.h)) * gm
      ok = false;
      return;
    end
  end
end
end
