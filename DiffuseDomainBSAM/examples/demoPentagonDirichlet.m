%demoPentagonDirichlet Diffuse-domain solve of the Helmholtz problem
%   u - lap u = f on a regular pentagon (exact signed-distance level set,
%   five corners) with Dirichlet data u = g on the boundary imposed by
%   the penalty formulation; manufactured u = r^2/4.
clear; clc;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'ddAddpaths.m'));
outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end

geom = ddGeometry('pentagon', 1.2);     % circumradius 1.2, point up
box = [-2.5 2.5 -2.5 2.5];
epsi = 0.01;
uex = @(x,y) (x.^2 + y.^2)/4;
fex = @(x,y) (x.^2 + y.^2)/4 - 1;        % f = u - lap u
data = struct('f', fex, 'g', uex, 'C', 1, 'D', 1, 'exact', uex);

sol = ddSolve(geom, epsi, box, 'dirichlet', data, ...
  'BaseCells', [64 64], 'MaxLevels', 4, 'TagWidth', 2.5, 'TagBuffer', 2, ...
  'Mu', 1, 'Verbose', 1, 'RelTol', 1e-10);

fprintf('\ninterior error: L2 = %.3e, Linf = %.3e (eps = %.4f)\n', ...
  sol.errInside.l2, sol.errInside.linf, epsi);
fprintf('mass conservation defect (relative): %.2e\n', sol.mass.rel);

ebPlot(sol.H, 'mesh', fullfile(outdir, 'demo_pentagon_mesh.png'), ...
  'Pentagon: adaptive mesh tracking the corners');
ebPlot(sol.H, 'solution', fullfile(outdir, 'demo_pentagon_solution.png'), ...
  'Pentagon: diffuse-domain solution');
fprintf('figures written to %s\n', outdir);
