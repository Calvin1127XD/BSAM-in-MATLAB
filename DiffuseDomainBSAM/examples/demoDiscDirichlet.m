%demoDiscDirichlet Diffuse-domain solve of  -lap u = 4  on the unit
%   disc with homogeneous Dirichlet boundary conditions imposed by the
%   penalty formulation; exact solution u = 1 - r^2.
clear; clc;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'ddAddpaths.m'));
outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end

geom = ddGeometry('circle', 1.0);
box = [-1.5 1.5 -1.5 1.5];
epsi = 0.005;
uex = @(x,y) 1 - x.^2 - y.^2;
data = struct('f', 4, 'g', 0, 'C', 0, 'D', 1, 'exact', uex);

sol = ddSolve(geom, epsi, box, 'dirichlet', data, ...
  'BaseCells', [64 64], 'MaxLevels', 4, 'TagWidth', 2.5, 'TagBuffer', 2, ...
  'Mu', 1, 'Verbose', 1, 'RelTol', 1e-10);

fprintf('\ninterior error: L2 = %.3e, Linf = %.3e (eps = %.4f)\n', ...
  sol.errInside.l2, sol.errInside.linf, epsi);

ebPlot(sol.H, 'mesh', fullfile(outdir, 'demo_disc_mesh.png'), ...
  'Adaptive mesh, disc with Dirichlet penalty');
ebPlot(sol.H, 'solution', fullfile(outdir, 'demo_disc_solution.png'), ...
  'Diffuse-domain solution');
fprintf('figures written to %s\n', outdir);
