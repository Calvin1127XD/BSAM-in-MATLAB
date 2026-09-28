%demoStarNeumann Diffuse-domain solve of  u - lap u = f  on the
%   star-shaped (flower) domain of Feng et al. (2018), sec. 6.4, with
%   inhomogeneous Neumann data imposed through the |grad phi| g source
%   term, and block-structured adaptive refinement of the diffuse layer.
clear; clc;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'ddAddpaths.m'));
outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end

geom = ddGeometry('star', 0.5);
box = [-2.5 2.5 -2.5 2.5];
epsi = 0.1;

uex = @(x,y) (x.^2 + y.^2)/4;          % sharp exact solution
fex = @(x,y) (x.^2 + y.^2)/4 - 1;      % f = u - lap u
gfun = @(X,Y) starG(geom, X, Y);      % g = grad u . n on Gamma (extended)
data = struct('f', fex, 'g', gfun, 'C', 1, 'D', 1, 'exact', uex);

sol = ddSolve(geom, epsi, box, 'neumann', data, ...
  'BaseCells', [64 64], 'MaxLevels', 4, 'TagWidth', 2.5, 'TagBuffer', 2, ...
  'Verbose', 1, 'RelTol', 1e-10);

fprintf('\ninterior error: L2 = %.3e, Linf = %.3e (eps = %.4f)\n', ...
  sol.errInside.l2, sol.errInside.linf, epsi);
fprintf('mass conservation defect (relative): %.2e\n', sol.mass.rel);

ebPlot(sol.H, 'mesh', fullfile(outdir, 'demo_star_mesh.png'), ...
  'Adaptive mesh tracking the diffuse interface');
ebPlot(sol.H, 'solution', fullfile(outdir, 'demo_star_solution.png'), ...
  'Diffuse-domain solution u_\epsilon');
fprintf('figures written to %s\n', outdir);

function g = starG(geom, X, Y)
gg = geom.gpsi(X, Y);
g = (X.*gg{1} + Y.*gg{2}) ./ (2 * max(geom.gpsimag(X, Y), 1e-12));
end
