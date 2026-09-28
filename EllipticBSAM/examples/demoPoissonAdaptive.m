%demoPoissonAdaptive Adaptive solve of a variable-coefficient elliptic
%   problem with three sharp interior bumps; the undivided-Laplacian
%   indicator places patches automatically.
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end

% manufactured: u = sum of Gaussians, f = -div(D grad u) + C u with D=C=1
s = 0.02;
xi = [0.30 0.66 0.49]; yi = [0.34 0.42 0.72]; Ai = [1 0.8 -0.9];
g  = @(x,y,a,b) exp(-((x-a).^2 + (y-b).^2)/(2*s^2));
lg = @(x,y,a,b) g(x,y,a,b).*(((x-a).^2 + (y-b).^2)/s^4 - 2/s^2);
uex = @(x,y) Ai(1)*g(x,y,xi(1),yi(1)) + Ai(2)*g(x,y,xi(2),yi(2)) + Ai(3)*g(x,y,xi(3),yi(3));
lap = @(x,y) Ai(1)*lg(x,y,xi(1),yi(1)) + Ai(2)*lg(x,y,xi(2),yi(2)) + Ai(3)*lg(x,y,xi(3),yi(3));
f = @(x,y) uex(x,y) - lap(x,y);

prob = ebProblem([0 1 0 1], 1, 1, f, {'dirichlet', uex}, uex);
% The undivided Laplacian scales like h^2, so divide the threshold by 4
% per level to keep tagging the same physical feature after refinement.
opts = ebOptions('BaseCells', [64 64], 'MaxLevels', 4, ...
  'Indicator', 'ulap', 'TagThreshold', @(L) 1e-2/4^(L-1), 'TagBuffer', 2, ...
  'Verbose', 1, 'RelTol', 1e-10);

sol = ebSolve(prob, opts);
fprintf('\nerror: L_inf = %.3e, L2 = %.3e, dof = %d\n', ...
  sol.err.linf, sol.err.l2, sol.err.dof);
fprintf('mass conservation defect (relative): %.2e\n', sol.mass.rel);

ebPlot(sol.H, 'mesh', fullfile(outdir, 'demo_adaptive_mesh.png'), ...
  'Adaptive mesh (4 levels)');
ebPlot(sol.H, 'solution', fullfile(outdir, 'demo_adaptive_solution.png'), ...
  'Composite solution');
fprintf('figures written to %s\n', outdir);
