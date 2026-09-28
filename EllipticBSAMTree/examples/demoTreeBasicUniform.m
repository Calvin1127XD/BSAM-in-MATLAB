%demoBasicUniform Minimal example: variable D(x,y), C(x,y), Dirichlet
%   data, single-level (uniform) solve with the geometric multigrid path.
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

% Manufactured solution: pick u, D, C, then build f = -div(D grad u) + C u
% by the product rule so the discrete error can be measured exactly.
uex = @(x,y) sin(pi*x).*cos(pi*y);
D = @(x,y) 1 + 0.5*sin(2*pi*x).*cos(pi*y);
C = @(x,y) 1 + x.*y;
Dx = @(x,y) pi*cos(2*pi*x).*cos(pi*y);
Dy = @(x,y) -0.5*pi*sin(2*pi*x).*sin(pi*y);
ux = @(x,y) pi*cos(pi*x).*cos(pi*y);
uy = @(x,y) -pi*sin(pi*x).*sin(pi*y);
f = @(x,y) -(Dx(x,y).*ux(x,y) + Dy(x,y).*uy(x,y) - 2*pi^2*D(x,y).*uex(x,y)) ...
         + C(x,y).*uex(x,y);

prob = ebtProblem([0 1 0 1], D, C, f, {'dirichlet', uex}, uex);
sol = ebtSolve(prob, 'BaseCells', [128 128], 'Verbose', 2, 'RelTol', 1e-11);

fprintf('\nerror vs exact: L_inf = %.3e, L2 = %.3e\n', sol.err.linf, sol.err.l2);
