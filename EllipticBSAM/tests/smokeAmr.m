%smokeAmr First integration test of the AMR machinery: paper-style
%   box-in-box hierarchy on the periodic Gaussian (paper eq. 5.1) with
%   exact homogeneous Neumann data, C = 1 (Helmholtz beta = 1).
clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

sg = 0.25;
gp  = @(t) -pi*sin(4*pi*t)/sg^2;          % g'(t), g = -sin(2 pi t)^2/(2 sg^2)
gpp = @(t) -4*pi^2*cos(4*pi*t)/sg^2;
uex = @(x,y) exp(-(sin(2*pi*x).^2 + sin(2*pi*y).^2) / (2*sg^2));
ffun = @(x,y) uex(x,y) .* (1 - (gpp(x) + gp(x).^2 + gpp(y) + gp(y).^2));

prob = ebProblem([-0.25 0.25 -0.25 0.25], 1, 1, ffun, {'neumann', 0}, uex);

% fixed box-in-box: level 2 over (-0.1875, 0.1875)^2, level 3 over
% (-0.125, 0.125)^2 (paper section 5.1)
tagfn = @(X, Y, L) ...
  (L == 1) .* (max(abs(X), abs(Y)) < 0.1875) + ...
  (L == 2) .* (max(abs(X), abs(Y)) < 0.125);

for m0 = [32 64]
  for nlv = [2 3]
    opts = ebOptions('BaseCells', [m0 m0], 'MaxLevels', nlv, ...
      'Indicator', 'fun', 'TagFunction', tagfn, 'TagBuffer', 0, ...
      'Verbose', 1, 'RelTol', 1e-10, 'MaxVCycles', 40);
    sol = ebSolve(prob, opts);
    st = sol.stats.stages{end};
    fprintf('m0=%3d  levels=%d : err_inf=%.4e  err_l2=%.4e  factor=%.3f  mass=%.2e\n', ...
      m0, nlv, sol.err.linf, sol.err.l2, st.avgFactor, sol.mass.rel);
  end
end
