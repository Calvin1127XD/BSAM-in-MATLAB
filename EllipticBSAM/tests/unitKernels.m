%unitKernels Unit tests for the numerical kernels (run before the suite).
%   1. QIF/LIF ghost interpolation reproduces global quadratics/linears.
%   2. Operator is 2nd-order on a manufactured solution.
%   3. Prolongation/restriction sanity.
%   4. Uniform-grid V-cycle (no AMR) matches the direct solve.

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
rng(7);
failures = 0;

% ---------- 1. ghost interpolation polynomial reproduction --------------
% Fine patch: mx x my cells of size h starting at (x0,y0); parent window
% cells of size 2h.  Fill interiors with u, check ghosts vs u exactly.
for method = ["qif", "lif"]
  if method == "qif"
    u = @(x, y) 1.5 + 0.3*x - 0.7*y + 0.25*x.^2 - 0.4*x.*y + 0.6*y.^2;
  else
    u = @(x, y) 1.5 + 0.3*x - 0.7*y;
  end
  mx = 8; my = 12; h = 0.05; x0 = 0.3; y0 = -0.2;
  xg = x0 + ((1:mx) - 0.5)*h;  yg = y0 + ((1:my) - 0.5)*h;
  [XF, YF] = ndgrid(x0 + ((0:mx+1) - 0.5)*h, y0 + ((0:my+1) - 0.5)*h);
  Qex = u(XF, YF);
  Q = Qex;
  Q(1, :) = 99; Q(end, :) = 99; Q(:, 1) = 99; Q(:, end) = 99;  % poison ghosts
  cmx = mx/2; cmy = my/2;
  [XC, YC] = ndgrid(x0 + ((0:cmx+1) - 1)*2*h + h, y0 + ((0:cmy+1) - 1)*2*h + h);
  QC = u(XC, YC);
  Q2 = ebInterpCf(Q, QC, true(1,4), true(1,4), char(method));
  ghostmask = false(size(Q)); ghostmask([1 end], :) = true; ghostmask(:, [1 end]) = true;
  errg = max(abs(Q2(ghostmask) - Qex(ghostmask)));
  fprintf('interp %s ghost reproduction error: %.3e\n', method, errg);
  if errg > 1e-12, failures = failures + 1; fprintf('  ** FAIL\n'); end
end

% ---------- 2. operator truncation -----------------------------------
K = ebKernels();
uex = @(x,y) sin(pi*x).*cos(pi*y);
Dfun = @(x,y) 1 + 0.5*sin(2*pi*x).*cos(pi*y);
Cfun = @(x,y) 1 + x.*y;
% f = -div(D grad u) + C u with D variable:
Dx = @(x,y) pi*cos(2*pi*x).*cos(pi*y);
Dy = @(x,y) -0.5*pi*sin(2*pi*x).*sin(pi*y);
ux = @(x,y) pi*cos(pi*x).*cos(pi*y);
uy = @(x,y) -pi*sin(pi*x).*sin(pi*y);
lap = @(x,y) -2*pi^2*uex(x,y);
ffun = @(x,y) -(Dx(x,y).*ux(x,y) + Dy(x,y).*uy(x,y) + Dfun(x,y).*lap(x,y)) ...
            + Cfun(x,y).*uex(x,y);
errs = [];
for n = [32 64 128]
  h = 1/n;
  [XF, YC] = ndgrid((0:n)*h, ((1:n)-0.5)*h);
  DeW = Dfun(XF, YC);
  [XC, YF] = ndgrid(((1:n)-0.5)*h, (0:n)*h);
  DnS = Dfun(XC, YF);
  [X, Y] = ndgrid(((0:n+1)-0.5)*h, ((0:n+1)-0.5)*h);
  Q = uex(X, Y);
  CCa = Cfun(X(2:end-1,2:end-1), Y(2:end-1,2:end-1));
  L = K.op(Q, DeW, DnS, CCa, h, h);
  F = ffun(X(2:end-1,2:end-1), Y(2:end-1,2:end-1));
  errs(end+1) = max(abs(L(:) - F(:))); %#ok<AGROW>
end
rates = log2(errs(1:end-1) ./ errs(2:end));
fprintf('operator truncation errors: %s rates: %s\n', mat2str(errs, 3), mat2str(rates, 3));
if any(rates < 1.8), failures = failures + 1; fprintf('  ** FAIL\n'); end

% ---------- 3. prolong/restrict --------------------------------------
% bilinear prolongation is exact on affine data, and conservative
% restriction of an affine field returns its cell averages, so
% restrict(prolong(W)) == W for affine W (catches any index swap).
[II, JJ] = ndgrid(0:5, 0:7);
W = 0.3 + 0.2*II - 0.45*JJ;           % padded coarse (4x6 interior)
Fp = K.prolong(W);
R = K.restrict2(Fp);
errpr = max(max(abs(R - W(2:end-1, 2:end-1))));
fprintf('restrict(prolong(affine)) - affine: %.3e\n', errpr);
if errpr > 1e-13, failures = failures + 1; fprintf('  ** FAIL\n'); end

% ---------- 4. uniform V-cycle vs direct -----------------------------
prob = ebProblem([0 1 0 1], Dfun, Cfun, ffun, {'dirichlet', uex}, uex);
opts = ebOptions('BaseCells', [64 64], 'MaxLevels', 1, 'Verbose', 0, ...
                  'RelTol', 1e-11);
sol = ebSolve(prob, opts);
fprintf('uniform 64^2: V-cycles = %d, rel res = %.2e, err_inf = %.3e\n', ...
  sol.stats.stages{1}.vcycles, sol.stats.stages{1}.finalRel, sol.err.linf);
if sol.stats.stages{1}.finalRel > 1e-11 || sol.err.linf > 5e-3
  failures = failures + 1; fprintf('  ** FAIL\n');
end

% ----------------------------------------------------------------------
if failures == 0
  fprintf('UNIT TESTS PASSED\n');
else
  error('unit_kernels: %d failure(s)', failures);
end
