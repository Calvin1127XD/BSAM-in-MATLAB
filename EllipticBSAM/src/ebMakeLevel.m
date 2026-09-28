function lev = ebMakeLevel(H, k, n, boxes)
%ebMakeLevel Create a level struct with patches [ilo ihi jlo jhi] (global
%   cell indices on the level-k grid of n = [nx ny] cells).  Evaluates the
%   analytic coefficients/source, builds red-black masks (global parity),
%   classifies physical faces and evaluates boundary data.

prob = H.prob;
dom = prob.domain;
lev.n = n;
lev.h = [(dom(2)-dom(1))/n(1), (dom(4)-dom(3))/n(2)];
lev.np = size(boxes, 1);
lev.box = boxes;
lev.sourceShift = 0; % constant removed from the originally sampled source

for p = 1:lev.np
  b = boxes(p, :);
  if any(b(1:2) < 1 | b(1:2) > n(1)) || any(b(3:4) < 1 | b(3:4) > n(2)) ...
      || b(1) > b(2) || b(3) > b(4)
    error('eb_make_level:box', 'invalid patch box.');
  end
  mx = b(2) - b(1) + 1;
  my = b(4) - b(3) + 1;
  lev.Q{p} = zeros(mx+2, my+2);

  % Physical coordinates: xf/yf at face positions, xc/yc at cell centers.
  % D is sampled at face midpoints (where the fluxes live), C and f at
  % cell centers.
  hx = lev.h(1); hy = lev.h(2);
  xf = dom(1) + (b(1)-1:b(2)) * hx;
  yc = dom(3) + ((b(3):b(4)) - 0.5) * hy;
  xc = dom(1) + ((b(1):b(2)) - 0.5) * hx;
  yf = dom(3) + (b(3)-1:b(4)) * hy;
  [XF, YCF] = ndgrid(xf, yc);
  lev.DeW{p} = evalgrid(prob.D, XF, YCF);
  [XCF, YF] = ndgrid(xc, yf);
  lev.DnS{p} = evalgrid(prob.D, XCF, YF);
  [XC, YC] = ndgrid(xc, yc);
  lev.CC{p}    = evalgrid(prob.C, XC, YC);
  lev.FTrue{p} = srcquad(prob.f, XC, YC, hx, hy, H.opts.SourceQuad);
  % F shares FTrue's buffer (copy-on-write); only patches whose F is
  % overwritten by the FAS coarse loading (i.e. parents) ever take a copy
  lev.F{p} = lev.FTrue{p};
  if any(lev.DeW{p}(:) <= 0) || any(lev.DnS{p}(:) <= 0)
    error('eb_make_level:D', 'D(x,y) must be strictly positive at faces.');
  end
  lev.par0(p) = mod(b(1) + b(3), 2);   % global checkerboard phase

  % A face is internal unless the patch touches the domain edge there
  % (order: west east south north).
  fint = [b(1) > 1, b(2) < n(1), b(3) > 1, b(4) < n(2)];
  lev.faceInternal{p} = fint;
  bcdat = struct('code', {0, 0, 0, 0}, 'g', {[], [], [], []});
  if ~fint(1), bcdat(1) = bcSide(prob, 1, dom(1), yc, true);  end
  if ~fint(2), bcdat(2) = bcSide(prob, 2, dom(2), yc, true);  end
  if ~fint(3), bcdat(3) = bcSide(prob, 3, xc, dom(3), false); end
  if ~fint(4), bcdat(4) = bcSide(prob, 4, xc, dom(4), false); end
  lev.bc{p} = bcdat;

  ID = H.K.invdiag(lev.DeW{p}, lev.DnS{p}, lev.CC{p}, hx, hy);
  lev.InvDiagS{p} = 1 ./ (1 ./ ID + ebBcFold(lev, p));

  lev.covered{p} = false(mx, my);
  lev.children{p} = [];
  lev.parents{p} = [];
  lev.cornerClass{p} = ones(1, 4, 'int8');
  lev.land{p} = struct('W', ones(my,1,'int8'), 'E', ones(my,1,'int8'), ...
                       'S', ones(mx,1,'int8'), 'N', ones(mx,1,'int8'));
  lev.sameOps{p} = [];
  lev.ops{p} = [];
  if k > 1
    cmx = mx/2; cmy = my/2;
    lev.QC{p} = zeros(cmx+2, cmy+2);
  end
end
end

% ------------------------------------------------------------------------
function s = bcSide(prob, side, xv, yv, isvert)
bcs = prob.bc(side);
if isvert
  g = bcs.g(xv * ones(size(yv)), yv);
  g = reshape(g, 1, []);
else
  g = bcs.g(xv, yv * ones(size(xv)));
  g = reshape(g, [], 1);
end
if isscalar(g)
  if isvert, g = g * ones(1, numel(yv)); else, g = g * ones(numel(xv), 1); end
end
if strcmp(bcs.type, 'dirichlet'), code = 1; else, code = 2; end
s = struct('code', code, 'g', g);
end

% ------------------------------------------------------------------------
function v = evalgrid(fh, X, Y)
v = fh(X, Y);
if isscalar(v), v = v * ones(size(X)); end
v = double(v);
if ~isequal(size(v), size(X))
  error('eb_make_level:eval', 'coefficient function returned wrong size.');
end
if any(~isfinite(v(:)))
  error('eb_make_level:eval', 'coefficient function returned non-finite values.');
end
end

% ------------------------------------------------------------------------
function F = srcquad(fh, XC, YC, hx, hy, nq)
%SRCQUAD  Cell value of the source: the midpoint value (nq <= 1) or the
%   nq-point tensor Gauss-Legendre CELL AVERAGE (nq > 1).  The cell average
%   makes the composite source quadrature  sum_visible h^2 F  equal the true
%   integral  \int_Omega f  to order 2*nq, because the visible cells tile
%   Omega exactly, so their averages sum to the exact integral.  That removes
%   the all-Neumann level drift at its origin (the drift is precisely the
%   composite-quadrature error of the source) while preserving exact discrete
%   conservation -- the scheme conserves whatever FTrue is, and now FTrue's
%   composite sum IS the true integral.  Differs from the midpoint value by
%   O(h^2), so the scheme stays second order.
if nq <= 1
  F = evalgrid(fh, XC, YC);
  return;
end
[xi, w] = gaussLegendre(nq);
w = w / 2;                              % normalize Gauss weights for an average
F = zeros(size(XC));
for i = 1:nq
  for j = 1:nq
    F = F + (w(i)*w(j)) * evalgrid(fh, XC + xi(i)*hx/2, YC + xi(j)*hy/2);
  end
end
end

% ------------------------------------------------------------------------
function [x, w] = gaussLegendre(n)
%gaussLegendre  Nodes/weights on [-1,1] by Golub-Welsch; sum(w) = 2.
if n == 1, x = 0; w = 2; return; end
beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
T = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(T);
[x, idx] = sort(diag(D));
w = 2 * (V(1, idx).^2).';
x = x(:);
end
