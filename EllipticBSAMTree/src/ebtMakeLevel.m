function lev = ebtMakeLevel(H, k, n, boxes, parents)
%ebtMakeLevel Create a level struct with patches [ilo ihi jlo jhi] (global
%   cell indices on the level-k grid of n = [nx ny] cells) and an EXPLICIT
%   single-parent lineage: parents(p) is the index of THE one level-(k-1)
%   patch that contains patch p (0 on level 1).
%
%   This is the tree-based variant (BSAM 2.0 treeops design): every patch
%   stores two index systems, exactly like the Fortran code:
%     lev.box(p,:)     "mglobal"  - global cell range at the patch's own
%                                   level (drives all SAME-LEVEL work:
%                                   sibling exchange, owner maps, plots);
%     lev.mbounds(p,:) "mbounds"  - the patch footprint in ITS PARENT's
%                                   local cell space, 1-based
%                                   (drives ALL coarse-fine transfer).
%   Invariants enforced here (Fortran MakeNewGrid / NewSubGrids checks):
%     * global box is odd-start/even-end in both directions (k >= 2), so
%       the patch coarsens exactly onto parent cells;
%     * the coarsened footprint lies inside the single parent:
%       1 <= mbounds(1) <= mbounds(2) <= parent mx (and same in y) -- a
%       child may touch its parent's edge but never extend past it.
%
%   Evaluates the analytic coefficients/source, builds red-black masks
%   (global parity), classifies physical faces and evaluates boundary data.

prob = H.prob;
dom = prob.domain;
lev.n = n;
lev.h = [(dom(2)-dom(1))/n(1), (dom(4)-dom(3))/n(2)];
lev.np = size(boxes, 1);
lev.box = boxes;
lev.sourceShift = 0; % constant removed from the originally sampled source

if nargin < 5 || isempty(parents)
  parents = zeros(lev.np, 1);
end
lev.parent = parents(:);
lev.mbounds = nan(lev.np, 4);

for p = 1:lev.np
  b = boxes(p, :);
  if any(b(1:2) < 1 | b(1:2) > n(1)) || any(b(3:4) < 1 | b(3:4) > n(2)) ...
      || b(1) > b(2) || b(3) > b(4)
    error('ebt_make_level:box', 'invalid patch box.');
  end

  % ---- lineage record (k >= 2): footprint in the parent's local cells
  if k > 1
    if mod(b(1),2) ~= 1 || mod(b(2),2) ~= 0 || mod(b(3),2) ~= 1 || mod(b(4),2) ~= 0
      error('ebt_make_level:align', ...
        'patch %d of level %d not aligned to parent cells.', p, k);
    end
    q = lev.parent(p);
    if q < 1 || q > H.lev{k-1}.np
      error('ebt_make_level:parent', ...
        'patch %d of level %d has invalid parent index %d.', p, k, q);
    end
    qb = H.lev{k-1}.box(q, :);
    % coarsened footprint in level-(k-1) GLOBAL cells, then parent-local
    fp = [(b(1)+1)/2, b(2)/2, (b(3)+1)/2, b(4)/2];
    mb = [fp(1)-qb(1)+1, fp(2)-qb(1)+1, fp(3)-qb(3)+1, fp(4)-qb(3)+1];
    pmx = qb(2) - qb(1) + 1; pmy = qb(4) - qb(3) + 1;
    if mb(1) < 1 || mb(2) > pmx || mb(3) < 1 || mb(4) > pmy
      error('ebt_make_level:containment', ...
        ['patch %d of level %d (footprint [%d %d %d %d]) is not contained ' ...
         'in its parent patch %d ([%d %d %d %d]) - single-parent lineage violated.'], ...
        p, k, fp, q, qb);
    end
    lev.mbounds(p, :) = mb;
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
    error('ebt_make_level:D', 'D(x,y) must be strictly positive at faces.');
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
  lev.InvDiagS{p} = 1 ./ (1 ./ ID + ebtBcFold(lev, p));

  lev.covered{p} = false(mx, my);
  lev.children{p} = [];
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
  error('ebt_make_level:eval', 'coefficient function returned wrong size.');
end
if any(~isfinite(v(:)))
  error('ebt_make_level:eval', 'coefficient function returned non-finite values.');
end
end

% ------------------------------------------------------------------------
function F = srcquad(fh, XC, YC, hx, hy, nq)
%SRCQUAD  Cell value of the source: the midpoint value (nq <= 1) or the
%   nq-point tensor Gauss-Legendre CELL AVERAGE (nq > 1).  The cell average
%   makes the composite source quadrature  sum_visible h^2 F  equal the true
%   integral  \int_Omega f  to order 2*nq, because the visible cells tile
%   Omega exactly, so their averages sum to the exact integral.  Differs
%   from the midpoint value by O(h^2), so the scheme stays second order.
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
