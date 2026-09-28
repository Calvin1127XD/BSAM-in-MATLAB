function err = ebtErrorNorms(H, ufun)
%ebtErrorNorms Composite error norms versus the exact solution, measured
%   on visible cells with the composite inner product (paper eq. 3.45):
%   ||e||_2 = sqrt( sum h_k^2 e^2 ),  ||e||_inf = max |e|.
%   Only the SINGULAR all-Neumann, C==0 problem is compared modulo an
%   additive constant.  Positive reaction fixes the solution's mean, so
%   its absolute error must be retained.  Source quadrature does not
%   change the operator's nullspace.  err.meanShift records the constant
%   removed from the error (zero for nonsingular problems).

if nargin < 2 || isempty(ufun)
  ufun = H.prob.exact;
end
if isempty(ufun)
  error('ebt_error_norms:noexact', 'no exact solution available.');
end
dom = H.prob.domain;

% Per-patch pointwise error at the cell centers (ghosts stripped).
E = cell(1, H.nlev);
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  for p = 1:lev.np
    b = lev.box(p, :);
    xc = dom(1) + ((b(1):b(2)) - 0.5) * lev.h(1);
    yc = dom(3) + ((b(3):b(4)) - 0.5) * lev.h(2);
    [XC, YC] = ndgrid(xc, yc);
    E{k}{p} = lev.Q{p}(2:end-1, 2:end-1) - ufun(XC, YC);
  end
end

% Gauge alignment is determined by the operator, independently of SourceQuad.
m = 0;
if H.coarsest.singular
  m = ebtCompositeMean(H, E);
  for k = H.rootIdx:H.nlev
    for p = 1:H.lev{k}.np
      E{k}{p} = E{k}{p} - m;
    end
  end
end

% Accumulate composite norms over visible cells, keeping a per-level
% breakdown for convergence tables.
linf = 0; l2sq = 0; dof = 0;
byLevel = [];
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  da = lev.h(1) * lev.h(2);
  li = 0; l2s = 0; nv = 0;
  for p = 1:lev.np
    vis = ~lev.covered{p};
    e = E{k}{p}(vis);
    if isempty(e), continue; end
    li = max(li, max(abs(e)));
    l2s = l2s + sum(e.^2) * da;
    nv = nv + numel(e);
  end
  linf = max(linf, li);
  l2sq = l2sq + l2s;
  dof = dof + nv;
  byLevel(end+1, :) = [k - H.rootIdx + 1, li, sqrt(l2s), nv]; %#ok<AGROW>
end

err.linf = linf;
err.l2 = sqrt(l2sq);
err.dof = dof;
err.byLevel = byLevel;
err.meanShift = m;
end
