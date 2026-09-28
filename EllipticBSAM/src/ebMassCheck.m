function md = ebMassCheck(H)
%ebMassCheck Discrete conservation identity of the composite scheme.
%   Summing the flux-form equations over all visible cells, all interior
%   fluxes (including coarse-fine, thanks to the flux-balance corrections)
%   telescope, leaving (paper Theorem 3.3 / eq. 3.47):
%
%       sum_chi (C u - f) h^2  =  oint D du/dn ds   (discrete),
%
%   where the boundary integral uses the ghost-based one-sided fluxes.
%   Returns md.defect (lhs - rhs) and a normalized md.rel.  At the solved
%   fixed point of the conservative (QIF/LIF + corrections) scheme this is
%   zero to roundoff; without the corrections it is O(h) at the interface.

% Left side: composite integral of (C u - f) over the visible cells.
lhs = 0; scale = 0;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  da = lev.h(1) * lev.h(2);
  for p = 1:lev.np
    vis = ~lev.covered{p};
    Qi = lev.Q{p}(2:end-1, 2:end-1);
    v = lev.CC{p}(vis) .* Qi(vis) - lev.FTrue{p}(vis);
    lhs = lhs + sum(v) * da;
    scale = scale + sum(abs(lev.FTrue{p}(vis))) * da;
  end
end

% Right side: outward boundary flux through physical faces, computed from
% the ghost closure (one-sided ghost-minus-interior differences), skipping
% boundary cells covered by a finer level.
rhs = 0;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  hx = lev.h(1); hy = lev.h(2);
  for p = 1:lev.np
    Q = lev.Q{p};
    mx = size(Q, 1) - 2; my = size(Q, 2) - 2;
    cov = lev.covered{p};
    bcd = lev.bc{p};
    if bcd(1).code > 0
      m = ~cov(1, :);
      fl = lev.DeW{p}(1, :) .* (Q(1, 2:my+1) - Q(2, 2:my+1)) / hx;
      rhs = rhs + sum(fl(m)) * hy;
    end
    if bcd(2).code > 0
      m = ~cov(mx, :);
      fl = lev.DeW{p}(mx+1, :) .* (Q(mx+2, 2:my+1) - Q(mx+1, 2:my+1)) / hx;
      rhs = rhs + sum(fl(m)) * hy;
    end
    if bcd(3).code > 0
      m = ~cov(:, 1);
      fl = lev.DnS{p}(:, 1) .* (Q(2:mx+1, 1) - Q(2:mx+1, 2)) / hy;
      rhs = rhs + sum(fl(m)) * hx;
    end
    if bcd(4).code > 0
      m = ~cov(:, my);
      fl = lev.DnS{p}(:, my+1) .* (Q(2:mx+1, my+2) - Q(2:mx+1, my+1)) / hy;
      rhs = rhs + sum(fl(m)) * hx;
    end
  end
end

md.lhs = lhs;
md.flux = rhs;
md.defect = lhs - rhs;
md.rel = abs(md.defect) / max([scale, abs(lhs), abs(rhs), realmin]);
end
