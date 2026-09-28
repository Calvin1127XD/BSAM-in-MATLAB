function [m, area] = ebtCompositeMean(H, E)
%ebtCompositeMean Mean over visible cells of the composite solution, or
%   of the per-patch cell arrays E{k}{p} if given.
tot = 0; area = 0;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  % Cell area at this level; each level contributes only its visible
  % (not-covered-by-finer) cells, so the sum is a true composite integral.
  da = lev.h(1) * lev.h(2);
  for p = 1:lev.np
    vis = ~lev.covered{p};
    if nargin >= 2
      v = E{k}{p}(vis);
    else
      % Default field is the solution; strip the ghost frame first.
      Qi = lev.Q{p}(2:end-1, 2:end-1);
      v = Qi(vis);
    end
    tot = tot + sum(v) * da;
    area = area + nnz(vis) * da;
  end
end
m = tot / max(area, realmin);
end
