function [PDeW, PDnS, PCC] = ebParentCoeffs(H, k, p)
%ebParentCoeffs Transiently materialize the parent-level coefficients
%   over the footprint of patch p of level k (used by the FAS coarse
%   loading; not stored -- the values live in the parent patches).

% Footprint of the patch in parent-level cells (half the fine extent).
b = H.lev{k}.box(p, :);
cmx = (b(2) - b(1) + 1) / 2;
cmy = (b(4) - b(3) + 1) / 2;
ops = H.lev{k}.ops{p};
levC = H.lev{k-1};

% Face-coefficient arrays sized for the coarse footprint: E/W faces are
% (cmx+1) x cmy, N/S faces cmx x (cmy+1), cell-centered term cmx x cmy.
% NaN fill makes any gap in the precomputed slice tables loudly visible.
PDeW = nan(cmx+1, cmy);
PDnS = nan(cmx, cmy+1);
PCC  = nan(cmx, cmy);
for t = 1:numel(ops.pdew)
  s = ops.pdew(t).sr; d = ops.pdew(t).dr;
  PDeW(d(1):d(2), d(3):d(4)) = levC.DeW{ops.pdew(t).sp}(s(1):s(2), s(3):s(4));
end
for t = 1:numel(ops.pdns)
  s = ops.pdns(t).sr; d = ops.pdns(t).dr;
  PDnS(d(1):d(2), d(3):d(4)) = levC.DnS{ops.pdns(t).sp}(s(1):s(2), s(3):s(4));
end
for t = 1:numel(ops.pcc)
  s = ops.pcc(t).sr; d = ops.pcc(t).dr;
  PCC(d(1):d(2), d(3):d(4)) = levC.CC{ops.pcc(t).sp}(s(1):s(2), s(3):s(4));
end
end
