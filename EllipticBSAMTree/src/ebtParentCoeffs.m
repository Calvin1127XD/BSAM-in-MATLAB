function [PDeW, PDnS, PCC] = ebtParentCoeffs(H, k, p)
%ebtParentCoeffs Transiently materialize the parent-level coefficients
%   over the footprint of patch p of level k (used by the FAS coarse
%   loading; not stored -- the values live in the parent patch).
%   TREE VERSION: single slices of the ONE parent's arrays (descriptors
%   ops.pcc/pdew/pdns from ebtBuildOps).

ops = H.lev{k}.ops{p};
levC = H.lev{k-1};
q = ops.sp;

s = ops.pdew.sr;
PDeW = levC.DeW{q}(s(1):s(2), s(3):s(4));
s = ops.pdns.sr;
PDnS = levC.DnS{q}(s(1):s(2), s(3):s(4));
s = ops.pcc.sr;
PCC = levC.CC{q}(s(1):s(2), s(3):s(4));
end
