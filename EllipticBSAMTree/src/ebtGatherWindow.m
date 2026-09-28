function W = ebtGatherWindow(H, k, p)
%ebtGatherWindow Materialize the (cmx+2)x(cmy+2) parent window of patch p
%   of level k.  TREE VERSION: the window -- interior AND ring -- is one
%   contiguous slice of the single parent's padded solution array
%   (ops.win, computed in ebtBuildOps).  Ring cells that leave the
%   parent's interior read the parent's own ghost cells, which carry the
%   synchronized coarse-level values (sibling data, coarser-level
%   interpolant, or physical closure).

ops = H.lev{k}.ops{p};
s = ops.win;
W = H.lev{k-1}.Q{ops.sp}(s(1):s(2), s(3):s(4));
end
