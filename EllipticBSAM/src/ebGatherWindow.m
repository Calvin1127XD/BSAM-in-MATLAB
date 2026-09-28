function W = ebGatherWindow(H, k, p)
%ebGatherWindow Materialize the (cmx+2)x(cmy+2) parent window of patch p
%   of level k from the parent level: interior via O(1) rect slices,
%   ring via the per-cell gather list (parent ghosts where the ring exits
%   the domain).

ops = H.lev{k}.ops{p};
b = H.lev{k}.box(p, :);
W = zeros((b(2)-b(1)+1)/2 + 2, (b(4)-b(3)+1)/2 + 2);
% Ring first (per-cell gathers), then the interior rect blocks on top.
for t = 1:numel(ops.ring)
  W(ops.ring(t).dst) = H.lev{k-1}.Q{ops.ring(t).sp}(ops.ring(t).src);
end
for t = 1:numel(ops.winrect)
  s = ops.winrect(t).sr; d = ops.winrect(t).dr;
  W(d(1):d(2), d(3):d(4)) = ...
    H.lev{k-1}.Q{ops.winrect(t).sp}(s(1):s(2), s(3):s(4));
end
end
