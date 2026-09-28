function [U, X, Y] = ddToUniform(H, N0)
%ddToUniform Conservative projection of the composite solution onto a
%   uniform N0 x N0 grid of cell averages. Visible cells tile the domain
%   and every level is 2^k-aligned, so cells coarser than the target
%   replicate into whole blocks and finer cells average exactly; the
%   result is the L2-orthogonal projection onto piecewise constants on
%   the N0 grid. Used for Cauchy-difference (self-convergence) studies
%   between runs of different resolution.

dom = H.prob.domain;
U = zeros(N0, N0);
W = zeros(N0, N0);
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  n = lev.n(1);
  if lev.n(2) ~= n
    error('dd_to_uniform:square', 'square grids only.');
  end
  for p = 1:lev.np
    b = lev.box(p, :);
    Qi = lev.Q{p}(2:end-1, 2:end-1);
    vis = ~lev.covered{p};
    if ~any(vis(:)), continue; end
    if n >= N0
      % Level finer than the target: each target cell receives the mean
      % of its f x f fine sub-cells (accumulate values and weights).
      if mod(n, N0) ~= 0
        error('dd_to_uniform:align', 'N0 incompatible with level size %d.', n);
      end
      f = n / N0;
      [gi, gj] = ndgrid(b(1):b(2), b(3):b(4));
      sub1 = ceil(gi(vis) / f);
      sub2 = ceil(gj(vis) / f);
      lin = sub1 + (sub2 - 1) * N0;
      U(:) = U(:) + accumarray(lin, Qi(vis) / f^2, [N0*N0, 1]);
      W(:) = W(:) + accumarray(lin, ones(nnz(vis), 1) / f^2, [N0*N0, 1]);
    else
      % Level coarser than the target: replicate each visible cell into
      % its f x f block of target cells (piecewise-constant injection).
      if mod(N0, n) ~= 0
        error('dd_to_uniform:align', 'N0 incompatible with level size %d.', n);
      end
      f = N0 / n;
      R = repelem(Qi, f, f);
      M = repelem(vis, f, f);
      rows = (b(1)-1)*f + 1 : b(2)*f;
      cols = (b(3)-1)*f + 1 : b(4)*f;
      Ublk = U(rows, cols);
      Wblk = W(rows, cols);
      Ublk(M) = Ublk(M) + R(M);
      Wblk(M) = Wblk(M) + 1;
      U(rows, cols) = Ublk;
      W(rows, cols) = Wblk;
    end
  end
end
% Every target cell must have been written with total weight exactly 1 --
% anything else means the visible cells failed to tile the domain.
if max(abs(W(:) - 1)) > 1e-10
  error('dd_to_uniform:tile', 'visible cells do not tile the domain.');
end
U = U ./ W;
hx0 = (dom(2) - dom(1)) / N0;
hy0 = (dom(4) - dom(3)) / N0;
[X, Y] = ndgrid(dom(1) + ((1:N0) - 0.5)*hx0, dom(3) + ((1:N0) - 0.5)*hy0);
end
