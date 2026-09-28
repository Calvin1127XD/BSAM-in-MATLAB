function vals = ebSample(H, xq, yq)
%ebSample Bilinear point sampling of the composite solution.
%   vals = ebSample(H, xq, yq) evaluates at the query points using the
%   finest patch whose box contains each point; interpolation is bilinear
%   on the ghost-padded patch array (ghosts carry coarse-fine/physical
%   closures, so sampling is well-defined up to the patch edges).

dom = H.prob.domain;
vals = nan(size(xq));
assigned = false(size(xq));

% Finest level first: once a point is assigned, coarser levels skip it.
for k = H.nlev:-1:H.rootIdx
  lev = H.lev{k};
  hx = lev.h(1); hy = lev.h(2);
  for p = 1:lev.np
    b = lev.box(p, :);
    x0 = dom(1) + (b(1)-1)*hx;  x1 = dom(1) + b(2)*hx;
    y0 = dom(3) + (b(3)-1)*hy;  y1 = dom(3) + b(4)*hy;
    in = ~assigned & xq >= x0 & xq <= x1 & yq >= y0 & yq <= y1;
    if ~any(in(:)), continue; end
    Q = lev.Q{p};
    [m2, n2] = size(Q);
    s = (xq(in) - x0)/hx + 1.5;          % fractional padded index
    t = (yq(in) - y0)/hy + 1.5;
    i0 = min(max(floor(s), 1), m2-1);
    j0 = min(max(floor(t), 1), n2-1);
    fr = s - i0;
    ft = t - j0;
    v = (1-fr).*(1-ft).*Q(sub2ind([m2 n2], i0,   j0))   ...
      +    fr .*(1-ft).*Q(sub2ind([m2 n2], i0+1, j0))   ...
      + (1-fr).*  ft .*Q(sub2ind([m2 n2], i0,   j0+1)) ...
      +    fr .*  ft .*Q(sub2ind([m2 n2], i0+1, j0+1));
    vals(in) = v;
    assigned(in) = true;
  end
end

if any(~assigned(:))
  error('eb_sample:cover', 'query points outside the domain.');
end
end
