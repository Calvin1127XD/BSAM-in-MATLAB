function err = ddErrorInside(H, uexact, aux)
%ddErrorInside Composite error norms restricted to the physical domain
%   Omega_1 = {psi < 0}, measured on visible cells (the diffuse-domain
%   solution approximates the sharp solution there; outside Omega_1 it is
%   an extension with no accuracy claim).

dom = H.prob.domain;
linf = 0; l2sq = 0; ncell = 0;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  da = lev.h(1) * lev.h(2);
  for p = 1:lev.np
    b = lev.box(p, :);
    xc = dom(1) + ((b(1):b(2)) - 0.5) * lev.h(1);
    yc = dom(3) + ((b(3):b(4)) - 0.5) * lev.h(2);
    [XC, YC] = ndgrid(xc, yc);
    inside = aux.psi(XC, YC) < 0;
    vis = ~lev.covered{p} & inside;
    if ~any(vis(:)), continue; end
    E = lev.Q{p}(2:end-1, 2:end-1) - uexact(XC, YC);
    e = E(vis);
    linf = max(linf, max(abs(e)));
    l2sq = l2sq + sum(e.^2) * da;
    ncell = ncell + numel(e);
  end
end
err.linf = linf;
err.l2 = sqrt(l2sq);
err.ncell = ncell;
end
