function H = ebCoeffFilldown(H)
%ebCoeffFilldown Restrict coefficients down the hierarchy so that every
%   coarse cell/face covered by a finer level carries the averaged fine
%   coefficient (paper eq. 4.15 compatibility): face D = mean of the two
%   fine faces, cell C = mean of the four fine cells.  Recomputes the
%   relaxation diagonals afterwards.

for k = H.nlev:-1:2
  levF = H.lev{k};
  levC = H.lev{k-1};
  for p = 1:levF.np
    b = levF.box(p, :);
    Il = (b(1)+1)/2; Ih = b(2)/2; Jl = (b(3)+1)/2; Jh = b(4)/2;
    mx = b(2) - b(1) + 1; my = b(4) - b(3) + 1;

    % Coarse x-face = mean of the two fine faces it contains (every other
    % fine face column, adjacent rows averaged); y-faces symmetrically;
    % cell term = 4-cell average.
    cDeW = 0.5 * (levF.DeW{p}(1:2:mx+1, 1:2:my-1) + levF.DeW{p}(1:2:mx+1, 2:2:my));
    cDnS = 0.5 * (levF.DnS{p}(1:2:mx-1, 1:2:my+1) + levF.DnS{p}(2:2:mx, 1:2:my+1));
    cCC  = H.K.restrict2(levF.CC{p});

    for q = 1:levC.np
      qb = levC.box(q, :);
      % cell-centered C over footprint
      r = rectIsect([Il Ih Jl Jh], qb);
      if ~isempty(r)
        H.lev{k-1}.CC{q}(r(1)-qb(1)+1 : r(2)-qb(1)+1, r(3)-qb(3)+1 : r(4)-qb(3)+1) ...
          = cCC(r(1)-Il+1 : r(2)-Il+1, r(3)-Jl+1 : r(4)-Jl+1);
      end
      % x-faces: global face index range [Il, Ih+1] x cells [Jl, Jh]
      rf = rectIsect([Il, Ih+1, Jl, Jh], [qb(1), qb(2)+1, qb(3), qb(4)]);
      if ~isempty(rf)
        H.lev{k-1}.DeW{q}(rf(1)-qb(1)+1 : rf(2)-qb(1)+1, rf(3)-qb(3)+1 : rf(4)-qb(3)+1) ...
          = cDeW(rf(1)-Il+1 : rf(2)-Il+1, rf(3)-Jl+1 : rf(4)-Jl+1);
      end
      % y-faces: cells [Il, Ih] x global face range [Jl, Jh+1]
      rf = rectIsect([Il, Ih, Jl, Jh+1], [qb(1), qb(2), qb(3), qb(4)+1]);
      if ~isempty(rf)
        H.lev{k-1}.DnS{q}(rf(1)-qb(1)+1 : rf(2)-qb(1)+1, rf(3)-qb(3)+1 : rf(4)-qb(3)+1) ...
          = cDnS(rf(1)-Il+1 : rf(2)-Il+1, rf(3)-Jl+1 : rf(4)-Jl+1);
      end
    end
  end
end

% refresh relaxation diagonals everywhere (with the boundary fold)
% InvDiagS stores 1/diag for Gauss-Seidel; ebBcFold adds the extra diagonal
% terms from eliminating physical-boundary ghosts.
for k = 1:H.nlev
  lev = H.lev{k};
  for p = 1:lev.np
    ID = H.K.invdiag(lev.DeW{p}, lev.DnS{p}, lev.CC{p}, lev.h(1), lev.h(2));
    H.lev{k}.InvDiagS{p} = 1 ./ (1 ./ ID + ebBcFold(H.lev{k}, p));
  end
end
end

function r = rectIsect(a, b)
r = [max(a(1), b(1)), min(a(2), b(2)), max(a(3), b(3)), min(a(4), b(4))];
if r(1) > r(2) || r(3) > r(4), r = []; end
end
