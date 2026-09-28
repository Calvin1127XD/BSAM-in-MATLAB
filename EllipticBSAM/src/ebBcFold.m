function add = ebBcFold(lev, p)
%ebBcFold Diagonal increment that folds the physical-boundary ghost
%   elimination into point relaxation, making boundary-cell updates exact
%   Gauss-Seidel for the ghost-eliminated equations:
%     quadratic Dirichlet ghost  u_g = (8/3)g - 2 u_1 + (1/3) u_2
%        -> diagonal += 2 D_face / h^2   (the -2 u_1 term),
%     Neumann ghost              u_g = u_1 + h g
%        -> diagonal -= D_face / h^2     (no u_1 coupling remains).
%   The remaining ghost parts ((8/3)g + (1/3)u_2, resp. h g) are written
%   into the ghost cells by ebRelaxLevel before each colored pass.

bcd = lev.bc{p};
mx = size(lev.CC{p}, 1); my = size(lev.CC{p}, 2);
hx2 = lev.h(1)^2; hy2 = lev.h(2)^2;
add = zeros(mx, my);

% Sides in order west/east/south/north; code 1 = Dirichlet, 2 = Neumann.
% West side touches interior column 1 through face coefficient DeW(1,:).
if bcd(1).code == 1
  add(1, :) = add(1, :) + 2 * lev.DeW{p}(1, :) / hx2;
elseif bcd(1).code == 2
  add(1, :) = add(1, :) - lev.DeW{p}(1, :) / hx2;
end
% East side: interior column mx, face DeW(mx+1,:).
if bcd(2).code == 1
  add(mx, :) = add(mx, :) + 2 * lev.DeW{p}(mx+1, :) / hx2;
elseif bcd(2).code == 2
  add(mx, :) = add(mx, :) - lev.DeW{p}(mx+1, :) / hx2;
end
% South side: interior row 1, face DnS(:,1).
if bcd(3).code == 1
  add(:, 1) = add(:, 1) + 2 * lev.DnS{p}(:, 1) / hy2;
elseif bcd(3).code == 2
  add(:, 1) = add(:, 1) - lev.DnS{p}(:, 1) / hy2;
end
% North side: interior row my, face DnS(:,my+1).
if bcd(4).code == 1
  add(:, my) = add(:, my) + 2 * lev.DnS{p}(:, my+1) / hy2;
elseif bcd(4).code == 2
  add(:, my) = add(:, my) - lev.DnS{p}(:, my+1) / hy2;
end
end
