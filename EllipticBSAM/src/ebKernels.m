function K = ebKernels()
%ebKernels Returns a struct of function handles for the small numerical
%   kernels shared across the solver (operator, smoother, transfer ops).
%   Kept in one file so the discretization is visible in one place.
K.op        = @op2d;
K.invdiag   = @invdiag2d;
K.relax     = @relax2d;
K.restrict2 = @restrict2;
K.prolong   = @prolongBilin;
end

% ------------------------------------------------------------------------
function L = op2d(Q, DeW, DnS, CC, hx, hy)
%OP2D  L(u) = -div(D grad u) + C u, conservative flux differencing.
%   Q padded (mx+2)x(my+2); DeW (mx+1)x(my) on x-faces; DnS (mx)x(my+1)
%   on y-faces; CC (mx)x(my).  Returns (mx)x(my).
mx = size(CC, 1); my = size(CC, 2);
% Fluxes on every face (differences of adjacent cells scaled by D), then
% the divergence as the difference of fluxes across each cell.
feW = DeW .* (Q(2:mx+2, 2:my+1) - Q(1:mx+1, 2:my+1));
fnS = DnS .* (Q(2:mx+1, 2:my+2) - Q(2:mx+1, 1:my+1));
L = -(feW(2:mx+1, :) - feW(1:mx, :)) / (hx*hx) ...
    -(fnS(:, 2:my+1) - fnS(:, 1:my)) / (hy*hy) ...
    + CC .* Q(2:mx+1, 2:my+1);
end

% ------------------------------------------------------------------------
function ID = invdiag2d(DeW, DnS, CC, hx, hy)
%INVDIAG2D Inverse of the operator diagonal (for point relaxation).
mx = size(CC, 1); my = size(CC, 2);
den = (DeW(2:mx+1, :) + DeW(1:mx, :)) / (hx*hx) ...
    + (DnS(:, 2:my+1) + DnS(:, 1:my)) / (hy*hy) + CC;
ID = 1 ./ den;
end

% ------------------------------------------------------------------------
function Q = relax2d(Q, F, DeW, DnS, InvDiag, bcd, par0, color, hx, hy, omega)
%RELAX2D One colored Gauss-Seidel pass with STRIDED checkerboard updates
%   (no mask arrays).  par0 = mod(ilo+jlo, 2) carries the patch's global
%   parity so simultaneous patch updates of one color form a true
%   half-sweep of the level union; color is 1 (red: global parity 0) or
%   2 (black).  Cells are updated by exactly solving their own equation
%   with neighbors frozen.
mx = size(F, 1); my = size(F, 2);
% Partial physical ghosts (quadratic Dirichlet: (8/3)g + (1/3)u_2; Neumann:
% h*g).  The u_1 ghost coupling is folded into InvDiag (ebBcFold), so the
% boundary cells perform exact ghost-eliminated Gauss-Seidel updates.  Set
% in-place here -- rather than on a separate copy passed in -- so the
% smoother owns the array it mutates: one allocation per call instead of
% two.  code 0 marks an internal (non-physical) face and is skipped.
if bcd(1).code == 1
  Q(1, 2:my+1) = (8/3)*bcd(1).g + (1/3)*Q(3, 2:my+1);
elseif bcd(1).code == 2
  Q(1, 2:my+1) = hx * bcd(1).g;
end
if bcd(2).code == 1
  Q(mx+2, 2:my+1) = (8/3)*bcd(2).g + (1/3)*Q(mx, 2:my+1);
elseif bcd(2).code == 2
  Q(mx+2, 2:my+1) = hx * bcd(2).g;
end
if bcd(3).code == 1
  Q(2:mx+1, 1) = (8/3)*bcd(3).g + (1/3)*Q(2:mx+1, 3);
elseif bcd(3).code == 2
  Q(2:mx+1, 1) = hy * bcd(3).g;
end
if bcd(4).code == 1
  Q(2:mx+1, my+2) = (8/3)*bcd(4).g + (1/3)*Q(2:mx+1, my);
elseif bcd(4).code == 2
  Q(2:mx+1, my+2) = hy * bcd(4).g;
end
hx2 = hx*hx; hy2 = hy*hy;
% Jacobi-style candidate values for ALL cells (rhs plus neighbor fluxes,
% divided by the diagonal); the strided writes below take only one color,
% which is what makes the pass Gauss-Seidel rather than Jacobi.
tmp = (F ...
    + Q(3:mx+2, 2:my+1) .* DeW(2:mx+1, :) / hx2 ...
    + Q(1:mx,   2:my+1) .* DeW(1:mx,   :) / hx2 ...
    + Q(2:mx+1, 3:my+2) .* DnS(:, 2:my+1) / hy2 ...
    + Q(2:mx+1, 1:my  ) .* DnS(:, 1:my  ) / hy2) .* InvDiag;
% local rows li in column lj are updated iff mod(li+lj+par0,2) == color-1
m = color - 1;
s1 = 2 - mod(m + par0 + 1, 2);    % start row for odd local columns
s2 = 2 - mod(m + par0, 2);        % start row for even local columns
inter = Q(2:mx+1, 2:my+1);
if omega == 1
  inter(s1:2:mx, 1:2:my) = tmp(s1:2:mx, 1:2:my);
  inter(s2:2:mx, 2:2:my) = tmp(s2:2:mx, 2:2:my);
else
  inter(s1:2:mx, 1:2:my) = omega*tmp(s1:2:mx, 1:2:my) ...
                         + (1-omega)*inter(s1:2:mx, 1:2:my);
  inter(s2:2:mx, 2:2:my) = omega*tmp(s2:2:mx, 2:2:my) ...
                         + (1-omega)*inter(s2:2:mx, 2:2:my);
end
Q(2:mx+1, 2:my+1) = inter;
end

% ------------------------------------------------------------------------
function R = restrict2(A)
%RESTRICT2 Conservative 4-point average (2m x 2n) -> (m x n).
R = 0.25 * (A(1:2:end-1, 1:2:end-1) + A(2:2:end, 1:2:end-1) ...
          + A(1:2:end-1, 2:2:end)   + A(2:2:end, 2:2:end));
end

% ------------------------------------------------------------------------
function Fout = prolongBilin(W)
%prolongBilin Cell-centered bilinear prolongation of a padded coarse
%   array W ((cm+2)x(cn+2), one ghost ring) to the fine grid (2cm x 2cn).
%   Weights 9/16, 3/16, 3/16, 1/16.
cm = size(W, 1) - 2; cn = size(W, 2) - 2;
% Shifted views of the coarse array: center, +-x neighbors, +-y neighbors,
% and the four diagonals; each fine cell blends its three nearest coarse
% neighbors with the center.
C  = W(2:cm+1, 2:cn+1);
Wm = W(1:cm,   2:cn+1);  Wp = W(3:cm+2, 2:cn+1);
Sm = W(2:cm+1, 1:cn  );  Sp = W(2:cm+1, 3:cn+2);
Dmm = W(1:cm,   1:cn  ); Dpm = W(3:cm+2, 1:cn  );
Dmp = W(1:cm,   3:cn+2); Dpp = W(3:cm+2, 3:cn+2);
Fout = zeros(2*cm, 2*cn);
Fout(1:2:end, 1:2:end) = (9*C + 3*Wm + 3*Sm + Dmm) / 16;
Fout(2:2:end, 1:2:end) = (9*C + 3*Wp + 3*Sm + Dpm) / 16;
Fout(1:2:end, 2:2:end) = (9*C + 3*Wm + 3*Sp + Dmp) / 16;
Fout(2:2:end, 2:2:end) = (9*C + 3*Wp + 3*Sp + Dpp) / 16;
end
