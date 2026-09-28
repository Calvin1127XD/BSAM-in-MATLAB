function Q = ebtInterpCf(Q, QC, faceInternal, cornerCF, method)
%ebtInterpCf Coarse-fine ghost interpolation for one fine patch (QIF/LIF).
%
%   Q  = ebtInterpCf(Q, QC, faceInternal, cornerCF, method)
%
%   Q  : (mx+2)x(my+2) padded fine-patch solution.  Ghost cells on the
%        faces marked internal are overwritten.
%   QC : (cmx+2)x(cmy+2) parent-level window, cmx = mx/2, cmy = my/2.
%        Interior (2:cmx+1, 2:cmy+1) holds the parent values under the
%        patch footprint, the ring holds the parent values just outside
%        (parent physical ghosts where the ring exits the domain).
%   faceInternal : logical [W E S N]  - face is a coarse-fine/same-level
%        interface (not a physical boundary).  Formulas are applied to the
%        whole face; same-level segments are overwritten afterwards by the
%        same-level exchange, exactly as in Feng et al. (2018) Fig. 4.
%   cornerCF : logical [SW SE NW NE]  - corner ghost cell lies over the
%        coarse level (fill with the diagonal stencil).
%   method : 'qif' (quadratic, paper eq. 3.12) or 'lif' (linear, eq. 3.19).
%
%   QIF: the ghost value is the quadratic interpolant along the 45-degree
%   line through the diagonally adjacent coarse cell center and the two
%   interior fine cells on that line (the "Lambda" stencil of Fig. 4a):
%
%       u_ghost = (8/15) u_coarse + (2/3) u_f1 - (1/5) u_f2 .
%
%   At the two end cells of every face the diagonal line would leave the
%   patch, so the coarse value at the ghost's tangential station is first
%   built by tangential quadratic interpolation with weights
%   (5/32, 15/16, -3/32) and the normal-direction quadratic is applied
%   along the grid line (the "rule-exception" cells of Fig. 4a).
%
%   LIF: linear interpolation along the same diagonal,
%       u_ghost = (2/3) u_coarse + (2/3) u_f(same row) - (1/3) u_f(partner row),
%   which needs no end exceptions (paper eq. 4.13-4.14).

mx = size(Q, 1) - 2;
my = size(Q, 2) - 2;
cmx = mx / 2;
cmy = my / 2;

a1 = -1/5;  a2 = 2/3;  a3 = 8/15;      % normal-direction quadratic
b1 = 5/32;  b2 = 15/16; b3 = -3/32;    % tangential quadratic (end cells)
qif = strcmpi(method, 'qif');

% ---------------------------------------------------------------- WEST
if faceInternal(1)
  if qif
    rho = 3:2:my-1;   % fine rows r = 2:2:my-2 (upper cell of parent pair)
    Q(1, rho) = a1*Q(3, rho+2) + a2*Q(2, rho+1) + a3*QC(1, (rho+1)/2);
    rho = 4:2:my;     % fine rows r = 3:2:my-1 (lower cell of parent pair)
    Q(1, rho) = a1*Q(3, rho-2) + a2*Q(2, rho-1) + a3*QC(1, (rho+2)/2);
    Q(1, my+1) = a1*Q(3, my+1) + a2*Q(2, my+1) ...
               + a3*(b1*QC(1, cmy+2) + b2*QC(1, cmy+1) + b3*QC(1, cmy));
    Q(1, 2)    = a1*Q(3, 2) + a2*Q(2, 2) ...
               + a3*(b1*QC(1, 1) + b2*QC(1, 2) + b3*QC(1, 3));
  else
    rho = 3:2:my+1;   % r even: partner row below
    Q(1, rho) = (2/3)*QC(1, (rho+1)/2) + (2/3)*Q(2, rho) - (1/3)*Q(2, rho-1);
    rho = 2:2:my;     % r odd: partner row above
    Q(1, rho) = (2/3)*QC(1, (rho+2)/2) + (2/3)*Q(2, rho) - (1/3)*Q(2, rho+1);
  end
end

% ---------------------------------------------------------------- EAST
if faceInternal(2)
  if qif
    % Mirror of WEST with the diagonal read from the east ring column.
    rho = 3:2:my-1;
    Q(mx+2, rho) = a1*Q(mx, rho+2) + a2*Q(mx+1, rho+1) + a3*QC(cmx+2, (rho+1)/2);
    rho = 4:2:my;
    Q(mx+2, rho) = a1*Q(mx, rho-2) + a2*Q(mx+1, rho-1) + a3*QC(cmx+2, (rho+2)/2);
    Q(mx+2, my+1) = a1*Q(mx, my+1) + a2*Q(mx+1, my+1) ...
                  + a3*(b1*QC(cmx+2, cmy+2) + b2*QC(cmx+2, cmy+1) + b3*QC(cmx+2, cmy));
    Q(mx+2, 2)    = a1*Q(mx, 2) + a2*Q(mx+1, 2) ...
                  + a3*(b1*QC(cmx+2, 1) + b2*QC(cmx+2, 2) + b3*QC(cmx+2, 3));
  else
    rho = 3:2:my+1;
    Q(mx+2, rho) = (2/3)*QC(cmx+2, (rho+1)/2) + (2/3)*Q(mx+1, rho) - (1/3)*Q(mx+1, rho-1);
    rho = 2:2:my;
    Q(mx+2, rho) = (2/3)*QC(cmx+2, (rho+2)/2) + (2/3)*Q(mx+1, rho) - (1/3)*Q(mx+1, rho+1);
  end
end

% ---------------------------------------------------------------- SOUTH
if faceInternal(3)
  if qif
    % Same construction rotated 90 degrees: tangential index runs along x.
    kap = 3:2:mx-1;
    Q(kap, 1) = a1*Q(kap+2, 3) + a2*Q(kap+1, 2) + a3*QC((kap+1)/2, 1);
    kap = 4:2:mx;
    Q(kap, 1) = a1*Q(kap-2, 3) + a2*Q(kap-1, 2) + a3*QC((kap+2)/2, 1);
    Q(mx+1, 1) = a1*Q(mx+1, 3) + a2*Q(mx+1, 2) ...
               + a3*(b1*QC(cmx+2, 1) + b2*QC(cmx+1, 1) + b3*QC(cmx, 1));
    Q(2, 1)    = a1*Q(2, 3) + a2*Q(2, 2) ...
               + a3*(b1*QC(1, 1) + b2*QC(2, 1) + b3*QC(3, 1));
  else
    kap = 3:2:mx+1;
    Q(kap, 1) = (2/3)*QC((kap+1)/2, 1) + (2/3)*Q(kap, 2) - (1/3)*Q(kap-1, 2);
    kap = 2:2:mx;
    Q(kap, 1) = (2/3)*QC((kap+2)/2, 1) + (2/3)*Q(kap, 2) - (1/3)*Q(kap+1, 2);
  end
end

% ---------------------------------------------------------------- NORTH
if faceInternal(4)
  if qif
    kap = 3:2:mx-1;
    Q(kap, my+2) = a1*Q(kap+2, my) + a2*Q(kap+1, my+1) + a3*QC((kap+1)/2, cmy+2);
    kap = 4:2:mx;
    Q(kap, my+2) = a1*Q(kap-2, my) + a2*Q(kap-1, my+1) + a3*QC((kap+2)/2, cmy+2);
    Q(mx+1, my+2) = a1*Q(mx+1, my) + a2*Q(mx+1, my+1) ...
                  + a3*(b1*QC(cmx+2, cmy+2) + b2*QC(cmx+1, cmy+2) + b3*QC(cmx, cmy+2));
    Q(2, my+2)    = a1*Q(2, my) + a2*Q(2, my+1) ...
                  + a3*(b1*QC(1, cmy+2) + b2*QC(2, cmy+2) + b3*QC(3, cmy+2));
  else
    kap = 3:2:mx+1;
    Q(kap, my+2) = (2/3)*QC((kap+1)/2, cmy+2) + (2/3)*Q(kap, my+1) - (1/3)*Q(kap-1, my+1);
    kap = 2:2:mx;
    Q(kap, my+2) = (2/3)*QC((kap+2)/2, cmy+2) + (2/3)*Q(kap, my+1) - (1/3)*Q(kap+1, my+1);
  end
end

% ---------------------------------------------------------------- CORNERS
% Diagonal stencil through the diagonally adjacent coarse cell and the two
% interior fine cells on the patch diagonal.
if cornerCF(1)   % SW
  if qif
    Q(1, 1) = a1*Q(3, 3) + a2*Q(2, 2) + a3*QC(1, 1);
  else
    Q(1, 1) = (2/3)*QC(1, 1) + (1/3)*Q(2, 2);
  end
end
if cornerCF(2)   % SE
  if qif
    Q(mx+2, 1) = a1*Q(mx, 3) + a2*Q(mx+1, 2) + a3*QC(cmx+2, 1);
  else
    Q(mx+2, 1) = (2/3)*QC(cmx+2, 1) + (1/3)*Q(mx+1, 2);
  end
end
if cornerCF(3)   % NW
  if qif
    Q(1, my+2) = a1*Q(3, my) + a2*Q(2, my+1) + a3*QC(1, cmy+2);
  else
    Q(1, my+2) = (2/3)*QC(1, cmy+2) + (1/3)*Q(2, my+1);
  end
end
if cornerCF(4)   % NE
  if qif
    Q(mx+2, my+2) = a1*Q(mx, my) + a2*Q(mx+1, my+1) + a3*QC(cmx+2, cmy+2);
  else
    Q(mx+2, my+2) = (2/3)*QC(cmx+2, cmy+2) + (1/3)*Q(mx+1, my+1);
  end
end
end
