function H = ebtFillGhosts(H, k)
%ebtFillGhosts Fill all ghost cells of level k, in the canonical order
%   (the Fortran SetGhost priority):
%   1. coarse-fine interpolation from the SINGLE parent (QIF/LIF),
%   2. same-level exchange (overwrites shared segments),
%   3. physical boundary conditions (quadratic Dirichlet / Neumann).
%
%   TREE VERSION of step 1: the coarse ring QC is refreshed by four
%   sliced copies from the one parent's padded array (ops.ring).  Where
%   the ring leaves the parent's interior, the parent's own ghost cells
%   are read -- they hold sibling/physical/coarser closures synchronized
%   at the parent's level, exactly like Fortran InterpolateCoarseFine
%   reading parent%q only.

method = H.opts.GhostInterp;
% Hoist the level's cell arrays into locals; mutate the locals and write the
% whole cell array back once at the end, instead of indexing H.lev{k}.Q{p}
% (cell-index -> field -> cell-index) on every patch of every ghost pass.
np = H.lev{k}.np;
Q  = H.lev{k}.Q;
faceInternal = H.lev{k}.faceInternal;  cornerClass = H.lev{k}.cornerClass;
boxk = H.lev{k}.box;  opsk = H.lev{k}.ops;  sameAll = H.lev{k}.sameOps;
bcAll = H.lev{k}.bc;  hx = H.lev{k}.h(1);  hy = H.lev{k}.h(2);

% ---- 1. coarse-fine interpolation (single-parent ring)
if k > 1
  QCk = H.lev{k}.QC;
  Qpar = H.lev{k-1}.Q;          % parent solutions, read only here
  for p = 1:np
    fint = faceInternal{p};
    cCF = (cornerClass{p} == 3);
    if ~any(fint) && ~any(cCF), continue; end
    QC = QCk{p};
    if isempty(QC)   % QC is released at solve end; reallocate on demand
      b = boxk(p, :);
      QC = zeros((b(2)-b(1)+1)/2 + 2, (b(4)-b(3)+1)/2 + 2);
    end
    ops = opsk{p};
    Qq = Qpar{ops.sp};
    for t = 1:4
      s = ops.ring(t).sr; d = ops.ring(t).dr;
      QC(d(1):d(2), d(3):d(4)) = Qq(s(1):s(2), s(3):s(4));
    end
    QCk{p} = QC;
    Q{p} = ebtInterpCf(Q{p}, QC, fint, cCF, method);
  end
  H.lev{k}.QC = QCk;
end

% ---- 2. same-level exchange (reads sibling interiors only)
for p = 1:np
  same = sameAll{p};
  if isempty(same), continue; end
  Qp = Q{p};
  for t = 1:numel(same)
    Qp(same(t).dst) = Q{same(t).sp}(same(t).src);
  end
  Q{p} = Qp;
end

% ---- 3. physical boundaries
for p = 1:np
  bcd = bcAll{p};
  if all([bcd.code] == 0) && all(cornerClass{p} ~= 1), continue; end
  Qp = Q{p};
  mx = size(Qp, 1) - 2; my = size(Qp, 2) - 2;

  if bcd(1).code == 1      % west Dirichlet (quadratic ghost)
    Qp(1, 2:my+1) = (8/3)*bcd(1).g - 2*Qp(2, 2:my+1) + (1/3)*Qp(3, 2:my+1);
  elseif bcd(1).code == 2  % west Neumann: du/dn = g, n = -x
    Qp(1, 2:my+1) = Qp(2, 2:my+1) + hx*bcd(1).g;
  end
  % east side
  if bcd(2).code == 1
    Qp(mx+2, 2:my+1) = (8/3)*bcd(2).g - 2*Qp(mx+1, 2:my+1) + (1/3)*Qp(mx, 2:my+1);
  elseif bcd(2).code == 2
    Qp(mx+2, 2:my+1) = Qp(mx+1, 2:my+1) + hx*bcd(2).g;
  end
  % south side
  if bcd(3).code == 1
    Qp(2:mx+1, 1) = (8/3)*bcd(3).g - 2*Qp(2:mx+1, 2) + (1/3)*Qp(2:mx+1, 3);
  elseif bcd(3).code == 2
    Qp(2:mx+1, 1) = Qp(2:mx+1, 2) + hy*bcd(3).g;
  end
  % north side
  if bcd(4).code == 1
    Qp(2:mx+1, my+2) = (8/3)*bcd(4).g - 2*Qp(2:mx+1, my+1) + (1/3)*Qp(2:mx+1, my);
  elseif bcd(4).code == 2
    Qp(2:mx+1, my+2) = Qp(2:mx+1, my+1) + hy*bcd(4).g;
  end

  % corner ghosts classified physical: bilinear-consistent extrapolation
  cc = cornerClass{p};
  if cc(1) == 1, Qp(1,1)       = Qp(2,1) + Qp(1,2) - Qp(2,2); end
  if cc(2) == 1, Qp(mx+2,1)    = Qp(mx+1,1) + Qp(mx+2,2) - Qp(mx+1,2); end
  if cc(3) == 1, Qp(1,my+2)    = Qp(2,my+2) + Qp(1,my+1) - Qp(2,my+1); end
  if cc(4) == 1, Qp(mx+2,my+2) = Qp(mx+1,my+2) + Qp(mx+2,my+1) - Qp(mx+1,my+1); end

  Q{p} = Qp;
end

H.lev{k}.Q = Q;            % single write-back of the level's solution array
end
