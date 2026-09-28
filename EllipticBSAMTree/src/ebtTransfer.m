function H = ebtTransfer(H, k, what)
%ebtTransfer Inter-level data motion for the AFAS cycle.
%   ebtTransfer(H, k, 'restrict')  fine solution -> QC interior + parent Q
%   ebtTransfer(H, k, 'ring')      parent Q -> QC ring (sliced copies)
%   ebtTransfer(H, k, 'correct')   coarse-grid correction:
%                                   Q += P(window - QC), bilinear prolong
%   ebtTransfer(H, k, 'filldown')  cascade restrict + ghost fills k..2
%
%   TREE VERSION: every patch has exactly one parent, so 'restrict' is a
%   single sliced write into the parent (Fortran RestrictCCSolution:
%   parent%q(mbounds) = qc), 'ring' copies the four window-perimeter
%   slices from the parent's padded array, and 'correct' reads the whole
%   window as one parent slice.

switch what
  case 'restrict'
    restrict2 = H.K.restrict2;
    np = H.lev{k}.np;
    Qk = H.lev{k}.Q;  QCk = H.lev{k}.QC;  Qpar = H.lev{k-1}.Q;
    opsk = H.lev{k}.ops;
    for p = 1:np
      % 2x2 cell average of the patch interior; goes to the patch's own
      % coarse shadow QC and, as ONE sliced write, to the single parent.
      R = restrict2(Qk{p}(2:end-1, 2:end-1));
      qc = QCk{p};
      if isempty(qc)   % QC is released at solve end; reallocate on demand
        qc = zeros(size(R) + 2);
      end
      qc(2:end-1, 2:end-1) = R;  QCk{p} = qc;
      ops = opsk{p};
      d = ops.usc.dr;
      Qpar{ops.sp}(d(1):d(2), d(3):d(4)) = R;
    end
    H.lev{k}.QC = QCk;
    H.lev{k-1}.Q = Qpar;

  case 'ring'
    np = H.lev{k}.np;
    QCk = H.lev{k}.QC;  Qpar = H.lev{k-1}.Q;
    opsk = H.lev{k}.ops;  boxk = H.lev{k}.box;
    for p = 1:np
      QC = QCk{p};
      if isempty(QC)   % QC is released at solve end; reallocate on demand
        b = boxk(p, :);
        QC = zeros((b(2)-b(1)+1)/2 + 2, (b(4)-b(3)+1)/2 + 2);
      end
      % The four perimeter slices (W/E columns incl. corners, S/N rows)
      % all come from the one parent's padded array -- proper nesting
      % guarantees the parent's own ghosts cover any overhang.
      ops = opsk{p};
      Qq = Qpar{ops.sp};
      for t = 1:4
        s = ops.ring(t).sr; d = ops.ring(t).dr;
        QC(d(1):d(2), d(3):d(4)) = Qq(s(1):s(2), s(3):s(4));
      end
      QCk{p} = QC;
    end
    H.lev{k}.QC = QCk;

  case 'correct'
    prolong = H.K.prolong;
    np = H.lev{k}.np;
    Qk = H.lev{k}.Q;  QCk = H.lev{k}.QC;  Qpar = H.lev{k-1}.Q;
    opsk = H.lev{k}.ops;
    for p = 1:np
      % FAS correction: the whole (footprint + ring) window is one slice
      % of the parent; subtract the pre-descent shadow QC and prolong the
      % remaining coarse-grid correction onto the fine interior.
      ops = opsk{p};
      s = ops.win;
      W = Qpar{ops.sp}(s(1):s(2), s(3):s(4));
      E = W - QCk{p};
      Qp = Qk{p};
      Qp(2:end-1, 2:end-1) = Qp(2:end-1, 2:end-1) + prolong(E);
      Qk{p} = Qp;
    end
    H.lev{k}.Q = Qk;

  case 'filldown'
    % Cascade fine -> coarse so every covered coarse cell holds the average
    % of its children, refreshing ghosts as each level becomes consistent.
    for kk = H.nlev:-1:2
      H = ebtTransfer(H, kk, 'restrict');
      H = ebtFillGhosts(H, kk-1);
    end

  otherwise
    error('ebt_transfer:what', 'unknown transfer "%s"', what);
end
end
