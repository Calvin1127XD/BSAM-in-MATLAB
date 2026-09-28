function H = ebTransfer(H, k, what)
%ebTransfer Inter-level data motion for the AFAS cycle.
%   ebTransfer(H, k, 'restrict')  fine solution -> QC interior + parent Q
%   ebTransfer(H, k, 'ring')      parent Q -> QC ring (gather)
%   ebTransfer(H, k, 'correct')   coarse-grid correction:
%                                  Q += P(window - QC), bilinear prolongation
%   ebTransfer(H, k, 'filldown')  cascade restrict + ghost fills k..2

switch what
  case 'restrict'
    restrict2 = H.K.restrict2;
    np = H.lev{k}.np;
    Qk = H.lev{k}.Q;  QCk = H.lev{k}.QC;  Qpar = H.lev{k-1}.Q;
    opsk = H.lev{k}.ops;
    for p = 1:np
      % 2x2 cell average of the patch interior; goes to two places: the
      % patch's own coarse shadow QC and the covered footprint in the parent.
      R = restrict2(Qk{p}(2:end-1, 2:end-1));
      qc = QCk{p};
      if isempty(qc)   % QC is released at solve end; reallocate on demand
        qc = zeros(size(R) + 2);
      end
      qc(2:end-1, 2:end-1) = R;  QCk{p} = qc;
      % Scatter the restriction into the parent-level patches it overlaps
      % (precomputed 'usc' slice pairs: source range in R -> dest in parent).
      ops = opsk{p}.usc;
      for t = 1:numel(ops)
        s = ops(t).sr; d = ops(t).dr;
        Qpar{ops(t).dp}(d(1):d(2), d(3):d(4)) = R(s(1):s(2), s(3):s(4));
      end
    end
    H.lev{k}.QC = QCk;
    H.lev{k-1}.Q = Qpar;

  case 'ring'
    np = H.lev{k}.np;
    QCk = H.lev{k}.QC;  Qpar = H.lev{k-1}.Q;
    opsk = H.lev{k}.ops;  boxk = H.lev{k}.box;
    for p = 1:np
      QC = QCk{p};
      if isempty(QC)
        % Coarse shadow is half the patch size in each direction plus a
        % one-cell ghost frame.
        b = boxk(p, :);
        QC = zeros((b(2)-b(1)+1)/2 + 2, (b(4)-b(3)+1)/2 + 2);
      end
      % Copy the ring of coarse cells surrounding this patch's footprint
      % out of the parent-level patches (linear-index gather ops).
      ops = opsk{p}.ring;
      for t = 1:numel(ops)
        QC(ops(t).dst) = Qpar{ops(t).sp}(ops(t).src);
      end
      QCk{p} = QC;
    end
    H.lev{k}.QC = QCk;

  case 'correct'
    prolong = H.K.prolong;
    np = H.lev{k}.np;
    Qk = H.lev{k}.Q;  QCk = H.lev{k}.QC;
    for p = 1:np
      % FAS correction: gather the updated coarse window over this patch,
      % subtract the pre-descent shadow QC (leaving only the coarse-grid
      % correction), prolong it bilinearly, and add to the fine interior.
      W = ebGatherWindow(H, k, p);
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
      H = ebTransfer(H, kk, 'restrict');
      H = ebFillGhosts(H, kk-1);
    end

  otherwise
    error('eb_transfer:what', 'unknown transfer "%s"', what);
end
end
