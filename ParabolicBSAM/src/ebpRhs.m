function H = ebpRhs(H, pprob, t1, Dt, scheme)
%EBPRHS Assemble the implicit step's right-hand side into FTrue on every
%   level (visible-cell values are what matter; covered cells are
%   overwritten by the FAS loading).  With op_eff = L_phys + (alpha/Dt) I
%   (the hierarchy's current operator):
%
%     be   :  f(t1)           + Uo/Dt
%     bdf2 :  f(t1)           + (2 Uo - Uoo/2)/Dt
%     cn   :  2 f(t1 - Dt/2) + (2/Dt) Uo - L_phys(t0)(u^n)
%
%   (CN written in its doubled form so the implicit operator is
%   L_phys(t1) + (2/Dt) I. The explicit OLD composite operator includes
%   coarse-fine flux corrections. Call BEFORE ebpSetTime refreshes to t1.)

dom = pprob.domain;
K = H.K;
if strcmp(scheme, 'cn')
  H = ebTransfer(H, 0, 'filldown');
  for k = 1:H.nlev
    H = ebFillGhosts(H, k);
  end
end
% sub-root levels are skipped: their FTrue is irrelevant (the FAS loading
% overwrites their full-domain footprints, and they carry no visible cells)
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  for p = 1:lev.np
    b = lev.box(p, :);
    [X, Y] = ndgrid(dom(1) + ((b(1):b(2)) - 0.5)*lev.h(1), ...
                    dom(3) + ((b(3):b(4)) - 0.5)*lev.h(2));
    Uo = H.par.Uo{k}{p};
    switch scheme
      case 'be'
        FT = pprob.f(X, Y, t1) + Uo / Dt;
      case 'bdf2'
        FT = pprob.f(X, Y, t1) + (2*Uo - 0.5*H.par.Uoo{k}{p}) / Dt;
      case 'cn'
        Le = K.op(lev.Q{p}, lev.DeW{p}, lev.DnS{p}, lev.CC{p}, ...
                  lev.h(1), lev.h(2));
        FT = 2 * pprob.f(X, Y, t1 - Dt/2) + ...
             (2/Dt + H.par.aDt) * Uo - Le;
      otherwise
        error('ebp_rhs:scheme', 'unknown scheme %s', scheme);
    end
    if isscalar(FT), FT = FT * ones(size(X)); end
    H.lev{k}.FTrue{p} = FT;
  end
end
% The conservative composite operator is L_patch - correction. Its
% explicit negative therefore ADDS each old-state correction to the RHS.
if strcmp(scheme, 'cn') && H.opts.MassCorrection
  for k = H.rootIdx+1:H.nlev
    H = ebTransfer(H, k, 'ring');
    for p = 1:H.lev{k}.np
      corr = ebFaceCorrections(H, k, p);
      for f = 1:4
        cs = corr{f};
        for j = 1:numel(cs.ops)
          op = cs.ops(j);
          FT = H.lev{k-1}.FTrue{op.dp};
          FT(op.dst) = FT(op.dst) + cs.vals(op.src);
          H.lev{k-1}.FTrue{op.dp} = FT;
        end
      end
    end
  end
end
end
