function H = ebCoarseLoad(H, k)
%ebCoarseLoad Build the FAS coarse-level force on level k-1 from level k
%   (paper step a4, eq. 4.25):
%     on covered cells:   f_{k-1} = R(f_k - L_k u_k) + L_{k-1}(R u_k)
%     adjacent outside:   f_{k-1} += C(u)  (flux-balance mass correction,
%                                           eq. 4.20, when enabled)
%   Assumes: level-k ghosts current, restriction already pushed
%   (ebTransfer 'restrict'), parent ghosts filled, QC ring gathered.
%   Parent coefficients over the footprint are materialized transiently
%   (ebParentCoeffs) -- they are not stored on the patch.

lev = H.lev{k};
hx = lev.h(1); hy = lev.h(2);
hcx = H.lev{k-1}.h(1); hcy = H.lev{k-1}.h(2);
withCorr = H.opts.MassCorrection;

for p = 1:lev.np
  % Fine residual f - L u, restricted, plus the coarse operator applied to
  % the restricted solution (the FAS "tau" combination on covered cells).
  Q = lev.Q{p};
  RF = lev.F{p} - H.K.op(Q, lev.DeW{p}, lev.DnS{p}, lev.CC{p}, hx, hy);
  [PDeW, PDnS, PCC] = ebParentCoeffs(H, k, p);
  loadC = H.K.restrict2(RF) ...
        + H.K.op(lev.QC{p}, PDeW, PDnS, PCC, hcx, hcy);

  % Scatter the coarse load into the parent patches' F over the footprint.
  ops = lev.ops{p};
  for t = 1:numel(ops.fsc)
    s = ops.fsc(t).sr; d = ops.fsc(t).dr;
    H.lev{k-1}.F{ops.fsc(t).dp}(d(1):d(2), d(3):d(4)) = ...
      loadC(s(1):s(2), s(3):s(4));
  end

  % Flux-balance correction: add the fine/coarse flux mismatch to the
  % parent cells just OUTSIDE the footprint, so the composite scheme
  % conserves mass across the coarse-fine interface.
  if ~withCorr, continue; end
  corr = ebFaceCorrections(H, k, p);
  for f = 1:4
    cs = corr{f};
    for t = 1:numel(cs.ops)
      H.lev{k-1}.F{cs.ops(t).dp}(cs.ops(t).dst) = ...
        H.lev{k-1}.F{cs.ops(t).dp}(cs.ops(t).dst) + cs.vals(cs.ops(t).src);
    end
  end
end
end
