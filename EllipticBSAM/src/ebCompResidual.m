function [res, H] = ebCompResidual(H)
%ebCompResidual Residual of the COMPOSITE discretization, measured on
%   visible cells (cells not covered by a finer level), from the root up.
%   On coarse cells adjacent to a coarse-fine interface the flux-balance
%   correction is added to the force (when MassCorrection is on), so this
%   is the residual of the conservative composite system whose fixed point
%   the V-cycle computes (paper Theorems 3.11 / 4.3).
%
%   Memory note: the per-patch residual array is created, normed and
%   discarded one patch at a time; the corrections are pre-collected as
%   perimeter-sized (dst, vals) chunks addressed to the receiving patch,
%   applied sequentially (chunk-internal indices are unique, so repeated
%   targets across chunks -- reentrant corners -- accumulate correctly).

r0 = H.rootIdx;
K = H.K;

chunks = cell(1, H.nlev);
for k = r0:H.nlev
  chunks{k} = cell(1, H.lev{k}.np);
end

if H.opts.MassCorrection
  for k = max(2, r0+1):H.nlev
    H = ebTransfer(H, k, 'ring');
    % Refresh the QC interior (= R u_fine) as well: the corrections are
    % functions of the state (u_f, R u_f, u_c), so a stale interior --
    % zeros right after ebAddLevel, or after the solve-end QC release --
    % would make them inconsistent and inflate the measured residual by
    % orders of magnitude.  On the normal cycle path ('filldown' just ran)
    % this rewrites the identical values, so the V-cycle is unaffected.
    QCk = H.lev{k}.QC;  Qk = H.lev{k}.Q;
    for p = 1:H.lev{k}.np
      qc = QCk{p};
      qc(2:end-1, 2:end-1) = K.restrict2(Qk{p}(2:end-1, 2:end-1));
      QCk{p} = qc;
    end
    H.lev{k}.QC = QCk;
  end
  for k = r0+1:H.nlev
    for p = 1:H.lev{k}.np
      corr = ebFaceCorrections(H, k, p);
      for f = 1:4
        cs = corr{f};
        for t = 1:numel(cs.ops)
          chunks{k-1}{cs.ops(t).dp}{end+1} = struct( ...
            'dst', cs.ops(t).dst, 'vals', cs.vals(cs.ops(t).src)); %#ok<AGROW>
        end
      end
    end
  end
end

linf = 0; l2sq = 0; fnorm = 0; ncell = 0;
for k = r0:H.nlev
  lev = H.lev{k};
  area = lev.h(1) * lev.h(2);
  for p = 1:lev.np
    % Plain per-level residual, then apply the pre-collected flux-balance
    % chunks addressed to this patch (from the finer level's interfaces).
    r = lev.FTrue{p} - K.op(lev.Q{p}, lev.DeW{p}, lev.DnS{p}, ...
                            lev.CC{p}, lev.h(1), lev.h(2));
    ch = chunks{k}{p};
    for c = 1:numel(ch)
      r(ch{c}.dst) = r(ch{c}.dst) + ch{c}.vals;
    end
    % Norm only over visible cells (covered cells belong to finer levels).
    vis = ~lev.covered{p};
    e = r(vis);
    if ~isempty(e)
      linf = max(linf, max(abs(e)));
      l2sq = l2sq + sum(e.^2) * area;
      fv = lev.FTrue{p}(vis);
      fnorm = max(fnorm, max(abs(fv)));
      ncell = ncell + numel(e);
    end
  end
end

res.linf = linf;
res.l2 = sqrt(l2sq);
res.fnorm = max(fnorm, realmin);
res.rel = linf / res.fnorm;
res.ncell = ncell;
end
