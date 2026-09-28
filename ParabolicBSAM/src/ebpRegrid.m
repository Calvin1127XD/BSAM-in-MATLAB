function Hn = ebpRegrid(Hold, pprob, t, aDtNext, opts)
%EBPREGRID Rebuild the hierarchy around the current solution and
%   transfer the state CONSERVATIVELY:
%     - the root is copied (identical grids);
%     - each new level is initialized by conservative-linear prolongation
%       from its (already transferred) parent, then overwritten by direct
%       copies wherever the OLD hierarchy carried the same resolution;
%     - regions that lost resolution inherit the old fill-down values,
%       i.e. conservative averages of the old fine data.
%   Net effect: the composite "mass" sum is preserved exactly across the
%   regrid. The BDF2 history Uoo rides the same transfer (conservative-
%   linear prolongation with one-sided slopes where freshly refined).
%   Tagging uses the transferred data progressively, so the same
%   indicator machinery (ulap/rte/fun) drives moving refinement.

probE = ebpEprob(pprob, t, aDtNext);
Hn = ebSetup(probE, opts.eb);
haveUoo = ~isempty(Hold.par.Uoo);

% ---- root transfer (single full-domain patch on both sides)
r = Hn.rootIdx;
Hn.lev{r}.Q{1}(2:end-1, 2:end-1) = Hold.lev{r}.Q{1}(2:end-1, 2:end-1);
UooN = cell(1, Hn.nlev);
if haveUoo
  UooN{r}{1} = Hold.par.Uoo{r}{1};
end
% sub-root levels: restriction of the root (consistent coarse states)
for k = r-1:-1:1
  Hn.lev{k}.Q{1}(2:end-1, 2:end-1) = ...
    Hn.K.restrict2(Hn.lev{k+1}.Q{1}(2:end-1, 2:end-1));
end
for k = 1:r
  Hn = ebFillGhosts(Hn, k);
end

% ---- progressive refinement with transferred data
while (Hn.nlev - Hn.rootIdx + 1) < opts.eb.MaxLevels
  boxes = ebTagCluster(Hn);
  if isempty(boxes), break; end
  Hn = ebAddLevel(Hn, boxes);
  kN = Hn.nlev;
  [Hn, UooN] = transferLevel(Hn, Hold, kN, UooN, haveUoo);
  Hn = ebFillGhosts(Hn, kN);
end
Hn = ebTransfer(Hn, 0, 'filldown');

% ---- history
Hn.par = struct();
Hn.par.Uo = interiors(Hn);
if haveUoo
  % fill any levels beyond the loop (none expected) defensively
  Hn.par.Uoo = UooN;
else
  Hn.par.Uoo = [];
end
Hn.par.stateKey = [t * (pprob.tdD || pprob.tdC || pprob.tdBC), aDtNext];
Hn.par.aDt = aDtNext;
Hn.stats = Hold.stats;          % keep the accumulated stage history
end

% ========================================================================
function [Hn, UooN] = transferLevel(Hn, Hold, kN, UooN, haveUoo)
lev = Hn.lev{kN};
for p = 1:lev.np
  % conservative-linear prolongation from the transferred parent level
  W = ebGatherWindow(Hn, kN, p);
  Qf = conslin(W);
  if haveUoo
    Pf = parentFootprint(Hn, kN, p, UooN{kN-1});
    Uf = conslinInterior(Pf);   % conservative-linear (one-sided at edges)
  end
  % overwrite with old same-level data wherever it exists
  b = lev.box(p, :);
  if kN <= Hold.nlev
    for q = 1:Hold.lev{kN}.np
      qb = Hold.lev{kN}.box(q, :);
      r = [max(b(1), qb(1)), min(b(2), qb(2)), ...
           max(b(3), qb(3)), min(b(4), qb(4))];
      if r(1) > r(2) || r(3) > r(4), continue; end
      di = r(1)-b(1)+1:r(2)-b(1)+1;
      dj = r(3)-b(3)+1:r(4)-b(3)+1;
      si = r(1)-qb(1)+1:r(2)-qb(1)+1;
      sj = r(3)-qb(3)+1:r(4)-qb(3)+1;
      OldQ = Hold.lev{kN}.Q{q};
      Qf(di, dj) = OldQ(si+1, sj+1);
      if haveUoo
        Uf(di, dj) = Hold.par.Uoo{kN}{q}(si, sj);
      end
    end
  end
  Hn.lev{kN}.Q{p}(2:end-1, 2:end-1) = Qf;
  if haveUoo
    UooN{kN}{p} = Uf;
  end
end
end

% ------------------------------------------------------------------------
function F = conslin(W)
%CONSLIN Conservative-linear prolongation of a padded coarse window:
%   per coarse cell, value + centered slopes; the four children average
%   exactly to the parent value (mass-preserving, 2nd order).
C0 = W(2:end-1, 2:end-1);
sx = (W(3:end, 2:end-1) - W(1:end-2, 2:end-1)) / 2;
sy = (W(2:end-1, 3:end) - W(2:end-1, 1:end-2)) / 2;
[cm, cn] = size(C0);
F = zeros(2*cm, 2*cn);
F(1:2:end, 1:2:end) = C0 - sx/4 - sy/4;
F(2:2:end, 1:2:end) = C0 + sx/4 - sy/4;
F(1:2:end, 2:2:end) = C0 - sx/4 + sy/4;
F(2:2:end, 2:2:end) = C0 + sx/4 + sy/4;
end

% ------------------------------------------------------------------------
function F = conslinInterior(P)
%CONSLININTERIOR Conservative-linear prolongation of an interior-only
%   coarse array: centered slopes inside, one-sided at the edges (still
%   exactly conservative; used for the history field, which carries no
%   ghost ring).
[cm, cn] = size(P);
sx = zeros(cm, cn); sy = zeros(cm, cn);
if cm > 1
  sx(2:end-1, :) = (P(3:end, :) - P(1:end-2, :)) / 2;
  sx(1, :) = P(2, :) - P(1, :);
  sx(end, :) = P(end, :) - P(end-1, :);
end
if cn > 1
  sy(:, 2:end-1) = (P(:, 3:end) - P(:, 1:end-2)) / 2;
  sy(:, 1) = P(:, 2) - P(:, 1);
  sy(:, end) = P(:, end) - P(:, end-1);
end
F = zeros(2*cm, 2*cn);
F(1:2:end, 1:2:end) = P - sx/4 - sy/4;
F(2:2:end, 1:2:end) = P + sx/4 - sy/4;
F(1:2:end, 2:2:end) = P - sx/4 + sy/4;
F(2:2:end, 2:2:end) = P + sx/4 + sy/4;
end

% ------------------------------------------------------------------------
function P = parentFootprint(Hn, kN, p, UooParentLevel)
%PARENTFOOTPRINT Gather the parent-level history over the footprint of
%   patch p (uses the fsc rect descriptors, which address the parents'
%   UNPADDED arrays).
b = Hn.lev{kN}.box(p, :);
cmx = (b(2)-b(1)+1)/2; cmy = (b(4)-b(3)+1)/2;
P = zeros(cmx, cmy);
ops = Hn.lev{kN}.ops{p}.fsc;
for t = 1:numel(ops)
  s = ops(t).sr; d = ops(t).dr;
  P(s(1):s(2), s(3):s(4)) = UooParentLevel{ops(t).dp}(d(1):d(2), d(3):d(4));
end
end

% ------------------------------------------------------------------------
function U = interiors(H)
U = cell(1, H.nlev);
for k = 1:H.nlev
  for p = 1:H.lev{k}.np
    U{k}{p} = H.lev{k}.Q{p}(2:end-1, 2:end-1);
  end
end
end
