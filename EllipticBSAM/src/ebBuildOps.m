function H = ebBuildOps(H, k)
%ebBuildOps Precompute all index-based transfer operations for level k
%   (k >= 2) against its parent level k-1, and the same-level ghost
%   exchange within level k.  Everything in the V-cycle hot path then
%   reduces to vectorized indexed copies.
%
%   Built per patch p:
%     ops.ring : parent Q -> QC ring cells (gather)
%     ops.winrect : parent Q -> window interior rect ranges (gather)
%     ops.usc  : restricted fine solution -> parent Q interior (scatter)
%     ops.fsc  : FAS coarse loading -> parent F (scatter)
%     ops.corrW/E/S/N : flux-balance mass-correction targets (scatter-add
%                       into parent F just outside the footprint)
%   plus landscape codes (1 phys, 2 same-level, 3 coarse-fine) for every
%   ghost cell, corner classes, and the same-level exchange list.

levF = H.lev{k};
levC = H.lev{k-1};
nF = levF.n; nC = levC.n;

% Owner maps: for every level-global cell, which patch (if any) owns it.
% These turn all neighbor/coverage questions below into O(1) lookups.
OF = zeros(nF(1), nF(2), 'uint32');
for p = 1:levF.np
  b = levF.box(p, :);
  OF(b(1):b(2), b(3):b(4)) = p;
end
OC = zeros(nC(1), nC(2), 'uint32');
for q = 1:levC.np
  b = levC.box(q, :);
  OC(b(1):b(2), b(3):b(4)) = q;
end

for p = 1:levF.np
  b = levF.box(p, :);
  ilo = b(1); ihi = b(2); jlo = b(3); jhi = b(4);
  mx = ihi - ilo + 1; my = jhi - jlo + 1;
  % A patch must start on an odd fine index and end on an even one so its
  % footprint is a whole number of parent cells.
  if mod(ilo,2) ~= 1 || mod(ihi,2) ~= 0 || mod(jlo,2) ~= 1 || mod(jhi,2) ~= 0
    error('eb_build_ops:align', 'patch not aligned to parent cells.');
  end
  % Footprint of the patch in parent-level (coarse) indices.
  Il = (ilo+1)/2; Ih = ihi/2; Jl = (jlo+1)/2; Jh = jhi/2;
  cmx = Ih - Il + 1; cmy = Jh - Jl + 1;

  % ---------------- landscape (per ghost cell) + same-level exchange ----
  % Classify every ghost cell: 1 = physical boundary, 2 = same-level
  % sibling supplies it, 3 = coarse-fine (interpolated from the parent).
  % Default 1 stands when the neighbor cell is outside the domain.
  land = struct('W', ones(my,1,'int8'), 'E', ones(my,1,'int8'), ...
                'S', ones(mx,1,'int8'), 'N', ones(mx,1,'int8'));
  same = struct('sp', {}, 'src', {}, 'dst', {});
  gj = (jlo:jhi).';
  gi = (ilo:ihi).';

  if ilo > 1
    nb = OF(ilo-1, jlo:jhi).';
    land.W(nb > 0) = 2;
    land.W(nb == 0) = 3;
    same = addFaceOps(same, levF, nb, ilo-1, gj, ...
      sub2ind([mx+2, my+2], ones(my,1), (2:my+1).'));
  end
  if ihi < nF(1)
    nb = OF(ihi+1, jlo:jhi).';
    land.E(nb > 0) = 2;
    land.E(nb == 0) = 3;
    same = addFaceOps(same, levF, nb, ihi+1, gj, ...
      sub2ind([mx+2, my+2], (mx+2)*ones(my,1), (2:my+1).'));
  end
  if jlo > 1
    nb = OF(ilo:ihi, jlo-1);
    land.S(nb > 0) = 2;
    land.S(nb == 0) = 3;
    same = addFaceOps(same, levF, nb, gi, jlo-1, ...
      sub2ind([mx+2, my+2], (2:mx+1).', ones(mx,1)));
  end
  if jhi < nF(2)
    nb = OF(ilo:ihi, jhi+1);
    land.N(nb > 0) = 2;
    land.N(nb == 0) = 3;
    same = addFaceOps(same, levF, nb, gi, jhi+1, ...
      sub2ind([mx+2, my+2], (2:mx+1).', (my+2)*ones(mx,1)));
  end

  % corners: SW SE NW NE
  % Each row: level-global (i,j) of the corner ghost cell and its local
  % subscripts in the padded (mx+2 x my+2) patch array.
  cdef = [ilo-1, jlo-1, 1, 1; ihi+1, jlo-1, mx+2, 1; ...
          ilo-1, jhi+1, 1, my+2; ihi+1, jhi+1, mx+2, my+2];
  cclass = ones(1, 4, 'int8');
  for c = 1:4
    cgi = cdef(c,1); cgj = cdef(c,2);
    if cgi < 1 || cgi > nF(1) || cgj < 1 || cgj > nF(2)
      cclass(c) = 1;
    else
      s = OF(cgi, cgj);
      if s > 0
        cclass(c) = 2;
        sb = levF.box(s, :);
        msx = sb(2) - sb(1) + 1; msy = sb(4) - sb(3) + 1;
        same(end+1) = struct('sp', double(s), ...
          'src', sub2ind([msx+2, msy+2], cgi-sb(1)+2, cgj-sb(3)+2), ...
          'dst', sub2ind([mx+2, my+2], cdef(c,3), cdef(c,4))); %#ok<AGROW>
      else
        cclass(c) = 3;
      end
    end
  end

  H.lev{k}.land{p} = land;
  H.lev{k}.cornerClass{p} = cclass;
  H.lev{k}.sameOps{p} = same;

  % ---------------- parent window gathers --------------------------------
  % ring cells (around footprint) - per-cell resolver
  % The ring is the one-cell frame of parent cells around the footprint:
  % full west and east columns (incl. corners), then the remaining south
  % and north rows.
  ringGi = [ (Il-1)*ones(cmy+2,1); (Ih+1)*ones(cmy+2,1); (Il:Ih).'; (Il:Ih).' ];
  ringGj = [ (Jl-1:Jh+1).';        (Jl-1:Jh+1).';        (Jl-1)*ones(cmx,1); (Jh+1)*ones(cmx,1) ];
  ops.ring = resolveCells(levC, OC, ringGi, ringGj, Il, Jl, cmx, cmy);

  % full window = ring (per-cell list) + interior rect blocks; rects are
  % stored as O(1) range descriptors and applied with sliced assignment
  ops.winrect = struct('sp', {}, 'sr', {}, 'dr', {});
  ncov = 0;
  for t = 1:numel(ops.ring), ncov = ncov + numel(ops.ring(t).dst); end
  for q = 1:levC.np
    qb = levC.box(q, :);
    r = rectIsect([Il Ih Jl Jh], qb);
    if isempty(r), continue; end
    ops.winrect(end+1) = struct('sp', q, ...
      'sr', [r(1)-qb(1)+2, r(2)-qb(1)+2, r(3)-qb(3)+2, r(4)-qb(3)+2], ...
      'dr', [r(1)-Il+2,    r(2)-Il+2,    r(3)-Jl+2,    r(4)-Jl+2]);
    ncov = ncov + (r(2)-r(1)+1) * (r(4)-r(3)+1);
  end
  if ncov ~= (cmx+2)*(cmy+2)
    error('eb_build_ops:nesting', ...
      'patch %d at level %d: parent window not fully covered (%d of %d) - proper nesting violated.', ...
      p, k, ncov, (cmx+2)*(cmy+2));
  end

  % ---------------- scatters to parent (rect ranges) ---------------------
  ops.usc = struct('dp', {}, 'sr', {}, 'dr', {});
  ops.fsc = struct('dp', {}, 'sr', {}, 'dr', {});
  for q = 1:levC.np
    qb = levC.box(q, :);
    r = rectIsect([Il Ih Jl Jh], qb);
    if isempty(r), continue; end
    sr = [r(1)-Il+1, r(2)-Il+1, r(3)-Jl+1, r(4)-Jl+1];     % in cmx x cmy
    ops.usc(end+1) = struct('dp', q, 'sr', sr, ...
      'dr', [r(1)-qb(1)+2, r(2)-qb(1)+2, r(3)-qb(3)+2, r(4)-qb(3)+2]);
    ops.fsc(end+1) = struct('dp', q, 'sr', sr, ...
      'dr', [r(1)-qb(1)+1, r(2)-qb(1)+1, r(3)-qb(3)+1, r(4)-qb(3)+1]);
  end

  % ---------------- mass-correction scatter targets ----------------------
  % west/east: one entry per coarse face j=1..cmy where the fine ghosts are
  % coarse-fine; target = parent cell just outside the footprint.
  % land is per FINE ghost cell; sampling 2:2:end takes one flag per
  % coarse face (both fine subfaces of a coarse face share the same class).
  ops.corrW = corrFace(levC, OC, land.W(2:2:end) == 3, Il-1, (Jl:Jh).');
  ops.corrE = corrFace(levC, OC, land.E(2:2:end) == 3, Ih+1, (Jl:Jh).');
  ops.corrS = corrFace(levC, OC, land.S(2:2:end) == 3, (Il:Ih).', Jl-1);
  ops.corrN = corrFace(levC, OC, land.N(2:2:end) == 3, (Il:Ih).', Jh+1);

  H.lev{k}.ops{p} = ops;
end

% ---------------- children/parents bookkeeping + covered masks -----------
% For every coarse patch record which fine patches overlap it and mark the
% overlapped cells as covered (excluded from residual norms and plots).
for q = 1:levC.np
  H.lev{k-1}.children{q} = [];
  qb = levC.box(q, :);
  cov = false(qb(2)-qb(1)+1, qb(4)-qb(3)+1);
  for p = 1:levF.np
    b = levF.box(p, :);
    % fp = the fine patch's footprint in coarse indices.
    fp = [(b(1)+1)/2, b(2)/2, (b(3)+1)/2, b(4)/2];
    r = rectIsect(fp, qb);
    if ~isempty(r)
      H.lev{k-1}.children{q}(end+1) = p;
      cov(r(1)-qb(1)+1 : r(2)-qb(1)+1, r(3)-qb(3)+1 : r(4)-qb(3)+1) = true;
    end
  end
  H.lev{k-1}.covered{q} = cov;
end
for p = 1:levF.np
  H.lev{k}.parents{p} = [];
  b = levF.box(p, :);
  fp = [(b(1)+1)/2, b(2)/2, (b(3)+1)/2, b(4)/2];
  for q = 1:levC.np
    if ~isempty(rectIsect(fp, levC.box(q, :)))
      H.lev{k}.parents{p}(end+1) = q;
    end
  end
end
end

% ========================================================================
function same = addFaceOps(same, levF, nb, gi, gj, dstAll)
%addFaceOps Append same-level copy ops for ghost cells whose neighbor
%   cell (gi,gj) is owned by a sibling patch.  gi or gj may be scalar.
idx = find(nb > 0);
if isempty(idx), return; end
if isscalar(gi), gi = gi * ones(size(gj)); end
if isscalar(gj), gj = gj * ones(size(gi)); end
sps = unique(nb(idx)).';
for s = sps
  ii = idx(nb(idx) == s);
  sb = levF.box(s, :);
  msx = sb(2) - sb(1) + 1; msy = sb(4) - sb(3) + 1;
  src = sub2ind([msx+2, msy+2], gi(ii) - sb(1) + 2, gj(ii) - sb(3) + 2);
  same(end+1) = struct('sp', double(s), 'src', src, 'dst', dstAll(ii)); %#ok<AGROW>
end
end

% ------------------------------------------------------------------------
function ops = resolveCells(levC, OC, gis, gjs, Il, Jl, cmx, cmy)
%resolveCells Build gather ops parent Q -> window for arbitrary parent
%   cells (gis,gjs), allowing one-cell overhang outside the domain (then
%   the owning patch's physical ghost cell is read).
nC = levC.n;
np = numel(gis);
sp = zeros(np, 1); src = zeros(np, 1); dst = zeros(np, 1);
for t = 1:np
  gi = gis(t); gj = gjs(t);
  gic = min(max(gi, 1), nC(1));
  gjc = min(max(gj, 1), nC(2));
  q = OC(gic, gjc);
  if q == 0
    error('eb_build_ops:ring', ...
      'parent cell (%d,%d) not covered - proper nesting violated.', gi, gj);
  end
  qb = levC.box(q, :);
  li = gi - qb(1) + 2; lj = gj - qb(3) + 2;
  mqx = qb(2) - qb(1) + 1; mqy = qb(4) - qb(3) + 1;
  if li < 1 || li > mqx+2 || lj < 1 || lj > mqy+2
    error('eb_build_ops:overhang', 'ring cell overhang exceeds ghost ring.');
  end
  sp(t) = q;
  src(t) = sub2ind([mqx+2, mqy+2], li, lj);
  dst(t) = sub2ind([cmx+2, cmy+2], gi - Il + 2, gj - Jl + 2);
end
ops = struct('sp', {}, 'src', {}, 'dst', {});
for q = unique(sp).'
  ii = (sp == q);
  ops(end+1) = struct('sp', q, 'src', src(ii), 'dst', dst(ii)); %#ok<AGROW>
end
end

% ------------------------------------------------------------------------
function cs = corrFace(levC, OC, mask, gis, gjs)
%corrFace Mass-correction scatter for one face.  mask(cmy|cmx x 1) marks
%   coarse interface faces; target parent cell is (gis(j), gjs) or
%   (gis, gjs(j)) just OUTSIDE the child footprint.
cs.mask = mask;
cs.ops = struct('dp', {}, 'src', {}, 'dst', {});
idx = find(mask);
if isempty(idx), return; end
if isscalar(gis), gis = gis * ones(size(gjs)); end
if isscalar(gjs), gjs = gjs * ones(size(gis)); end
dp = zeros(numel(idx), 1); dst = zeros(numel(idx), 1);
for t = 1:numel(idx)
  j = idx(t);
  q = OC(gis(j), gjs(j));
  if q == 0
    error('eb_build_ops:corr', 'mass-correction target cell uncovered.');
  end
  qb = levC.box(q, :);
  mqx = qb(2) - qb(1) + 1; mqy = qb(4) - qb(3) + 1;
  dp(t) = q;
  dst(t) = sub2ind([mqx, mqy], gis(j) - qb(1) + 1, gjs(j) - qb(3) + 1);
end
for q = unique(dp).'
  ii = (dp == q);
  cs.ops(end+1) = struct('dp', q, 'src', idx(ii), 'dst', dst(ii)); %#ok<AGROW>
end
end

% ------------------------------------------------------------------------
function r = rectIsect(a, b)
%rectIsect Intersection of [ilo ihi jlo jhi] boxes ([] if empty).
r = [max(a(1), b(1)), min(a(2), b(2)), max(a(3), b(3)), min(a(4), b(4))];
if r(1) > r(2) || r(3) > r(4), r = []; end
end
