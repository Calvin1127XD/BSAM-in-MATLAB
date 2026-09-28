function H = ebtBuildOps(H, k)
%ebtBuildOps Precompute all index-based transfer operations for level k
%   (k >= 2) against its parent level k-1, and the same-level ghost
%   exchange within level k.  Everything in the V-cycle hot path then
%   reduces to vectorized sliced copies.
%
%   TREE VERSION: because every patch has exactly ONE parent that fully
%   contains it (lev.parent / lev.mbounds, see ebtMakeLevel), all
%   parent-side transfers are single contiguous slices of the parent's
%   arrays -- the (cmx+2)x(cmy+2) parent window INCLUDING its ghost ring
%   is one rectangle of the parent's padded array.  Where the ring leaves
%   the parent's interior it lands on the parent's own ghost cells, whose
%   values were synchronized at the coarser level (coarse-fine interp,
%   sibling exchange, or physical BC) -- exactly the Fortran BSAM
%   InterpolateCoarseFine contract: the child reads ONLY parent%q.
%
%   Built per patch p (parent q, mbounds [Il Ih Jl Jh] in parent cells):
%     ops.win     : padded-parent slice [Il Ih+2 Jl Jh+2] = whole window
%     ops.ring    : 4 rect copies parent Q -> QC perimeter (W/E/S/N)
%     ops.usc     : restricted fine solution -> parent Q interior (1 rect)
%     ops.fsc     : FAS coarse loading -> parent F (1 rect)
%     ops.pcc/pdew/pdns : parent coefficient windows (1 rect each)
%     ops.corrW/E/S/N   : flux-balance mass-correction targets (scatter-
%                         add into level-(k-1) F just outside the
%                         footprint; targets may belong to the parent's
%                         SIBLINGS when the child touches the parent's
%                         edge -- resolved by owner map, the precomputed
%                         analogue of the Fortran masscorrlayer strips)
%   plus landscape codes (1 phys, 2 same-level, 3 coarse-fine) for every
%   ghost cell, corner classes, and the same-level exchange list.

levF = H.lev{k};
levC = H.lev{k-1};
nF = levF.n; nC = levC.n;

% owner maps (same-level neighbor resolution + correction routing)
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
  q  = levF.parent(p);
  mb = levF.mbounds(p, :);                 % [Il Ih Jl Jh], parent-local
  Il = mb(1); Ih = mb(2); Jl = mb(3); Jh = mb(4);
  cmx = Ih - Il + 1; cmy = Jh - Jl + 1;
  if cmx ~= mx/2 || cmy ~= my/2
    error('ebt_build_ops:mbounds', 'mbounds inconsistent with patch size.');
  end

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

  % ---------------- single-parent window (interior + ring, one slice) ---
  % Parent cell (i,j) lives at padded index (i+1, j+1); the window spans
  % parent cells Il-1 .. Ih+1 (0 and pmx+1 are the parent's ghost cells),
  % i.e. padded rows Il .. Ih+2 -- always within the parent's padded array
  % because 1 <= Il and Ih <= pmx (containment, ebtMakeLevel).
  ops.sp  = q;                              % the ONE parent patch
  ops.win = [Il, Ih+2, Jl, Jh+2];           % padded-parent source range

  % ring: perimeter of the window, as 4 sliced copies parent Q -> QC.
  % W/E columns include the corners; S/N rows cover the remainder.
  ops.ring = struct( ...
    'sr', {[Il,    Il,    Jl,   Jh+2], ...   % west column
           [Ih+2,  Ih+2,  Jl,   Jh+2], ...   % east column
           [Il+1,  Ih+1,  Jl,   Jl  ], ...   % south row (no corners)
           [Il+1,  Ih+1,  Jh+2, Jh+2]}, ...  % north row (no corners)
    'dr', {[1,     1,     1,    cmy+2], ...
           [cmx+2, cmx+2, 1,    cmy+2], ...
           [2,     cmx+1, 1,    1    ], ...
           [2,     cmx+1, cmy+2, cmy+2]});

  % ---------------- scatters to the parent (single rects) ---------------
  % usc: QC interior (restricted fine solution) -> parent padded interior
  ops.usc = struct('sr', [1, cmx, 1, cmy], 'dr', [Il+1, Ih+1, Jl+1, Jh+1]);
  % fsc: FAS coarse load -> parent F (unpadded)
  ops.fsc = struct('sr', [1, cmx, 1, cmy], 'dr', [Il, Ih, Jl, Jh]);

  % ---------------- parent coefficient windows (single slices) ----------
  % Cell C over the footprint; x-faces Il..Ih+1; y-faces Jl..Jh+1.
  ops.pcc  = struct('sr', [Il, Ih,   Jl, Jh  ]);
  ops.pdew = struct('sr', [Il, Ih+1, Jl, Jh  ]);
  ops.pdns = struct('sr', [Il, Ih,   Jl, Jh+1]);

  % ---------------- mass-correction scatter targets ----------------------
  % One entry per coarse interface face; the target coarse cell just
  % outside the footprint may belong to the parent OR to one of the
  % parent's siblings (when the child touches the parent's edge); the
  % owner map routes it, like the Fortran global masscorrlayer strips.
  % land is per FINE ghost cell; sampling 2:2:end takes one flag per
  % coarse face (both fine subfaces of a coarse face share the class).
  fpg = [(ilo+1)/2, ihi/2, (jlo+1)/2, jhi/2];   % footprint, level-(k-1) global
  ops.corrW = corrFace(levC, OC, land.W(2:2:end) == 3, fpg(1)-1, (fpg(3):fpg(4)).');
  ops.corrE = corrFace(levC, OC, land.E(2:2:end) == 3, fpg(2)+1, (fpg(3):fpg(4)).');
  ops.corrS = corrFace(levC, OC, land.S(2:2:end) == 3, (fpg(1):fpg(2)).', fpg(3)-1);
  ops.corrN = corrFace(levC, OC, land.N(2:2:end) == 3, (fpg(1):fpg(2)).', fpg(4)+1);

  H.lev{k}.ops{p} = ops;
end

% ---------------- children bookkeeping + covered masks -------------------
% Direct tree lineage: children{q} = the patches whose parent is q; the
% covered mask of q is the union of its children's footprints (children
% never overlap and lie fully inside q, so plain rect stamping suffices).
for q = 1:levC.np
  H.lev{k-1}.children{q} = [];
  qb = levC.box(q, :);
  H.lev{k-1}.covered{q} = false(qb(2)-qb(1)+1, qb(4)-qb(3)+1);
end
for p = 1:levF.np
  q = levF.parent(p);
  mb = levF.mbounds(p, :);
  H.lev{k-1}.children{q}(end+1) = p;
  H.lev{k-1}.covered{q}(mb(1):mb(2), mb(3):mb(4)) = true;
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
function cs = corrFace(levC, OC, mask, gis, gjs)
%corrFace Mass-correction scatter for one face.  mask(cmy|cmx x 1) marks
%   coarse interface faces; target coarse cell is (gis(j), gjs) or
%   (gis, gjs(j)) just OUTSIDE the child footprint, owned by the parent
%   or one of its siblings.
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
    error('ebt_build_ops:corr', 'mass-correction target cell uncovered.');
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
