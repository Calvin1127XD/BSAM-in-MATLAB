function stats = ebtTreeCheck(H)
%ebtTreeCheck Validate the strict single-parent tree invariants of the
%   hierarchy (the MATLAB analogue of the Fortran NewSubGrids/MakeNewGrid
%   STOP checks).  Errors out on the first violation; returns summary
%   statistics when everything holds.
%
%   Checked, for every level k >= 2 and patch p:
%     1. parent index valid; mbounds within the parent's extent
%        (single-parent containment: a child may touch but never exceed
%        its parent);
%     2. global box odd-start/even-end (even-aligned on the ratio-2
%        lattice) and consistent with parent.box + mbounds;
%     3. patches of one level have disjoint interiors (owner-map count);
%     4. children lists match the parent field, and the parent's covered
%        mask equals the union of its children's footprints;
%     5. proper nesting: every ring cell of a patch is covered by the
%        parent level's patch union or lies outside the domain (so the
%        QIF/LIF stencil always reads synchronized coarse data).
%
%   stats fields: nlev, rootIdx, perLevel (np, cells), maxChildren,
%   maxDepth (AMR levels above the root).

stats.nlev = H.nlev;
stats.rootIdx = H.rootIdx;
stats.maxChildren = 0;
stats.perLevel = zeros(H.nlev, 2);

for k = 1:H.nlev
  lev = H.lev{k};
  ncell = 0;
  for p = 1:lev.np, ncell = ncell + numel(lev.F{p}); end
  stats.perLevel(k, :) = [lev.np, ncell];

  % ---- 3. disjoint interiors via owner map
  O = zeros(lev.n(1), lev.n(2), 'uint8');
  for p = 1:lev.np
    b = lev.box(p, :);
    O(b(1):b(2), b(3):b(4)) = O(b(1):b(2), b(3):b(4)) + 1;
  end
  if any(O(:) > 1)
    error('ebt_tree_check:overlap', 'level %d has overlapping patches.', k);
  end

  if k == 1, continue; end
  levC = H.lev{k-1};

  % rebuild the coarse cover mask for the nesting check
  coverC = false(levC.n(1), levC.n(2));
  for q = 1:levC.np
    qb = levC.box(q, :);
    coverC(qb(1):qb(2), qb(3):qb(4)) = true;
  end

  childCount = zeros(levC.np, 1);
  covChk = cell(levC.np, 1);
  for q = 1:levC.np
    qb = levC.box(q, :);
    covChk{q} = false(qb(2)-qb(1)+1, qb(4)-qb(3)+1);
  end

  for p = 1:lev.np
    b = lev.box(p, :);
    q = lev.parent(p);
    mb = lev.mbounds(p, :);

    % ---- 1. lineage containment
    if q < 1 || q > levC.np || any(~isfinite(mb))
      error('ebt_tree_check:parent', 'level %d patch %d: bad parent record.', k, p);
    end
    qb = levC.box(q, :);
    pmx = qb(2)-qb(1)+1; pmy = qb(4)-qb(3)+1;
    if mb(1) < 1 || mb(2) > pmx || mb(3) < 1 || mb(4) > pmy || ...
       mb(1) > mb(2) || mb(3) > mb(4)
      error('ebt_tree_check:containment', ...
        'level %d patch %d: mbounds outside parent %d.', k, p, q);
    end

    % ---- 2. alignment + box/mbounds consistency
    if mod(b(1),2) ~= 1 || mod(b(2),2) ~= 0 || mod(b(3),2) ~= 1 || mod(b(4),2) ~= 0
      error('ebt_tree_check:align', 'level %d patch %d: box not 2-aligned.', k, p);
    end
    fp = [(b(1)+1)/2, b(2)/2, (b(3)+1)/2, b(4)/2];
    if ~isequal(mb, [fp(1)-qb(1)+1, fp(2)-qb(1)+1, fp(3)-qb(3)+1, fp(4)-qb(3)+1])
      error('ebt_tree_check:mbounds', ...
        'level %d patch %d: mbounds inconsistent with global box.', k, p);
    end

    childCount(q) = childCount(q) + 1;
    covChk{q}(mb(1):mb(2), mb(3):mb(4)) = true;

    % ---- 5. proper nesting of the ring (one coarse cell around fp)
    ri = max(fp(1)-1, 1) : min(fp(2)+1, levC.n(1));
    rj = max(fp(3)-1, 1) : min(fp(4)+1, levC.n(2));
    if ~all(coverC(ri, rj), 'all')
      error('ebt_tree_check:nesting', ...
        'level %d patch %d: ring cell not covered by level %d.', k, p, k-1);
    end
  end

  % ---- 4. children lists and covered masks
  for q = 1:levC.np
    kids = levC.children{q};
    if numel(kids) ~= childCount(q) || (~isempty(kids) ...
        && ~isequal(sort(lev.parent(kids)).', q*ones(1, numel(kids))))
      error('ebt_tree_check:children', ...
        'level %d patch %d: children list inconsistent.', k-1, q);
    end
    if ~isequal(levC.covered{q}, covChk{q})
      error('ebt_tree_check:covered', ...
        'level %d patch %d: covered mask mismatch.', k-1, q);
    end
    stats.maxChildren = max(stats.maxChildren, numel(kids));
  end
end

stats.maxDepth = H.nlev - H.rootIdx + 1;
end
