function H = ebtCoeffFilldown(H)
%ebtCoeffFilldown Restrict coefficients down the hierarchy so that every
%   coarse cell/face covered by a finer level carries the averaged fine
%   coefficient (paper eq. 4.15 compatibility): face D = mean of the two
%   fine faces, cell C = mean of the four fine cells.  Recomputes the
%   relaxation diagonals afterwards.
%
%   TREE VERSION: each patch writes its coarsened coefficients into its
%   SINGLE parent over mbounds -- three sliced assignments, no
%   intersection search (Fortran FillDown over the parent pointer).

for k = H.nlev:-1:2
  levF = H.lev{k};
  levC = H.lev{k-1};

  % owner map of the parent level, for routing shared-boundary faces
  OC = zeros(levC.n(1), levC.n(2), 'uint32');
  for q = 1:levC.np
    qb = levC.box(q, :);
    OC(qb(1):qb(2), qb(3):qb(4)) = q;
  end

  for p = 1:levF.np
    b = levF.box(p, :);
    mx = b(2) - b(1) + 1; my = b(4) - b(3) + 1;
    q  = levF.parent(p);
    qb = levC.box(q, :);
    pmx = qb(2) - qb(1) + 1; pmy = qb(4) - qb(3) + 1;
    mb = levF.mbounds(p, :);               % [Il Ih Jl Jh] parent-local
    Il = mb(1); Ih = mb(2); Jl = mb(3); Jh = mb(4);

    % coarsened coefficients of this patch (footprint resolution)
    cDeW = 0.5 * (levF.DeW{p}(1:2:mx+1, 1:2:my-1) + levF.DeW{p}(1:2:mx+1, 2:2:my));
    cDnS = 0.5 * (levF.DnS{p}(1:2:mx-1, 1:2:my+1) + levF.DnS{p}(2:2:mx, 1:2:my+1));
    cCC  = H.K.restrict2(levF.CC{p});

    % write into the one parent: cells, x-faces (Il..Ih+1), y-faces (Jl..Jh+1)
    H.lev{k-1}.CC{q}(Il:Ih, Jl:Jh)      = cCC;
    H.lev{k-1}.DeW{q}(Il:Ih+1, Jl:Jh)   = cDeW;
    H.lev{k-1}.DnS{q}(Il:Ih, Jl:Jh+1)   = cDnS;

    % ---- shared-face duplication fix ----------------------------------
    % A physical face on the boundary between two coarse patches is
    % stored TWICE (in both patches' face arrays).  When the child's
    % footprint touches the parent's edge, the averaged face value must
    % also land in the abutting sibling's copy, or the two sides of that
    % face would carry different D and the flux telescoping (discrete
    % conservation) breaks.  Route per face cell through the owner map,
    % since several siblings may share the segment.
    if Il == 1 && qb(1) > 1              % child west edge on parent west edge
      H = pushFaceX(H, k-1, levC, OC, qb(1), qb(3)+Jl-1 : qb(3)+Jh-1, ...
                    qb(1)-1, cDeW(1, :));
    end
    if Ih == pmx && qb(2) < levC.n(1)    % east edge
      H = pushFaceX(H, k-1, levC, OC, qb(2)+1, qb(3)+Jl-1 : qb(3)+Jh-1, ...
                    qb(2)+1, cDeW(end, :));
    end
    if Jl == 1 && qb(3) > 1              % south edge
      H = pushFaceY(H, k-1, levC, OC, qb(1)+Il-1 : qb(1)+Ih-1, qb(3), ...
                    qb(3)-1, cDnS(:, 1));
    end
    if Jh == pmy && qb(4) < levC.n(2)    % north edge
      H = pushFaceY(H, k-1, levC, OC, qb(1)+Il-1 : qb(1)+Ih-1, qb(4)+1, ...
                    qb(4)+1, cDnS(:, end));
    end
  end
end

% refresh relaxation diagonals everywhere (with the boundary fold)
for k = 1:H.nlev
  lev = H.lev{k};
  for p = 1:lev.np
    ID = H.K.invdiag(lev.DeW{p}, lev.DnS{p}, lev.CC{p}, lev.h(1), lev.h(2));
    H.lev{k}.InvDiagS{p} = 1 ./ (1 ./ ID + ebtBcFold(H.lev{k}, p));
  end
end
end

% ========================================================================
function H = pushFaceX(H, kc, levC, OC, faceG, jGs, nbrI, vals)
%pushFaceX Copy averaged x-face values into the abutting neighbors' DeW.
%   faceG : global x-face index (face g = west face of cell g) of the
%           shared boundary; jGs = global cell rows of the segment;
%   nbrI  : global cell column on the OTHER side of the face;
%   vals  : 1 x numel(jGs) averaged values (already written to the parent).
%   Cells with no owner (a genuine coarse-fine boundary of level kc) keep
%   a single copy and are skipped.
own = zeros(numel(jGs), 1);
for t = 1:numel(jGs)
  own(t) = OC(nbrI, jGs(t));
end
for s = unique(own(own > 0)).'
  ii = find(own == s).';
  sb = levC.box(s, :);
  fl = faceG - sb(1) + 1;                 % local x-face index in patch s
  H.lev{kc}.DeW{s}(fl, jGs(ii) - sb(3) + 1) = vals(ii);
end
end

% ------------------------------------------------------------------------
function H = pushFaceY(H, kc, levC, OC, iGs, faceG, nbrJ, vals)
%pushFaceY Copy averaged y-face values into the abutting neighbors' DnS.
own = zeros(numel(iGs), 1);
for t = 1:numel(iGs)
  own(t) = OC(iGs(t), nbrJ);
end
for s = unique(own(own > 0)).'
  ii = find(own == s);
  sb = levC.box(s, :);
  fl = faceG - sb(3) + 1;                 % local y-face index in patch s
  H.lev{kc}.DnS{s}(iGs(ii) - sb(1) + 1, fl) = vals(ii);
end
end
