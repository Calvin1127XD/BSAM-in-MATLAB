function H = ebtSetup(prob, opts)
%ebtSetup Build the initial multigrid hierarchy (sub-root chain + root).
%
%   The hierarchy H.lev{1..nlev} is one uniform ladder: levels 1..rootIdx-1
%   are full-domain uniform grids used only as multigrid coarse levels;
%   level rootIdx is the root (BaseCells); adaptive levels are appended
%   above by ebtAddLevel.  Every level uses the same patch data layout,
%   so the AFAS V-cycle code makes no distinction between them.
%
%   Tree lineage: this mirrors the Fortran BSAM "below-seed" chain
%   (CreateBelowSeedLevels): each sub-root level holds exactly ONE patch
%   covering the whole domain, whose parent is the single patch of the
%   level below, with mbounds spanning that parent entirely.  The root
%   patch of level rootIdx is the (single) tree root; adaptive patches
%   above it each carry exactly one parent.

H.prob = prob;
H.opts = opts;
H.K = ebtKernels();

% ---- sub-root size chain (refinement ratio 2 throughout)
% Halve the root grid while it stays even, at least 8 wide, and above
% CoarsestCells; the chain (coarsest first) becomes levels 1..rootIdx.
n0 = opts.BaseCells;
sizes = n0;
ns = n0;
while all(mod(ns, 2) == 0) && min(ns)/2 >= 4 && min(ns) > opts.CoarsestCells
  ns = ns / 2;
  sizes = [ns; sizes]; %#ok<AGROW>
end
H.rootIdx = size(sizes, 1);
H.nlev = H.rootIdx;

for k = 1:H.nlev
  % whole-domain patch; parent = the single patch of the level below
  par = double(k > 1);                     % 0 for level 1, else patch 1
  H.lev{k} = ebtMakeLevel(H, k, sizes(k, :), [1 sizes(k,1) 1 sizes(k,2)], par);
end

% ---- coefficient fill-down (restrict face D / cell C down the ladder)
H = ebtCoeffFilldown(H);

% ---- transfer operators between consecutive levels
for k = 2:H.nlev
  H = ebtBuildOps(H, k);
end
H = ebtGatherParentCoeffs(H);

% ---- initial guess + ghosts
ig = opts.InitialGuess;
if ~isempty(ig)
  dom = prob.domain;
  for k = 1:H.nlev
    lev = H.lev{k};
    for p = 1:lev.np
      b = lev.box(p, :);
      xc = dom(1) + ((b(1):b(2)) - 0.5) * lev.h(1);
      yc = dom(3) + ((b(3):b(4)) - 0.5) * lev.h(2);
      [XC, YC] = ndgrid(xc, yc);
      H.lev{k}.Q{p}(2:end-1, 2:end-1) = ig(XC, YC);
    end
  end
end
for k = 1:H.nlev
  H = ebtFillGhosts(H, k);
end

% ---- coarsest-level direct solver
H = ebtCoarsest(H, 'assemble');

% singular pure-Neumann: project the source onto the compatible subspace
% (no-op otherwise; see ebtCompatProject for the multi-level stall this
% prevents)
H = ebtCompatProject(H);

H.stats = struct();
H.stats.stages = {};
end
