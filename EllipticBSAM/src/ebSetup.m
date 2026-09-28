function H = ebSetup(prob, opts)
%ebSetup Build the initial multigrid hierarchy (sub-root chain + root).
%
%   The hierarchy H.lev{1..nlev} is one uniform ladder: levels 1..rootIdx-1
%   are full-domain uniform grids used only as multigrid coarse levels;
%   level rootIdx is the root (BaseCells); adaptive levels are appended
%   above by ebAddLevel.  Every level uses the same patch data layout,
%   so the AFAS V-cycle code makes no distinction between them.

H.prob = prob;
H.opts = opts;
H.K = ebKernels();

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
  H.lev{k} = ebMakeLevel(H, k, sizes(k, :), [1 sizes(k,1) 1 sizes(k,2)]);
end

% ---- coefficient fill-down (restrict face D / cell C down the ladder)
H = ebCoeffFilldown(H);

% ---- transfer operators between consecutive levels
for k = 2:H.nlev
  H = ebBuildOps(H, k);
end
H = ebGatherParentCoeffs(H);

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
  H = ebFillGhosts(H, k);
end

% ---- coarsest-level direct solver
H = ebCoarsest(H, 'assemble');

% singular pure-Neumann: project the source onto the compatible subspace
% (no-op otherwise; see ebCompatProject for the multi-level stall this
% prevents)
H = ebCompatProject(H);

H.stats = struct();
H.stats.stages = {};
end
