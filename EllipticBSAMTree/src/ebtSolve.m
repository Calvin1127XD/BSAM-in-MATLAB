function sol = ebtSolve(prob, varargin)
%ebtSolve Tree-based block-structured adaptive FAS multigrid solve of
%     -div(D(x,y) grad u) + C(x,y) u = f(x,y)   on [a,b]x[c,d].
%
%   sol = ebtSolve(prob)                  with default options
%   sol = ebtSolve(prob, opts)            opts from ebtOptions
%   sol = ebtSolve(prob, 'Name', value)   options inline
%
%   Mass-conservative AFAS scheme of Feng, Guo, Lowengrub & Wise,
%   J. Comput. Phys. 352 (2018) 463-497, organized as a STRICT PATCH TREE
%   in the style of the Fortran BSAM 2.0 library (treeops.f90): tags are
%   clustered per parent patch, every child patch belongs to exactly one
%   parent that fully contains it, and all coarse-fine data motion
%   (window, ring, restriction, FAS loading, correction) is a direct
%   slice exchange with that single parent.
%
%   Returns sol with fields:
%     .H       full hierarchy (levels, patches, solution, lineage)
%     .stats   per-stage V-cycle histories / convergence factors
%     .err     composite error vs prob.exact (if provided)
%     .mass    discrete conservation check (ebtMassCheck)
%     .levels  per-level summary (boxes, cells)

if numel(varargin) == 1 && isstruct(varargin{1})
  opts = varargin{1};
else
  opts = ebtOptions(varargin{:});
end

H = ebtSetup(prob, opts);

if strcmpi(opts.Indicator, 'fun')
  % geometry-driven tagging needs no intermediate solutions: build the
  % whole hierarchy first, then solve the composite system once
  while (H.nlev - H.rootIdx + 1) < opts.MaxLevels
    boxes = ebtTagCluster(H);
    if isempty(boxes), break; end
    H = ebtAddLevel(H, boxes);
  end
  H = ebtSolveCycles(H);
else
  % solution-driven tagging: solve, tag, refine, re-solve
  H = ebtSolveCycles(H);
  while (H.nlev - H.rootIdx + 1) < opts.MaxLevels
    boxes = ebtTagCluster(H);
    if isempty(boxes)
      if opts.Verbose >= 1
        fprintf('[ebtsam] no cells tagged; stopping at %d AMR level(s)\n', ...
          H.nlev - H.rootIdx + 1);
      end
      break;
    end
    H = ebtAddLevel(H, boxes);
    H = ebtSolveCycles(H);
  end
end

% pin the mean for the singular pure-Neumann case
mz = opts.MeanZero;
if (ischar(mz) || isstring(mz)) && strcmpi(mz, 'auto')
  mz = H.coarsest.singular;
end
if mz
  [m, ~] = ebtCompositeMean(H);
  for k = 1:H.nlev
    for p = 1:H.lev{k}.np
      H.lev{k}.Q{p} = H.lev{k}.Q{p} - m;
    end
  end
end

% release the per-cycle restriction buffers (reallocated on demand if the
% hierarchy is solved again); the solved state keeps only Q, sources,
% coefficients, diagonals and O(perimeter) metadata
for k = 2:H.nlev
  for p = 1:H.lev{k}.np
    H.lev{k}.QC{p} = [];
  end
end

sol.H = H;
sol.stats = H.stats;
if ~isempty(H.prob.exact)
  sol.err = ebtErrorNorms(H);
else
  sol.err = [];
end
sol.mass = ebtMassCheck(H);

for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  ncell = 0;
  for p = 1:lev.np, ncell = ncell + numel(lev.F{p}); end
  sol.levels(k - H.rootIdx + 1) = struct('n', lev.n, 'h', lev.h, ...
    'np', lev.np, 'cells', ncell, 'box', lev.box, ...
    'parent', lev.parent, 'mbounds', lev.mbounds);
end
end
