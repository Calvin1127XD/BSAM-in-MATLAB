function sol = ebSolve(prob, varargin)
%ebSolve Block-structured adaptive FAS multigrid solve of
%     -div(D(x,y) grad u) + C(x,y) u = f(x,y)   on [a,b]x[c,d].
%
%   sol = ebSolve(prob)                  with default options
%   sol = ebSolve(prob, opts)            opts from ebOptions
%   sol = ebSolve(prob, 'Name', value)   options inline
%
%   Mass-conservative AFAS scheme of Feng, Guo, Lowengrub & Wise,
%   J. Comput. Phys. 352 (2018) 463-497: cell-centered flux-form finite
%   differences, QIF/LIF ghost interpolation at coarse-fine interfaces,
%   flux-balanced zombie corrections to the coarse-level force, red-black
%   Gauss-Seidel smoothing, Berger-Rigoutsos regridding.
%
%   Returns sol with fields:
%     .H       full hierarchy (levels, patches, solution)
%     .stats   per-stage V-cycle histories / convergence factors
%     .err     composite error vs prob.exact (if provided)
%     .mass    discrete conservation check (ebMassCheck)
%     .levels  per-level summary (boxes, cells)

if numel(varargin) == 1 && isstruct(varargin{1})
  opts = varargin{1};
else
  opts = ebOptions(varargin{:});
end

H = ebSetup(prob, opts);

if strcmpi(opts.Indicator, 'fun')
  % geometry-driven tagging needs no intermediate solutions: build the
  % whole hierarchy first, then solve the composite system once
  while (H.nlev - H.rootIdx + 1) < opts.MaxLevels
    boxes = ebTagCluster(H);
    if isempty(boxes), break; end
    H = ebAddLevel(H, boxes);
  end
  H = ebSolveCycles(H);
else
  % solution-driven tagging: solve, tag, refine, re-solve
  H = ebSolveCycles(H);
  while (H.nlev - H.rootIdx + 1) < opts.MaxLevels
    boxes = ebTagCluster(H);
    if isempty(boxes)
      if opts.Verbose >= 1
        fprintf('[ebsam] no cells tagged; stopping at %d AMR level(s)\n', ...
          H.nlev - H.rootIdx + 1);
      end
      break;
    end
    H = ebAddLevel(H, boxes);
    H = ebSolveCycles(H);
  end
end

% pin the mean for the singular pure-Neumann case
mz = opts.MeanZero;
if (ischar(mz) || isstring(mz)) && strcmpi(mz, 'auto')
  mz = H.coarsest.singular;
end
if mz
  [m, ~] = ebCompositeMean(H);
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
  sol.err = ebErrorNorms(H);
else
  sol.err = [];
end
sol.mass = ebMassCheck(H);

for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  ncell = 0;
  for p = 1:lev.np, ncell = ncell + numel(lev.F{p}); end
  sol.levels(k - H.rootIdx + 1) = struct('n', lev.n, 'h', lev.h, ...
    'np', lev.np, 'cells', ncell, 'box', lev.box);
end
end
