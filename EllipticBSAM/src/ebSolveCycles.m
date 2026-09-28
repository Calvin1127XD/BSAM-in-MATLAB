function H = ebSolveCycles(H)
%ebSolveCycles Run V-cycles on the current hierarchy until the composite
%   residual drops below RelTol * ||f||_inf (or MaxVCycles).  Each cycle:
%   reset the working force to the true source on every level (the FAS
%   coarse loading overwrites it during the descent), run one V-cycle from
%   the finest level, re-sync covered coarse data (fill-down), measure.

opts = H.opts;

% Make every level's ghost layer consistent with the current solution
% before measuring anything (physical BCs + same-level seams + CF interp).
for k = 1:H.nlev
  H = ebFillGhosts(H, k);
end

% Baseline composite residual (visible cells only) for the convergence test.
[res, H] = ebCompResidual(H);
hist = res.linf;
rel = res.rel;
if opts.Verbose >= 2
  fprintf('   it  0   |r|_inf = %.6e   rel = %.3e\n', res.linf, rel);
end

it = 0;
while it < opts.MaxVCycles && rel > opts.RelTol
  it = it + 1;
  % Restore the true right-hand side everywhere: the previous V-cycle's FAS
  % coarse loading replaced F on covered regions, so it must be reset before
  % the next descent.
  for k = 1:H.nlev
    for p = 1:H.lev{k}.np
      H.lev{k}.F{p} = H.lev{k}.FTrue{p};
    end
  end
  % One full V-cycle from the finest level, then push the fine solution back
  % down so covered coarse cells agree with their children before measuring.
  H = ebVcycle(H, H.nlev);
  H = ebTransfer(H, 0, 'filldown');
  [res, H] = ebCompResidual(H);
  hist(end+1) = res.linf; %#ok<AGROW>
  rel = res.rel;
  if opts.Verbose >= 2
    fprintf('   it %2d   |r|_inf = %.6e   rel = %.3e   factor = %.3f\n', ...
      it, res.linf, rel, hist(end) / max(hist(end-1), realmin));
  end
end

% Total cell count including covered cells (res.ncell counts visible only).
ncellTot = 0;
for k = H.rootIdx:H.nlev
  for p = 1:H.lev{k}.np
    ncellTot = ncellTot + numel(H.lev{k}.F{p});
  end
end

% Record this solve as one "stage" in the run statistics (one stage per
% regrid; resHist lets callers inspect the per-cycle contraction factors).
stage.amrLevels = H.nlev - H.rootIdx + 1;
stage.nlev = H.nlev;
stage.vcycles = it;
stage.resHist = hist;
stage.finalRel = rel;
stage.visibleCells = res.ncell;
stage.totalCells = ncellTot;
if numel(hist) >= 2
  stage.avgFactor = (hist(end) / max(hist(1), realmin))^(1 / (numel(hist)-1));
else
  stage.avgFactor = NaN;
end
H.stats.stages{end+1} = stage;

if opts.Verbose >= 1
  npmax = 0;
  for k = H.rootIdx:H.nlev, npmax = max(npmax, H.lev{k}.np); end
  fprintf(['[ebsam] stage %d: %d AMR level(s), %d V-cycle(s), ' ...
           '|r|/|f| = %.2e, avg factor %.3f, %d visible cells\n'], ...
    numel(H.stats.stages), stage.amrLevels, it, rel, stage.avgFactor, ...
    res.ncell);
end
end
