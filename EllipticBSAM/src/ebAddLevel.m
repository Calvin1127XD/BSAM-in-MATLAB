function H = ebAddLevel(H, boxes)
%ebAddLevel Append a new finest level whose patches are the 2x refinement
%   of the given parent-space boxes, rebuild the static machinery, and
%   initialize the new level by bilinear prolongation of the parent data.

kNew = H.nlev + 1;
nNew = H.lev{H.nlev}.n * 2;
% Convert parent-cell boxes to fine indices: coarse cell c covers fine
% cells 2c-1 and 2c, so [lo hi] -> [2*lo-1, 2*hi] in each direction.
fineBoxes = [2*boxes(:,1) - 1, 2*boxes(:,2), 2*boxes(:,3) - 1, 2*boxes(:,4)];

H.nlev = kNew;
H.lev{kNew} = ebMakeLevel(H, kNew, nNew, fineBoxes);

% Rebuild everything that depends on the mesh: averaged coefficients on
% the new level, the transfer slice tables, the parent-coefficient gathers,
% and the prefactored coarsest system.
H = ebCoeffFilldown(H);
H = ebBuildOps(H, kNew);
H = ebGatherParentCoeffs(H);
H = ebCoarsest(H, 'assemble');     % coarse coefficients changed
H = ebCompatProject(H);            % re-project the source for the new
                                   % visible mesh (singular Neumann only)

% initial guess: bilinear prolongation of the parent window
H = ebFillGhosts(H, kNew - 1);
lev = H.lev{kNew};
for p = 1:lev.np
  W = ebGatherWindow(H, kNew, p);
  H.lev{kNew}.Q{p}(2:end-1, 2:end-1) = H.K.prolong(W);
end
H = ebFillGhosts(H, kNew);

if H.opts.Verbose >= 1
  ncell = 0;
  for p = 1:lev.np, ncell = ncell + numel(H.lev{kNew}.F{p}); end
  fprintf('[ebsam] added level %d (%d patch(es), %d cells, h = %.3e)\n', ...
    kNew - H.rootIdx + 1, lev.np, ncell, lev.h(1));
end
end
