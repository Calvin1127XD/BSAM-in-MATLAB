function H = ebtAddLevel(H, boxes)
%ebtAddLevel Append a new finest level whose patches are the 2x refinement
%   of the given parent-space cluster boxes, rebuild the static machinery,
%   and initialize the new level by bilinear prolongation of the parent
%   window.
%
%   boxes : nb x 5 rows [parentIdx, ilo ihi jlo jhi], where parentIdx is
%           the patch of the CURRENT finest level that produced the box
%           (per-parent clustering, ebtTagCluster) and [ilo..jhi] is the
%           box in that level's GLOBAL cell space, guaranteed to lie
%           inside the parent patch.  The child is its exact 2x
%           refinement, hence automatically odd/even-aligned at the new
%           level (Fortran MakeNewGrid: child mglobal = [2*g-1, 2*g]).

kNew = H.nlev + 1;
nNew = H.lev{H.nlev}.n * 2;
parents = boxes(:, 1);
cb = boxes(:, 2:5);
% Coarse cell c covers fine cells 2c-1 and 2c: [lo hi] -> [2*lo-1, 2*hi].
fineBoxes = [2*cb(:,1) - 1, 2*cb(:,2), 2*cb(:,3) - 1, 2*cb(:,4)];

H.nlev = kNew;
H.lev{kNew} = ebtMakeLevel(H, kNew, nNew, fineBoxes, parents);

% Rebuild everything that depends on the mesh: averaged coefficients on
% the new level, the single-parent slice tables, the parent-coefficient
% windows, and the prefactored coarsest system.
H = ebtCoeffFilldown(H);
H = ebtBuildOps(H, kNew);
H = ebtGatherParentCoeffs(H);
H = ebtCoarsest(H, 'assemble');     % coarse coefficients changed
H = ebtCompatProject(H);            % re-project the source for the new
                                    % visible mesh (singular Neumann only)

% initial guess: bilinear prolongation of the single-parent window
H = ebtFillGhosts(H, kNew - 1);
lev = H.lev{kNew};
for p = 1:lev.np
  W = ebtGatherWindow(H, kNew, p);
  H.lev{kNew}.Q{p}(2:end-1, 2:end-1) = H.K.prolong(W);
end
H = ebtFillGhosts(H, kNew);

if H.opts.Verbose >= 1
  ncell = 0;
  for p = 1:lev.np, ncell = ncell + numel(H.lev{kNew}.F{p}); end
  fprintf('[ebtsam] added level %d (%d patch(es), %d cells, h = %.3e)\n', ...
    kNew - H.rootIdx + 1, lev.np, ncell, lev.h(1));
end
end
