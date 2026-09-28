function H = ebtGatherParentCoeffs(H)
%ebtGatherParentCoeffs Extract, per patch of every level k >= 2, the
%   persistent coarse interface-face diffusivities used by the
%   flux-balance mass corrections (perimeter-sized vectors):
%   cfD.W/E (cmy x 1) and cfD.S/N (cmx x 1).
%
%   TREE VERSION: the coefficient-window descriptors (ops.pcc/pdew/pdns)
%   are already single-parent slices built by ebtBuildOps, so this
%   routine only materializes them once (ebtParentCoeffs) and keeps the
%   boundary faces.  Call after ebtCoeffFilldown so the values match the
%   parent operator.

for k = 2:H.nlev
  levF = H.lev{k};
  for p = 1:levF.np
    mb = levF.mbounds(p, :);
    cmx = mb(2) - mb(1) + 1;
    cmy = mb(4) - mb(3) + 1;

    [PDeW, PDnS, ~] = ebtParentCoeffs(H, k, p);
    if ~isequal(size(PDeW), [cmx+1, cmy]) || ~isequal(size(PDnS), [cmx, cmy+1])
      error('ebt_gather_parent_coeffs:size', ...
        'parent coefficient window of patch %d level %d has wrong size.', p, k);
    end
    cfD.W = PDeW(1, :).';
    cfD.E = PDeW(cmx+1, :).';
    cfD.S = PDnS(:, 1);
    cfD.N = PDnS(:, cmy+1);
    H.lev{k}.cfD{p} = cfD;
  end
end
end
