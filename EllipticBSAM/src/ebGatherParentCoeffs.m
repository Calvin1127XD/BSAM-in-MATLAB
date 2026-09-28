function H = ebGatherParentCoeffs(H)
%ebGatherParentCoeffs Build, per patch of every level k >= 2:
%   (i)  O(1) rect descriptors for gathering the parent coefficients over
%        the footprint (PCC cells, PDeW x-faces, PDnS y-faces) -- the
%        arrays themselves are materialized TRANSIENTLY by
%        ebParentCoeffs inside the coarse loading, not stored;
%   (ii) the persistent coarse interface-face diffusivities used by the
%        flux-balance mass corrections (perimeter-sized vectors):
%        cfD.W/E (cmy x 1) and cfD.S/N (cmx x 1).
%   Call after ebCoeffFilldown so values match the parent operator.

for k = 2:H.nlev
  levF = H.lev{k};
  levC = H.lev{k-1};
  for p = 1:levF.np
    b = levF.box(p, :);
    Il = (b(1)+1)/2; Ih = b(2)/2; Jl = (b(3)+1)/2; Jh = b(4)/2;
    cmx = Ih - Il + 1; cmy = Jh - Jl + 1;

    % One rect descriptor per parent patch that overlaps the footprint;
    % cell arrays and the two face orientations use slightly different
    % index ranges (faces extend one past the last cell).
    pcc  = struct('sp', {}, 'sr', {}, 'dr', {});
    pdew = struct('sp', {}, 'sr', {}, 'dr', {});
    pdns = struct('sp', {}, 'sr', {}, 'dr', {});
    for q = 1:levC.np
      qb = levC.box(q, :);
      r = rectIsect([Il Ih Jl Jh], qb);
      if ~isempty(r)
        pcc(end+1) = struct('sp', q, ...
          'sr', [r(1)-qb(1)+1, r(2)-qb(1)+1, r(3)-qb(3)+1, r(4)-qb(3)+1], ...
          'dr', [r(1)-Il+1,    r(2)-Il+1,    r(3)-Jl+1,    r(4)-Jl+1]); %#ok<*AGROW>
      end
      rf = rectIsect([Il, Ih+1, Jl, Jh], [qb(1), qb(2)+1, qb(3), qb(4)]);
      if ~isempty(rf)
        pdew(end+1) = struct('sp', q, ...
          'sr', [rf(1)-qb(1)+1, rf(2)-qb(1)+1, rf(3)-qb(3)+1, rf(4)-qb(3)+1], ...
          'dr', [rf(1)-Il+1,    rf(2)-Il+1,    rf(3)-Jl+1,    rf(4)-Jl+1]);
      end
      rf = rectIsect([Il, Ih, Jl, Jh+1], [qb(1), qb(2), qb(3), qb(4)+1]);
      if ~isempty(rf)
        pdns(end+1) = struct('sp', q, ...
          'sr', [rf(1)-qb(1)+1, rf(2)-qb(1)+1, rf(3)-qb(3)+1, rf(4)-qb(3)+1], ...
          'dr', [rf(1)-Il+1,    rf(2)-Il+1,    rf(3)-Jl+1,    rf(4)-Jl+1]);
      end
    end
    H.lev{k}.ops{p}.pcc = pcc;
    H.lev{k}.ops{p}.pdew = pdew;
    H.lev{k}.ops{p}.pdns = pdns;

    % materialize once to extract the interface-face vectors
    [PDeW, PDnS, ~] = ebParentCoeffs(H, k, p);
    if any(isnan(PDeW(:))) || any(isnan(PDnS(:)))
      error('eb_gather_parent_coeffs:cover', ...
        'footprint of patch %d level %d not covered by parent level.', p, k);
    end
    cfD.W = PDeW(1, :).';
    cfD.E = PDeW(cmx+1, :).';
    cfD.S = PDnS(:, 1);
    cfD.N = PDnS(:, cmy+1);
    H.lev{k}.cfD{p} = cfD;
  end
end
end

function r = rectIsect(a, b)
r = [max(a(1), b(1)), min(a(2), b(2)), max(a(3), b(3)), min(a(4), b(4))];
if r(1) > r(2) || r(3) > r(4), r = []; end
end
