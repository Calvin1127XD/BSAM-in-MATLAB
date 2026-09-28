function H = ebCompatProject(H)
%ebCompatProject Enforce the DISCRETE solvability (compatibility) of the
%   singular pure-Neumann problem on the composite mesh.
%
%   For all-Neumann sides with C == 0 the composite operator has the
%   constants in its kernel, so the discrete system is solvable iff the
%   conservation identity holds exactly for the data:
%
%       -sum_visible h_k^2 F  =  sum_bdry D g h_t        (*)
%
%   (the C == 0 case of the ebMassCheck identity; the right side is the
%   ghost-closure boundary flux, which for Neumann faces equals D*g
%   independently of u).  With midpoint source quadrature the left side
%   matches the true integral only to O(h^2) PER LEVEL, and on a
%   multi-level mesh the interlevel mismatch cannot be removed by the
%   solver: the coarsest KKT solve projects out only the coarse
%   representation of the incompatibility, and the V-cycle then stalls
%   at the O(h^2) floor (factor -> 1) -- the failure mode documented in
%   STRESS_REPORT.md section 4.
%
%   The fix is the standard one for singular Neumann problems: project
%   the source onto the compatible subspace by subtracting the constant
%
%       c = (sum_visible h^2 F + sum_bdry D g h_t) / area_visible
%
%   from FTrue on every level, which makes (*) hold to roundoff for ANY
%   quadrature (including numeric forcing data).  For a continuously
%   compatible problem c = O(h^2), so the computed solution is a
%   second-order perturbation and convergence rates are unaffected.
%   Called automatically after every hierarchy change when the problem
%   is singular (H.coarsest.singular); a no-op otherwise.

if ~isfield(H, 'coarsest')
  return;
end
if ~H.coarsest.singular
  % A newly sampled positive reaction can remove a previous nullspace.
  % Restore the unprojected source in that case.
  for k = 1:H.nlev
    if ~isfield(H.lev{k},'sourceShift'), continue; end
    old = H.lev{k}.sourceShift;
    if old == 0, continue; end
    for p = 1:H.lev{k}.np
      H.lev{k}.FTrue{p} = H.lev{k}.FTrue{p} + old;
      H.lev{k}.F{p} = H.lev{k}.FTrue{p};
    end
    H.lev{k}.sourceShift = 0;
  end
  H.compatibility = struct('sourceShift',0,'unprojectedDefect',0);
  return;
end

% ---- defect ingredients over the VISIBLE composite mesh (mirrors
% ebMassCheck's masks and ghost-closure flux orientation)
S = 0; Av = 0; R = 0;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  old = 0;
  if isfield(lev,'sourceShift'), old = lev.sourceShift; end
  hx = lev.h(1); hy = lev.h(2);
  da = hx * hy;
  for p = 1:lev.np
    cov = lev.covered{p};
    vis = ~cov;
    % Integrate the ORIGINAL source, restoring the previous shift in the
    % sum.  A newly created level has sourceShift=0; old levels may have
    % already been projected.  Mixing those representations would add a
    % different source perturbation on each level after regridding.
    S = S + (sum(lev.FTrue{p}(vis)) + old*nnz(vis)) * da;
    Av = Av + nnz(vis) * da;
    [mx, my] = size(cov);
    bcd = lev.bc{p};
    if bcd(1).code == 2                       % west Neumann: flux = D*g
      m = ~cov(1, :);
      fl = lev.DeW{p}(1, :) .* bcd(1).g;
      R = R + sum(fl(m)) * hy;
    end
    if bcd(2).code == 2
      m = ~cov(mx, :);
      fl = lev.DeW{p}(mx+1, :) .* bcd(2).g;
      R = R + sum(fl(m)) * hy;
    end
    if bcd(3).code == 2
      m = ~cov(:, 1);
      fl = lev.DnS{p}(:, 1) .* bcd(3).g;
      R = R + sum(fl(m)) * hx;
    end
    if bcd(4).code == 2
      m = ~cov(:, my);
      fl = lev.DnS{p}(:, my+1) .* bcd(4).g;
      R = R + sum(fl(m)) * hx;
    end
  end
end

c = (S + R) / max(Av, realmin);

% ---- subtract the constant from the source on EVERY level (so
% restriction consistency R(F - c) = RF - c holds and the FAS coarse
% loads stay compatible); F re-aliases the projected FTrue
for k = 1:H.nlev
  old = 0;
  if isfield(H.lev{k},'sourceShift'), old = H.lev{k}.sourceShift; end
  for p = 1:H.lev{k}.np
    H.lev{k}.FTrue{p} = H.lev{k}.FTrue{p} - (c - old);
    H.lev{k}.F{p} = H.lev{k}.FTrue{p};
  end
  H.lev{k}.sourceShift = c;
end
H.compatibility = struct('sourceShift',c,'unprojectedDefect',S+R);
end
