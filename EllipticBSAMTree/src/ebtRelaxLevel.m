function H = ebtRelaxLevel(H, k, nsweeps)
%ebtRelaxLevel nsweeps red-black Gauss-Seidel sweeps on level k.
%   Colors use GLOBAL cell parity so that simultaneous patch updates of one
%   color are a true Gauss-Seidel half-sweep of the level union.  Ghosts
%   (coarse-fine, same-level, physical) are refreshed after each color.
%
%   Physical boundaries are relaxed EXACTLY: before each colored pass the
%   physical ghost cells are loaded with the partial ghost values
%   (quadratic Dirichlet: (8/3)g + (1/3)u_2; Neumann: h g) and the
%   u_1-coupling of the ghost is folded into the diagonal (InvDiagS, see
%   ebtBcFold), so boundary cells perform exact Gauss-Seidel updates of
%   the ghost-eliminated equations.  ebtFillGhosts then restores the true
%   ghost values for operator/residual evaluation.

% Hoist the smoother handle, the relax factor and the level's immutable
% per-patch coefficient/boundary arrays into locals once (rather than
% dereferencing H.K.relax and H.lev{k}.* per patch per colour per sweep).
% The partial physical ghosts are now set inside the kernel, so each patch
% update is a single in-place relax with no separate ghost-copy.
relax = H.K.relax;
om = H.opts.Omega;
hx = H.lev{k}.h(1); hy = H.lev{k}.h(2);
np  = H.lev{k}.np;
F   = H.lev{k}.F;          % force is fixed across the sweeps of this call
DeW = H.lev{k}.DeW;  DnS = H.lev{k}.DnS;  ID = H.lev{k}.InvDiagS;
bc  = H.lev{k}.bc;   par0 = H.lev{k}.par0;

for s = 1:nsweeps
  for color = 1:2
    Q = H.lev{k}.Q;        % COW alias; copied once on the first {p} write
    for p = 1:np
      Q{p} = relax(Q{p}, F{p}, DeW{p}, DnS{p}, ID{p}, bc{p}, par0(p), ...
                   color, hx, hy, om);
    end
    H.lev{k}.Q = Q;        % write the whole cell array back once per colour
    H = ebtFillGhosts(H, k);
  end
end
end
