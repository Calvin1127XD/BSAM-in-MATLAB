function corr = ebtFaceCorrections(H, k, p)
%ebtFaceCorrections Flux-balance mass corrections for patch p of level k
%   (Feng et al. 2018, eq. 4.20).  For each coarse face on the coarse-fine
%   interface the correction equals
%
%     C = [ sum of the two fine fluxes  D_f (u_in - u_ghost)
%           + coarse flux              D_c (u_out - [R u]_in) ] / h_c^2 ,
%
%   which is exactly  M (Z - [R u]) / h_c^2  with Z the flux-balanced
%   zombie value (eq. 4.10).  Added to the coarse force at the cell just
%   outside the footprint, it makes the coarse equation see the fine
%   fluxes, preserving the composite mass (Theorem 4.2).
%
%   Requires: fine ghosts current, QC ring gathered (u_out = QC ring,
%   [R u]_in = QC first interior row/col).
%   Returns corr{1..4} = W,E,S,N structs with .ops (scatter) and .vals.

lev = H.lev{k};
Q  = lev.Q{p};
QC = lev.QC{p};
DeW = lev.DeW{p}; DnS = lev.DnS{p};
cfD = lev.cfD{p};      % coarse interface-face diffusivities (vectors)
mx = size(Q, 1) - 2; my = size(Q, 2) - 2;
cmx = mx/2; cmy = my/2;
hcx2 = H.lev{k-1}.h(1)^2;
hcy2 = H.lev{k-1}.h(2)^2;
ops = lev.ops{p};

% WEST: fine ghost col 1, interior col 2; coarse out = QC(1,:), in = QC(2,:)
% Coarse face j spans fine rows 2j and 2j+1 in the padded array (interior
% rows start at 2), with fine face diffusivities DeW(:,2j-1) and DeW(:,2j).
j = (1:cmy).';
vals = ( DeW(1, 2*j-1).' .* (Q(2, 2*j).'   - Q(1, 2*j).') ...
       + DeW(1, 2*j).'   .* (Q(2, 2*j+1).' - Q(1, 2*j+1).') ...
       + cfD.W .* (QC(1, j+1).'  - QC(2, j+1).') ) / hcx2;
corr{1} = pack(ops.corrW, vals);

% EAST
vals = ( DeW(mx+1, 2*j-1).' .* (Q(mx+1, 2*j).'   - Q(mx+2, 2*j).') ...
       + DeW(mx+1, 2*j).'   .* (Q(mx+1, 2*j+1).' - Q(mx+2, 2*j+1).') ...
       + cfD.E .* (QC(cmx+2, j+1).' - QC(cmx+1, j+1).') ) / hcx2;
corr{2} = pack(ops.corrE, vals);

% SOUTH
i = (1:cmx).';
vals = ( DnS(2*i-1, 1) .* (Q(2*i, 2)   - Q(2*i, 1)) ...
       + DnS(2*i, 1)   .* (Q(2*i+1, 2) - Q(2*i+1, 1)) ...
       + cfD.S .* (QC(i+1, 1)  - QC(i+1, 2)) ) / hcy2;
corr{3} = pack(ops.corrS, vals);

% NORTH
vals = ( DnS(2*i-1, my+1) .* (Q(2*i, my+1)   - Q(2*i, my+2)) ...
       + DnS(2*i, my+1)   .* (Q(2*i+1, my+1) - Q(2*i+1, my+2)) ...
       + cfD.N .* (QC(i+1, cmy+2) - QC(i+1, cmy+1)) ) / hcy2;
corr{4} = pack(ops.corrN, vals);
end

function c = pack(cs, vals)
c.ops = cs.ops;
c.vals = vals;
end
