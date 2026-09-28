function ind = ebtRte(H, k, p)
%ebtRte Relative truncation error indicator (paper sec. 4.4, the RTE
%   test of [45,47]) for patch p of level k, returned per FINE cell:
%
%       tau = L_{2h}( R u_h ) - R( L_h u_h ),
%
%   i.e. apply the level's operator coarsened by 2 (coefficients
%   restricted by face/cell averaging) to the restricted solution and
%   subtract the restricted fine operator; |tau| estimates the local
%   truncation error and is upsampled to the fine cells for tagging.
%   The 2h ghost ring is obtained by averaging the (already filled) fine
%   ghost pairs, consistent to the order needed for an indicator.

lev = H.lev{k};
Q = lev.Q{p};
mx = size(Q, 1) - 2; my = size(Q, 2) - 2;
cmx = mx/2; cmy = my/2;
hx = lev.h(1); hy = lev.h(2);

% R(L_h u)
Lh = H.K.op(Q, lev.DeW{p}, lev.DnS{p}, lev.CC{p}, hx, hy);
RLh = H.K.restrict2(Lh);

% R u with a ghost ring from averaged fine ghosts
Rpad = zeros(cmx+2, cmy+2);
Rpad(2:end-1, 2:end-1) = H.K.restrict2(Q(2:end-1, 2:end-1));
j = 1:cmy;
Rpad(1,     j+1) = 0.5*(Q(1,    2*j) + Q(1,    2*j+1));
Rpad(cmx+2, j+1) = 0.5*(Q(mx+2, 2*j) + Q(mx+2, 2*j+1));
i = 1:cmx;
Rpad(i+1, 1)     = 0.5*(Q(2*i, 1)    + Q(2*i+1, 1));
Rpad(i+1, cmy+2) = 0.5*(Q(2*i, my+2) + Q(2*i+1, my+2));

% coarsened coefficients (same averaging as the coefficient fill-down)
cDeW = 0.5*(lev.DeW{p}(1:2:mx+1, 1:2:my-1) + lev.DeW{p}(1:2:mx+1, 2:2:my));
cDnS = 0.5*(lev.DnS{p}(1:2:mx-1, 1:2:my+1) + lev.DnS{p}(2:2:mx, 1:2:my+1));
cCC  = H.K.restrict2(lev.CC{p});

% tau on the 2h grid, then copy each coarse value to its four fine cells.
tau = H.K.op(Rpad, cDeW, cDnS, cCC, 2*hx, 2*hy) - RLh;
ind = repelem(abs(tau), 2, 2);
end
