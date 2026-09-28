function H = ebpInit(pprob, opts, aDt1)
%EBPINIT Build the initial hierarchy at t = 0, adapted to the initial
%   condition: set u0 on the root, tag/refine progressively (each new
%   level gets u0 evaluated ANALYTICALLY, not interpolated), and store
%   the history field Uo = u0.

% Elliptic hierarchy for the very first implicit solve (its coefficients
% embed the time discretization via aDt1).
probE = ebpEprob(pprob, 0, aDt1);
H = ebSetup(probE, opts.eb);
H = setU0(H, pprob);
H = ebFillGhostsAll(H);

% Grow the hierarchy against the initial condition: tag on u0, refine,
% re-evaluate u0 analytically on the new level, repeat.
while (H.nlev - H.rootIdx + 1) < opts.eb.MaxLevels
  boxes = ebTagCluster(H);
  if isempty(boxes), break; end
  H = ebAddLevel(H, boxes);
  H = setU0(H, pprob, H.nlev);
  H = ebFillGhosts(H, H.nlev);
end
H = ebTransfer(H, 0, 'filldown');

% Time-stepping history: Uo = u^n (here u0), Uoo = u^{n-1} (BDF2 only).
H.par = struct();
H.par.Uo = interiors(H);
H.par.Uoo = [];
H.par.stateKey = [0 * (pprob.tdD || pprob.tdC || pprob.tdBC), aDt1];
H.par.aDt = aDt1;
end

% ------------------------------------------------------------------------
function H = setU0(H, pprob, konly)
dom = pprob.domain;
if nargin < 3, ks = 1:H.nlev; else, ks = konly; end
for k = ks
  lev = H.lev{k};
  for p = 1:lev.np
    b = lev.box(p, :);
    [X, Y] = ndgrid(dom(1) + ((b(1):b(2)) - 0.5)*lev.h(1), ...
                    dom(3) + ((b(3):b(4)) - 0.5)*lev.h(2));
    H.lev{k}.Q{p}(2:end-1, 2:end-1) = pprob.u0(X, Y);
  end
end
end

function H = ebFillGhostsAll(H)
for k = 1:H.nlev
  H = ebFillGhosts(H, k);
end
end

function U = interiors(H)
U = cell(1, H.nlev);
for k = 1:H.nlev
  for p = 1:H.lev{k}.np
    U{k}{p} = H.lev{k}.Q{p}(2:end-1, 2:end-1);
  end
end
end
