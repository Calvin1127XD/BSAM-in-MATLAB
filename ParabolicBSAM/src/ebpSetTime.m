function H = ebpSetTime(H, pprob, t, aDt, force)
%EBPSETTIME Refresh the hierarchy's time-frozen state to time t and
%   implicit diagonal shift aDt = alpha/Dt:
%     CC  <- C(x,y,t) + aDt     (the implicit operator's reaction)
%     DeW/DnS <- D(x,y,t)       (only when D is time-dependent)
%     bc data <- g(x,y,t)       (only when the BC data are)
%   then coefficient fill-down, parent-coefficient/interface gathers,
%   relaxation diagonals and the coarsest factorization.  Cheap path: if
%   nothing is time-dependent and aDt is unchanged, this is a no-op.

if nargin < 5, force = false; end
if ~isfield(H, 'par'), H.par = struct(); end
key = [t * (pprob.tdD || pprob.tdC || pprob.tdBC), aDt];
if ~force && isfield(H.par, 'stateKey') && isequal(H.par.stateKey, key)
  return;
end
dom = pprob.domain;

for k = 1:H.nlev
  lev = H.lev{k};
  for p = 1:lev.np
    b = lev.box(p, :);
    hx = lev.h(1); hy = lev.h(2);
    xc = dom(1) + ((b(1):b(2)) - 0.5) * hx;
    yc = dom(3) + ((b(3):b(4)) - 0.5) * hy;
    if pprob.tdD || force
      xf = dom(1) + (b(1)-1:b(2)) * hx;
      yf = dom(3) + (b(3)-1:b(4)) * hy;
      [XF, YC] = ndgrid(xf, yc);
      H.lev{k}.DeW{p} = expand(pprob.D(XF, YC, t), XF);
      [XC, YF] = ndgrid(xc, yf);
      H.lev{k}.DnS{p} = expand(pprob.D(XC, YF, t), XC);
    end
    [X, Y] = ndgrid(xc, yc);
    H.lev{k}.CC{p} = expand(pprob.C(X, Y, t), X) + aDt;
    if pprob.tdBC || force
      bcd = H.lev{k}.bc{p};
      fint = lev.faceInternal{p};
      if ~fint(1), bcd(1).g = reshape(pprob.bc(1).g(dom(1)*ones(size(yc)), yc, t), 1, []); end
      if ~fint(2), bcd(2).g = reshape(pprob.bc(2).g(dom(2)*ones(size(yc)), yc, t), 1, []); end
      if ~fint(3), bcd(3).g = reshape(pprob.bc(3).g(xc, dom(3)*ones(size(xc)), t), [], 1); end
      if ~fint(4), bcd(4).g = reshape(pprob.bc(4).g(xc, dom(4)*ones(size(xc)), t), [], 1); end
      H.lev{k}.bc{p} = bcd;
    end
  end
end

% restriction of coefficients + InvDiagS refresh (does all levels)
H = ebCoeffFilldown(H);
H = ebGatherParentCoeffs(H);
H = ebCoarsest(H, 'assemble');

H.par.stateKey = key;
H.par.aDt = aDt;
end

function v = expand(v, X)
if isscalar(v), v = v * ones(size(X)); end
end
