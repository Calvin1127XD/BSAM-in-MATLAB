function probE = ebpEprob(pprob, t, aDt)
%EBPEPROB Freeze the parabolic problem at time t with implicit diagonal
%   shift aDt = alpha/Dt into an elliptic ebProblem:
%     -div(D(.,t) grad u) + (C(.,t) + alpha/Dt) u = <set per step>.
%   The source handle is a placeholder; ebpRhs overwrites FTrue each step.

Dh = @(x, y) pprob.D(x, y, t);
Ch = @(x, y) pprob.C(x, y, t) + aDt;
fh = @(x, y) zeros(size(x));

sides = {'west', 'east', 'south', 'north'};
bc = struct();
for s = 1:4
  g = pprob.bc(s).g;
  bc.(sides{s}) = {pprob.bc(s).type, @(x, y) g(x, y, t)};
end

probE = ebProblem(pprob.domain, Dh, Ch, fh, bc, []);
end
