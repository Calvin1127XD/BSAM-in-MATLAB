function pprob = ebpProblem(domain, D, C, f, bc, u0, exact)
%EBPPROBLEM Define the parabolic problem
%     du/dt - div(D grad u) + C u = f   on [a,b] x [c,d],
%   with initial condition u0 and Dirichlet/Neumann data per side.
%
%   pprob = ebpProblem(domain, D, C, f, bc, u0 [, exact])
%
%   D, C : @(X,Y) or @(X,Y,T) (or scalars); D > 0, C >= 0.
%   f    : @(X,Y,T) (or @(X,Y) autonomous, or scalar).
%   bc   : {type, g} for all sides or struct west/east/south/north with
%          {type, g}; g is @(x,y), @(x,y,t) or a scalar; Neumann data is
%          the OUTWARD normal derivative.
%   u0   : @(X,Y) initial condition.
%   exact: optional @(X,Y,T) exact solution (for error reporting).
%
%   Time-dependence of each ingredient is detected from its arity, so the
%   solver only re-evaluates what actually changes between steps.

if numel(domain) ~= 4 || domain(2) <= domain(1) || domain(4) <= domain(3)
  error('ebp_problem:domain', 'domain must be [a b c d].');
end
pprob.domain = domain(:).';

[pprob.D, pprob.tdD] = wrapT(D, 'D');
[pprob.C, pprob.tdC] = wrapT(C, 'C');
[pprob.f, pprob.tdF] = wrapT(f, 'f');

sides = {'west', 'east', 'south', 'north'};
if iscell(bc) && numel(bc) == 2
  one = bc; bc = struct();
  for s = 1:4, bc.(sides{s}) = one; end
end
pprob.tdBC = false;
for s = 1:4
  spec = bc.(sides{s});
  t = lower(spec{1});
  if ~ismember(t, {'dirichlet', 'neumann'})
    error('ebp_problem:bc', 'bc type must be dirichlet or neumann');
  end
  [gfun, td] = wrapT(spec{2}, ['g_' sides{s}]);
  pprob.bc(s).type = t;
  pprob.bc(s).g = gfun;
  pprob.tdBC = pprob.tdBC || td;
end

if ~isa(u0, 'function_handle')
  error('ebp_problem:u0', 'u0 must be @(X,Y).');
end
pprob.u0 = u0;

if nargin >= 7 && ~isempty(exact)
  pprob.exact = exact;          % @(X,Y,T)
else
  pprob.exact = [];
end
end

% ------------------------------------------------------------------------
function [fh, td] = wrapT(v, label)
%WRAPT Normalize to @(X,Y,T); td = true if genuinely time-dependent.
if isa(v, 'function_handle')
  na = nargin(v);
  if na == 3 || na < 0          % negative for varargs: assume t-dependent
    fh = v;
    td = true;
  elseif na == 2
    fh = @(X, Y, T) v(X, Y);
    td = false;
  else
    error('ebp_problem:arity', '%s must take (x,y) or (x,y,t).', label);
  end
elseif isnumeric(v) && isscalar(v)
  vv = double(v);
  fh = @(X, Y, T) vv * ones(size(X));
  td = false;
else
  error('ebp_problem:badinput', '%s must be a handle or scalar.', label);
end
end
