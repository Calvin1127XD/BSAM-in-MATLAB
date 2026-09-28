function prob = ebProblem(domain, D, C, f, bc, exact)
%ebProblem Define  -div(D(x,y) grad u) + C(x,y) u = f(x,y)  on [a,b]x[c,d].
%
%   prob = ebProblem(domain, D, C, f, bc [, exact])
%
%   domain : [a b c d]
%   D,C,f  : function handles @(X,Y) -> array (vectorized) or scalars.
%            D > 0 required; C >= 0 recommended (M-matrix / solvability).
%   bc     : boundary conditions. Either a single cell {type, g} applied
%            to all four sides, or a struct with fields west/east/south/
%            north, each {type, g}.  type is 'dirichlet' or 'neumann';
%            g is @(x,y) or a scalar.  Neumann data is the OUTWARD normal
%            derivative du/dn on that side (so the boundary diffusive
%            flux is -D*g).
%   exact  : optional @(X,Y) exact solution (for error reporting).
%
%   The solver discretizes with cell-centered finite differences in flux
%   form (paper eq. 4.4-4.5), with D evaluated at face midpoints.

if numel(domain) ~= 4 || domain(2) <= domain(1) || domain(4) <= domain(3)
  error('eb_problem:domain', 'domain must be [a b c d] with b>a, d>c.');
end
prob.domain = domain(:).';
prob.D = tofun(D, 'D');
prob.C = tofun(C, 'C');
prob.f = tofun(f, 'f');

sides = {'west', 'east', 'south', 'north'};
% Shorthand form {type, g}: replicate the one spec to all four sides.
if iscell(bc) && numel(bc) == 2
  one = bc; bc = struct();
  for s = 1:4, bc.(sides{s}) = one; end
end
for s = 1:4
  if ~isfield(bc, sides{s})
    error('eb_problem:bc', 'bc missing side %s', sides{s});
  end
  spec = bc.(sides{s});
  t = lower(spec{1});
  if ~ismember(t, {'dirichlet', 'neumann'})
    error('eb_problem:bc', 'bc type must be dirichlet or neumann');
  end
  prob.bc(s).type = t;
  prob.bc(s).g = tofun(spec{2}, ['g_' sides{s}]);
end

if nargin >= 6 && ~isempty(exact)
  prob.exact = tofun(exact, 'exact');
else
  prob.exact = [];
end
end

function fh = tofun(v, label)
if isa(v, 'function_handle')
  fh = v;
elseif isnumeric(v) && isscalar(v)
  vv = double(v);
  fh = @(X, Y) vv * ones(size(X));
else
  error('eb_problem:badinput', '%s must be a function handle or scalar.', label);
end
end
