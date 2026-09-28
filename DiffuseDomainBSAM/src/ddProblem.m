function [prob, aux] = ddProblem(geom, epsi, box, bctype, data, varargin)
%ddProblem Diffuse-domain reformulation of an elliptic problem posed on an
%   irregular domain Omega_1 (described by geom) into the general solver's
%   form  -div(Ddd grad u) + Cdd u = fdd  on the rectangular box.
%
%   [prob, aux] = ddProblem(geom, eps, box, bctype, data, 'Name', value...)
%
%   geom   : from ddGeometry; phase field phi = (1 - tanh(3 psi/eps))/2
%            (phi ~ 1 inside Omega_1, ~ 0 outside; layer width O(eps)).
%   box    : [a b c d] computational rectangle containing Omega_1.
%   bctype : 'neumann'   (Feng et al. eq. 6.18; Li-Lowengrub-Ratz-Voigt)
%               beta*phi*u - div(phi*D grad u) = phi*f + |grad phi| g,
%            with g = D du/dn on Gamma (n outward of Omega_1);
%            'dirichlet' penalty formulation
%               -div(phi*D grad u) + phi*C u + (mu/eps^3)(1-phi)(u-g) = phi*f.
%   data   : struct with fields f, g (handles or scalars), and optionally
%            D (default 1), C (default: beta=1 for neumann, 0 dirichlet),
%            exact (sharp exact solution for error reporting).
%   Options: 'PhiFloor' (1e-8)  lower bound on phi in Ddd (and Cdd for the
%            Neumann case) keeping the operator uniformly elliptic;
%            'Mu' (1) penalty strength;  'OuterBC' ({'neumann',0}).
%
%   aux: phi, gradphimag, psi handles; tagfun(width) returning a geometric
%   tag function |psi| < width*eps for ebOptions('Indicator','fun',...);
%   inside(X,Y) mask of Omega_1.

validateattributes(epsi, {'numeric'}, {'real','scalar','finite','positive'}, mfilename, 'eps');
p = inputParser;
p.addParameter('PhiFloor', 1e-8, @(v) isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>=0 && v<1);
p.addParameter('Mu', 1.0, @(v) isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v) && v>0);
p.addParameter('OuterBC', {'neumann', 0});
p.parse(varargin{:});
o = p.Results;

Dfun = tof(getdef(data, 'D', 1));
ffun = tof(getdef(data, 'f', 0));
gfun = tof(getdef(data, 'g', 0));

psi = geom.psi;
phi = @(X,Y) 0.5 * (1 - tanh(3 * psi(X,Y) / epsi));
gphimag = @(X,Y) (3/(2*epsi)) * sech(3 * psi(X,Y) / epsi).^2 .* geom.gpsimag(X,Y);
phifl = @(X,Y) max(phi(X,Y), o.PhiFloor);

switch lower(bctype)
  case 'neumann'
    Cfun = tof(getdef(data, 'C', 1));
    Ddd = @(X,Y) phifl(X,Y) .* Dfun(X,Y);
    Cdd = @(X,Y) phifl(X,Y) .* Cfun(X,Y);
    fdd = @(X,Y) phi(X,Y) .* ffun(X,Y) + gphimag(X,Y) .* gfun(X,Y);
  case 'dirichlet'
    Cfun = tof(getdef(data, 'C', 0));
    pen = o.Mu / epsi^3;
    Ddd = @(X,Y) phifl(X,Y) .* Dfun(X,Y);
    Cdd = @(X,Y) phi(X,Y) .* Cfun(X,Y) + pen * (1 - phi(X,Y));
    fdd = @(X,Y) phi(X,Y) .* ffun(X,Y) + pen * (1 - phi(X,Y)) .* gfun(X,Y);
  otherwise
    error('dd_problem:bctype', 'bctype must be neumann or dirichlet');
end

exact = getdef(data, 'exact', []);
prob = ebProblem(box, Ddd, Cdd, fdd, o.OuterBC, exact);

aux.phi = phi;
aux.gradphimag = gphimag;
aux.psi = psi;
aux.epsi = epsi;
aux.inside = @(X,Y) psi(X,Y) < 0;
aux.tagfun = @(width) @(X,Y) abs(psi(X,Y)) < width * epsi;
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function fh = tof(v)
if isa(v, 'function_handle'), fh = v;
else, vv = double(v); fh = @(X,Y) vv * ones(size(X)); end
end
