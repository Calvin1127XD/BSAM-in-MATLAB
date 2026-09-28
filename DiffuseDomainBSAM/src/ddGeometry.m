function geom = ddGeometry(kind, varargin)
%ddGeometry Level-set descriptions of internal domains Omega_1.
%   geom = ddGeometry('circle', R0 [, cx, cy])
%   geom = ddGeometry('ellipse', a, b [, cx, cy])
%   geom = ddGeometry('star', b)        flower of Feng et al. eq. (6.19):
%          R0(theta) = (1/10) e^{2b} cos(3 theta) + (1/50) e^{6b} cos(5 theta) + 1
%   geom = ddGeometry('polygon', Nv, R [, rot])   regular Nv-gon with
%          circumradius R (EXACT signed distance; |grad psi| = 1 a.e.,
%          gradient jumps along the corner bisectors)
%   geom = ddGeometry('pentagon', R [, rot])      = polygon with Nv = 5
%   geom = ddGeometry('polar', Rfun, Rpfun, name) generic simple closed
%          polar curve r = R(theta) > 0 (Rpfun = dR/dtheta)
%   geom = ddGeometry('flower3')   r = 1 + 0.25 cos(3 theta)
%   geom = ddGeometry('peanut')    r = 1 + 0.40 cos(2 theta)
%
%   Returns handles:
%     geom.psi(X,Y)     level set, psi < 0 inside Omega_1, psi = 0 on Gamma
%     geom.gpsi(X,Y)    -> {psix, psiy}
%     geom.gpsimag(X,Y) |grad psi|
%   The phase field is phi = (1 - tanh(3 psi / eps)) / 2.

switch lower(kind)
  case 'circle'
    R0 = varargin{1};
    if numel(varargin) >= 3, cx = varargin{2}; cy = varargin{3};
    else, cx = 0; cy = 0; end
    geom.psi = @(X,Y) hypot(X-cx, Y-cy) - R0;
    geom.gpsi = @(X,Y) circGrad(X-cx, Y-cy);
    geom.gpsimag = @(X,Y) ones(size(X));
    geom.name = sprintf('circle R=%g', R0);

  case 'ellipse'
    a = varargin{1}; b = varargin{2};
    if numel(varargin) >= 4, cx = varargin{3}; cy = varargin{4};
    else, cx = 0; cy = 0; end
    % psi = sqrt((x/a)^2 + (y/b)^2) - 1 (approximate distance scaling)
    geom.psi = @(X,Y) hypot((X-cx)/a, (Y-cy)/b) - 1;
    geom.gpsi = @(X,Y) ellGrad(X-cx, Y-cy, a, b);
    geom.gpsimag = @(X,Y) ellGmag(X-cx, Y-cy, a, b);
    geom.name = sprintf('ellipse a=%g b=%g', a, b);

  case 'star'
    if isempty(varargin), bpar = 0.5; else, bpar = varargin{1}; end
    A3 = exp(2*bpar)/10; A5 = exp(6*bpar)/50;
    Rf  = @(th) A3*cos(3*th) + A5*cos(5*th) + 1;
    Rpf = @(th) -3*A3*sin(3*th) - 5*A5*sin(5*th);
    geom = polarGeom(Rf, Rpf, sprintf('star b=%g', bpar));

  case 'polar'
    geom = polarGeom(varargin{1}, varargin{2}, varargin{3});

  case 'flower3'
    geom = polarGeom(@(th) 1 + 0.25*cos(3*th), @(th) -0.75*sin(3*th), ...
                      'flower3 r=1+0.25cos3t');

  case 'peanut'
    geom = polarGeom(@(th) 1 + 0.40*cos(2*th), @(th) -0.80*sin(2*th), ...
                      'peanut r=1+0.4cos2t');

  case {'polygon', 'pentagon'}
    if strcmpi(kind, 'pentagon')
      Nv = 5; R = varargin{1};
      if numel(varargin) >= 2, rot = varargin{2}; else, rot = pi/2; end
    else
      Nv = varargin{1}; R = varargin{2};
      if numel(varargin) >= 3, rot = varargin{3}; else, rot = pi/2; end
    end
    th = rot + 2*pi*(0:Nv-1)/Nv;     % CCW vertex ordering
    vx = R*cos(th); vy = R*sin(th);
    geom.psi = @(X,Y) polySdf(X, Y, vx, vy);
    geom.gpsi = @(X,Y) polyGrad(X, Y, vx, vy);
    geom.gpsimag = @(X,Y) ones(size(X));   % exact distance: |grad psi|=1 a.e.
    geom.name = sprintf('regular %d-gon R=%g', Nv, R);
    geom.vertices = [vx(:), vy(:)];

  otherwise
    error('dd_geometry:kind', 'unknown geometry %s', kind);
end
end

% ------------------------------------------------------------------------
function geom = polarGeom(Rf, Rpf, name)
%polarGeom Simple closed polar curve r = R(theta) > 0:
%   psi = r - R(theta);  grad theta = (-y, x)/r^2; grad r = (x, y)/r.
geom.psi = @(X,Y) polarPsi(X, Y, Rf);
geom.gpsi = @(X,Y) polarGrad(X, Y, Rpf);
geom.gpsimag = @(X,Y) polarGmag(X, Y, Rpf);
geom.name = name;
end

function v = polarPsi(X, Y, Rf)
r = hypot(X, Y);
th = atan2(Y, X);
v = r - Rf(th);
end

function g = polarGrad(X, Y, Rpf)
r = max(hypot(X, Y), 1e-10);
th = atan2(Y, X);
Rp = Rpf(th);
g = {X./r + Rp.*Y./r.^2, Y./r - Rp.*X./r.^2};
end

function m = polarGmag(X, Y, Rpf)
g = polarGrad(X, Y, Rpf);
m = hypot(g{1}, g{2});
end

% ------------------------------------------------------------------------
function g = circGrad(X, Y)
r = max(hypot(X, Y), 1e-12);
g = {X./r, Y./r};
end

function g = ellGrad(X, Y, a, b)
r = max(hypot(X/a, Y/b), 1e-12);
g = {X./(a^2*r), Y./(b^2*r)};
end

function m = ellGmag(X, Y, a, b)
g = ellGrad(X, Y, a, b);
m = hypot(g{1}, g{2});
end

% ------------------------------------------------------------------------
function [psi, gx, gy] = polyCore(X, Y, vx, vy)
%polyCore Exact signed distance to a CONVEX polygon (CCW vertices) and
%   its gradient (outward unit direction to/from the closest boundary
%   feature; jumps along corner bisectors, |grad| = 1 a.e.).
sz = size(X);
px = X(:); py = Y(:);
n = numel(vx);
dmin = inf(numel(px), 1);
dvx = zeros(numel(px), 1); dvy = zeros(numel(px), 1);
edgeNx = zeros(numel(px), 1); edgeNy = zeros(numel(px), 1);
inside = true(numel(px), 1);
for k = 1:n
  ax = vx(k); ay = vy(k);
  kb = mod(k, n) + 1;
  bx = vx(kb); by = vy(kb);
  ex = bx - ax; ey = by - ay;
  t = ((px - ax)*ex + (py - ay)*ey) / (ex^2 + ey^2);
  t = min(max(t, 0), 1);
  dx = px - (ax + t*ex);
  dy = py - (ay + t*ey);
  d = hypot(dx, dy);
  m = d < dmin;
  dmin(m) = d(m); dvx(m) = dx(m); dvy(m) = dy(m);
  edgeNx(m) = ey / hypot(ex, ey); edgeNy(m) = -ex / hypot(ex, ey);
  inside = inside & (ex*(py - ay) - ey*(px - ax) >= 0);
end
s = ones(numel(px), 1);
s(inside) = -1;
psi = reshape(s .* dmin, sz);
dn = max(dmin, realmin);
gx = reshape(s .* dvx ./ dn, sz);
gy = reshape(s .* dvy ./ dn, sz);
% On a smooth polygon edge the closest-point displacement vanishes,
% but the outward normal is well defined. Dividing by a distance floor
% incorrectly returned zero there (and attenuated nearby normals).
onEdge = dmin <= 64*eps(max(1,max(abs([vx(:);vy(:)]))));
gx(onEdge) = edgeNx(onEdge); gy(onEdge) = edgeNy(onEdge);
end

function psi = polySdf(X, Y, vx, vy)
[psi, ~, ~] = polyCore(X, Y, vx, vy);
end

function g = polyGrad(X, Y, vx, vy)
[~, gx, gy] = polyCore(X, Y, vx, vy);
g = {gx, gy};
end
