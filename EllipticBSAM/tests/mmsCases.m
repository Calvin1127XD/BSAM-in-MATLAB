function cs = mmsCases(name)
%mmsCases Manufactured test problems for the suite.
switch name
  case 'vardc'
    % variable D and C on [0,1]^2, Dirichlet
    uex = @(x,y) sin(pi*x).*cos(pi*y);
    D = @(x,y) 1 + 0.5*sin(2*pi*x).*cos(pi*y);
    C = @(x,y) 1 + x.*y;
    Dx = @(x,y) pi*cos(2*pi*x).*cos(pi*y);
    Dy = @(x,y) -0.5*pi*sin(2*pi*x).*sin(pi*y);
    ux = @(x,y) pi*cos(pi*x).*cos(pi*y);
    uy = @(x,y) -pi*sin(pi*x).*sin(pi*y);
    lap = @(x,y) -2*pi^2*uex(x,y);
    f = @(x,y) -(Dx(x,y).*ux(x,y) + Dy(x,y).*uy(x,y) + D(x,y).*lap(x,y)) ...
             + C(x,y).*uex(x,y);
    cs.prob = ebProblem([0 1 0 1], D, C, f, {'dirichlet', uex}, uex);

  case 'pureneumann'
    % singular pure-Neumann Poisson on [0,1]^2 (compatible f, zero flux)
    uex = @(x,y) cos(2*pi*x).*cos(2*pi*y);
    f = @(x,y) 8*pi^2*uex(x,y);
    cs.prob = ebProblem([0 1 0 1], 1, 0, f, {'neumann', 0}, uex);

  case 'gauss51'
    % periodic-symmetric Gaussian, paper eq. (5.1): u - lap(u) = f,
    % exact homogeneous Neumann on (-0.25,0.25)^2
    sg = 0.25;
    gp  = @(t) -pi*sin(4*pi*t)/sg^2;
    gpp = @(t) -4*pi^2*cos(4*pi*t)/sg^2;
    uex = @(x,y) exp(-(sin(2*pi*x).^2 + sin(2*pi*y).^2)/(2*sg^2));
    f = @(x,y) uex(x,y).*(1 - (gpp(x) + gp(x).^2 + gpp(y) + gp(y).^2));
    cs.prob = ebProblem([-0.25 0.25 -0.25 0.25], 1, 1, f, {'neumann', 0}, uex);

  case 'gauss52'
    % multi-peak periodic Gaussian, paper eq. (5.2), on (0,1)^2, Dirichlet
    sg = 0.25;
    xi = [0.4375 0.4375 0.5625]; yi = [0.5625 0.4375 0.4375];
    Ai = [1 -1 1];
    cs.prob = makeGauss52(sg, xi, yi, Ai);

  case 'threebumps'
    % three sharp isotropic Gaussians (for adaptive tagging tests)
    s = 0.02;
    xi = [0.30 0.66 0.49]; yi = [0.34 0.42 0.72];
    Ai = [1 0.8 -0.9];
    g  = @(x,y,a,b) exp(-((x-a).^2 + (y-b).^2)/(2*s^2));
    lg = @(x,y,a,b) g(x,y,a,b) .* (((x-a).^2 + (y-b).^2)/s^4 - 2/s^2);
    uex = @(x,y) Ai(1)*g(x,y,xi(1),yi(1)) + Ai(2)*g(x,y,xi(2),yi(2)) ...
               + Ai(3)*g(x,y,xi(3),yi(3));
    lap = @(x,y) Ai(1)*lg(x,y,xi(1),yi(1)) + Ai(2)*lg(x,y,xi(2),yi(2)) ...
               + Ai(3)*lg(x,y,xi(3),yi(3));
    f = @(x,y) uex(x,y) - lap(x,y);
    cs.prob = ebProblem([0 1 0 1], 1, 1, f, {'dirichlet', uex}, uex);

  otherwise
    error('mms_cases:name', 'unknown case %s', name);
end
end

function prob = makeGauss52(sg, xi, yi, Ai)
% u = sum Ai exp(-(sin^2(x-xi) + sin^2(y-yi))/(2 sg^2)),  f = u - lap u
  function v = uex(x, y)
    v = 0;
    for i = 1:3
      v = v + Ai(i)*exp(-(sin(x-xi(i)).^2 + sin(y-yi(i)).^2)/(2*sg^2));
    end
  end
  function v = ff(x, y)
    v = 0;
    for i = 1:3
      gi = exp(-(sin(x-xi(i)).^2 + sin(y-yi(i)).^2)/(2*sg^2));
      gpx = -sin(2*(x-xi(i)))/(2*sg^2);  gppx = -cos(2*(x-xi(i)))/sg^2;
      gpy = -sin(2*(y-yi(i)))/(2*sg^2);  gppy = -cos(2*(y-yi(i)))/sg^2;
      v = v + Ai(i)*gi.*(1 - (gppx + gpx.^2 + gppy + gpy.^2));
    end
  end
prob = ebProblem([0 1 0 1], 1, 1, @ff, {'dirichlet', @uex}, @uex);
end
