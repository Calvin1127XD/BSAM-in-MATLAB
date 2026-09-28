function runDdRegressions
%RUNDDREGRESSIONS Geometry normals, source signs, projection and validation.
here=fileparts(mfilename('fullpath'));
run(fullfile(here,'..','src','ddAddpaths.m'));
g=ddGeometry('polygon',5,1); v=g.vertices;
for k=1:5
  a=v(k,:); b=v(mod(k,5)+1,:); p=(a+b)/2; tangent=b-a;
  expected=[tangent(2),-tangent(1)]/norm(tangent);
  q=g.gpsi(p(1),p(2)); normal=[q{1},q{2}];
  assert(norm(normal-expected)<1e-12,'dd_regression:normal', ...
    'Polygon edge midpoint must have the outward unit normal.');
  % Signed-distance finite difference across a smooth edge.
  h=1e-6;
  fd=(g.psi(p(1)+h*expected(1),p(2)+h*expected(2))- ...
      g.psi(p(1)-h*expected(1),p(2)-h*expected(2)))/(2*h);
  assert(abs(fd-1)<1e-8,'dd_regression:normalDerivative');
end
fprintf('PASS polygon edge normals and independent signed-distance derivatives\n');
g=ddGeometry('circle',1); ep=0.2;
u=@(x,y) (x.^2+y.^2)/4;
data=struct('D',2,'C',1,'f',@(x,y) u(x,y)-2, ...
  'g',@(x,y) hypot(x,y),'exact',u);
[prob,aux]=ddProblem(g,ep,[-2 2 -2 2],'neumann',data,'PhiFloor',0);
x=0.8; y=0.3; h=1e-5;
% Independently differentiate flux Ddd * grad u; grad u=(x,y)/2.
fluxX=@(x,y) prob.D(x,y).*x/2; fluxY=@(x,y) prob.D(x,y).*y/2;
div=(fluxX(x+h,y)-fluxX(x-h,y)+fluxY(x,y+h)-fluxY(x,y-h))/(2*h);
r=prob.C(x,y)*u(x,y)-div-prob.f(x,y);
assert(abs(r)<2e-7,'dd_regression:fluxSign','DDM Neumann source sign is inconsistent.');
fprintf('PASS diffuse flux/source sign, finite-difference residual %.3e\n',r);
for bad={0,-1,NaN,Inf}
  rejected=false;
  try, ddProblem(g,bad{1},[-2 2 -2 2],'neumann',data); catch, rejected=true; end
  assert(rejected,'dd_regression:epsilon','Invalid epsilon accepted.');
end
% Projection preserves the composite integral on a dyadically aligned mesh.
p=ebProblem([0 1 0 1],1,1,1,{'neumann',0},@(x,y) ones(size(x)));
s=ebSolve(p,'BaseCells',[16 16],'MaxLevels',2,'Verbose',0, ...
  'Indicator','fun','TagFunction',@(x,y) x<0.5 & y<0.5);
[U,~,~]=ddToUniform(s.H,32); [m,vol]=ebCompositeMean(s.H);
assert(abs(sum(U,'all')/numel(U)-m*vol)<1e-12,'dd_regression:projection');
fprintf('ALL DIFFUSE-DOMAIN REGRESSIONS PASSED\n');
end
