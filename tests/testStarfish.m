function testStarfish
%TESTSTARFISH Independent geometry and transmission-sign checks.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); bsamSetup;
g=ddStarfishGeometry;
rng(1127); X=4*rand(300,1)-2; Y=4*rand(300,1)-2;
v=g.vertices; brute=inf(size(X));
for k=1:size(v,1)-1
    a=v(k,:); e=v(k+1,:)-a;
    t=max(0,min(1,((X-a(1))*e(1)+(Y-a(2))*e(2))/sum(e.^2)));
    brute=min(brute,hypot(X-a(1)-t*e(1),Y-a(2)-t*e(2)));
end
d=g.signedDistance(X,Y); inside=inpolygon(X,Y,v(:,1),v(:,2));
assert(max(abs(abs(d)-brute))<5e-14,'Bounding-box distance disagrees with brute force.');
assert(all((d>=0)==inside),'Distance sign disagrees with inpolygon.');
assert(g.signedDistance(0,0)>0 && g.signedDistance(2,2)<0);
assert(max(abs(g.signedDistance(v(:,1),v(:,2))))<1e-13);
fprintf('PASS starfish distance: 300 independent point-to-segment/sign checks.\n');
% With u=2, choose q=gamma*u, h=beta*u, g=-kappa*u. This tests
% the interface source SIGN against an exact constant solution.
data=struct('alpha',3,'beta',2,'gamma',1,'kappa',0.01, ...
    'q',2,'h',4,'g',-0.02,'exact',@(x,y) 2+zeros(size(x)));
p=ddTransmissionProblem(g,0.05,[-2 2 -2 2],data);
assert(max(abs(p.f(X,Y)-2*p.C(X,Y)))<1e-13);
s=ebSolve(p,'BaseCells',32,'MaxLevels',3,'Indicator','fun', ...
    'TagFunction',@(x,y) abs(g.signedDistance(x,y))<0.2, ...
    'RelTol',1e-11,'Verbose',0);
assert(s.err.linf<1e-9 && s.stats.stages{end}.finalRel<1e-11);
fprintf('PASS transmission constant solution: Linf %.3e, residual %.3e.\n', ...
    s.err.linf,s.stats.stages{end}.finalRel);
fprintf('ALL STARFISH REGRESSIONS PASSED\n');
end
