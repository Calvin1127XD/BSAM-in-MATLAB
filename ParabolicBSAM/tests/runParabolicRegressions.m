function runParabolicRegressions
%RUNPARABOLICREGRESSIONS Nonautonomous CN order and composite conservation.
here=fileparts(mfilename('fullpath'));
run(fullfile(here,'..','src','ebpAddpaths.m'));
dts=[0.1 0.05 0.025 0.0125];
% Spatially constant solution isolates time error for C(t)=1+t.
p=ebpProblem([0 1 0 1],1,@(x,y,t) 1+t,@(x,y,t) t*exp(-t), ...
  {'neumann',0},@(x,y) ones(size(x)),@(x,y,t) exp(-t)*ones(size(x)));
e=zeros(size(dts));
for j=1:numel(dts)
  s=ebpSolve(p,'Scheme','cn','Dt',dts(j),'TFinal',0.4, ...
    'BaseCells',[8 8],'MaxLevels',1,'CoarsestCells',4, ...
    'StepVerbose',0,'RegridEvery',0,'RelTol',1e-12);
  q=s.H.lev{s.H.nlev}.Q{1}(2:end-1,2:end-1);
  e(j)=max(abs(q-exp(-0.4)),[],'all');
end
rates=log2(e(1:end-1)./e(2:end));
assert(all(rates>1.8),'parabolic_regression:cnOrder', ...
  'CN must retain second order for a time-dependent reaction.');
fprintf('PASS CN time-dependent C: errors %s; rates %s\n',mat2str(e,6),mat2str(rates,5));
% D(t), nonzero time-dependent boundary data, exact quadratic spatial
% profile. Same-grid Cauchy differences isolate time convergence from
% the fixed spatial boundary truncation error.
u=@(x,y,t) exp(-t)*(1+x.^2+y.^2);
f=@(x,y,t) -u(x,y,t)-4*(1+t)*exp(-t);
p=ebpProblem([0 1 0 1],@(x,y,t) 1+t,0,f,{'dirichlet',u}, ...
  @(x,y) u(x,y,0),u);
U=cell(size(dts));
for j=1:numel(dts)
  s=ebpSolve(p,'Scheme','cn','Dt',dts(j),'TFinal',0.4, ...
    'BaseCells',[16 16],'MaxLevels',1,'CoarsestCells',4, ...
    'StepVerbose',0,'RegridEvery',0,'RelTol',1e-12);
  U{j}=s.H.lev{s.H.nlev}.Q{1}(2:end-1,2:end-1);
end
d=zeros(1,numel(dts)-1);
for j=1:numel(d), d(j)=max(abs(U{j}-U{j+1}),[],'all'); end
rates=log2(d(1:end-1)./d(2:end));
assert(rates(end)>1.8,'parabolic_regression:cnDiffusionOrder', ...
  'CN must retain second order with time-dependent D and boundary data.');
fprintf('PASS CN time-dependent D/BC: Cauchy %s; rates %s\n',mat2str(d,6),mat2str(rates,5));
% A fixed local fine patch crosses a nonconstant heat profile. Both the
% explicit and implicit diffusion operators must telescope globally.
p=ebpProblem([0 1 0 1],@(x,y) 1+0.2*x,0,0,{'neumann',0}, ...
  @(x,y) 1+0.2*cos(2*pi*x).*cos(2*pi*y));
for regrid=[0 2]
  s=ebpSolve(p,'Scheme','cn','Dt',0.002,'TFinal',0.02, ...
    'BaseCells',[16 16],'MaxLevels',2,'CoarsestCells',4, ...
    'Indicator','fun','TagFunction',@(x,y) x>0.25 & x<0.625 & y>0.25 & y<0.75, ...
    'TagBuffer',1,'StepVerbose',0,'RegridEvery',regrid,'RelTol',1e-11);
  assert(numel(s.levels)==2,'parabolic_regression:mesh','Expected a coarse-fine interface.');
  drift=max(abs(s.mass-s.mass(1)));
  assert(drift<2e-10,'parabolic_regression:cnMass', ...
    'CN heat mass drift exceeds solve tolerance on the composite mesh.');
  fprintf('PASS CN composite mass: regrid=%d drift=%.3e rel=%.3e\n',regrid,drift,max(s.relres));
end
for bad={0,-0.1,NaN,Inf}
  rejected=false;
  try, ebpOptions('Dt',bad{1},'TFinal',1); catch, rejected=true; end
  assert(rejected,'parabolic_regression:time','Invalid Dt was accepted.');
end
fprintf('ALL PARABOLIC REGRESSIONS PASSED\n');
end
