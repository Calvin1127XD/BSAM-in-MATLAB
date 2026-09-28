function [prob, aux] = ddTransmissionProblem(geom, epsilon, box, data)
%DDTRANSMISSIONPROBLEM Two-sided diffuse-domain transmission coefficients.
% -div(D grad u) + (c + kappa*delta)u = h(1-phi)+q*phi-g*delta,
% phi=(1+tanh(d/epsilon))/2, delta=sech(d/epsilon)^2/(2*epsilon).
% geom.signedDistance is positive INSIDE (the opposite of ddGeometry.psi).
% data: alpha>0, beta/gamma/kappa>=0, q/h/g scalar or @(x,y).
% Default outer condition is homogeneous Neumann. No exact solution is
% assumed; data.exact is optional and must match the actual forcing.
validateattributes(epsilon,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(data.alpha,{'numeric'},{'scalar','real','finite','positive'});
for f={'beta','gamma','kappa'}
    validateattributes(data.(f{1}),{'numeric'},{'scalar','real','finite','nonnegative'});
end
phi=@(X,Y) (1+tanh(geom.signedDistance(X,Y)/epsilon))/2;
delta=@(X,Y) sech(geom.signedDistance(X,Y)/epsilon).^2/(2*epsilon);
q=tofun(data.q); h=tofun(data.h); g=tofun(data.g);
D=@(X,Y) data.alpha+(1-data.alpha)*phi(X,Y);
C=@(X,Y) data.beta+(data.gamma-data.beta)*phi(X,Y)+data.kappa*delta(X,Y);
F=@(X,Y) h(X,Y).*(1-phi(X,Y))+q(X,Y).*phi(X,Y)-g(X,Y).*delta(X,Y);
bc={'neumann',0}; exact=[];
if isfield(data,'OuterBC'), bc=data.OuterBC; end
if isfield(data,'exact'), exact=data.exact; end
prob=ebProblem(box,D,C,F,bc,exact);
aux=struct('phi',phi,'delta',delta,'epsilon',epsilon, ...
    'signedDistance',geom.signedDistance,'vertices',geom.vertices);
end

function f=tofun(v)
if isa(v,'function_handle'), f=v; else, f=@(X,Y) v+zeros(size(X)); end
end
