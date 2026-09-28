function testEbNeumannSemantics
% Regression checks for reaction, gauge and source-projection semantics.
here = fileparts(mfilename('fullpath')); addpath(fullfile(here,'..','src'));
H = ebSetup(ebProblem([0 1 0 1],1,1,0,{'neumann',0},0), ...
    ebOptions('BaseCells',[8 8],'Verbose',0));
H.lev{1}.Q{1}(:) = 7;
e = ebErrorNorms(H);
assert(abs(e.linf-7)<1e-14 && e.meanShift==0, ...
    'Positive reaction must retain the absolute constant error.');
for nq = [1 3]
  H = ebSetup(ebProblem([0 1 0 1],1,0,0,{'neumann',0},7), ...
    ebOptions('BaseCells',[8 8],'SourceQuad',nq,'Verbose',0));
  e = ebErrorNorms(H);
  assert(e.linf<1e-14 && abs(e.meanShift+7)<1e-14, ...
      'Singular error must be gauge-aligned for every SourceQuad.');
end
% Small coefficients are well conditioned when scaled together.  They
% must not be classified as singular by an absolute reaction threshold.
prob = ebProblem([0 1 0 1],1e-18,1e-14,1e-14,{'neumann',0},1);
sol = ebSolve(prob,'BaseCells',[8 8],'Verbose',0,'RelTol',1e-10);
assert(~sol.H.coarsest.singular && sol.err.linf<1e-11, ...
    'A small positive reaction was treated as a nullspace.');

% u=x^4, -Delta u=-12x^2, outward g=4 on east and zero elsewhere.
% The midpoint quadrature defect changes as a region is refined.  The
% final projected source must differ from raw samples by ONE constant.
bc = struct('west',{{'neumann',0}},'east',{{'neumann',4}}, ...
    'south',{{'neumann',0}},'north',{{'neumann',0}});
prob = ebProblem([0 1 0 1],1,0,@(x,y)-12*x.^2,bc,@(x,y)x.^4);
H = ebSetup(prob,ebOptions('BaseCells',[16 16],'Verbose',0));
H = ebAddLevel(H,[5 12 5 12]);
H = ebAddLevel(H,[12 20 12 20]);
shifts = []; rawSum=0;
for k=H.rootIdx:H.nlev
  lev=H.lev{k};
  for p=1:lev.np
    b=lev.box(p,:); [x,y]=ndgrid(((b(1):b(2))-.5)*lev.h(1), ...
                               ((b(3):b(4))-.5)*lev.h(2));
    raw=prob.f(x,y); vis=~lev.covered{p};
    delta=raw-lev.FTrue{p}; shifts=[shifts;delta(vis)]; %#ok<AGROW>
    rawSum=rawSum+sum(raw(vis))*prod(lev.h);
  end
end
assert(max(shifts)-min(shifts)<1e-13,'Regridding mixed projected/raw sources.');
assert(abs(mean(shifts)-(rawSum+4))<1e-13,'Projection sign or volume is wrong.');
H2 = ebCompatProject(H);
assert(abs(H2.compatibility.sourceShift-H.compatibility.sourceShift)<1e-13, ...
    'Repeated compatibility projection must be idempotent.');
H=ebSolveCycles(H);
assert(H.stats.stages{end}.finalRel<1e-8,'Projected Neumann solve stalled.');
fprintf('testEbNeumannSemantics PASSED: gauge, scaled reaction, two regrids; finalRel=%.3e\n', ...
        H.stats.stages{end}.finalRel);
end

