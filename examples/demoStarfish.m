function result = demoStarfish(mode)
%DEMOSTARFISH Reproduce the Helmholtz starfish on an adaptive BSAM grid.
% demoStarfish('quick') : 64 -> 512, figures and solve diagnostics.
% demoStarfish('full')  : 64 -> 1024, plus uniform 512/1024 comparisons.
% The default 'full' run regenerates the README figures and JSON/CSV data.
if nargin<1, mode='full'; end
assert(ismember(string(mode),["quick","full"]),'Use quick or full.');
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); bsamSetup;
out=fullfile(root,'figures'); if ~exist(out,'dir'), mkdir(out); end
vr=fullfile(root,'verification'); if ~exist(vr,'dir'), mkdir(vr); end
g=ddStarfishGeometry(3000); epsilon=0.05;
data=struct('alpha',3,'beta',2,'gamma',1,'kappa',0.01, ...
    'q',@(x,y) 15-x.^2, 'h',@(x,y) 2.5*sin(x)+exp(cos(y)), 'g',4);
[prob,aux]=ddTransmissionProblem(g,epsilon,[-2 2 -2 2],data);
base=64; levels=4+strcmp(mode,'full');
tag=@(x,y,L) abs(g.signedDistance(x,y))<max(3*epsilon,3*4/(base*2^(L-1)));
timer=tic;
sol=ebSolve(prob,'BaseCells',base,'MaxLevels',levels, ...
    'Indicator','fun','TagFunction',tag,'TagBuffer',2, ...
    'RelTol',1e-10,'MaxVCycles',100,'Verbose',1);
elapsed=toc(timer); st=sol.stats.stages{end};
assert(st.finalRel<=1e-10,'Starfish solve did not meet its residual tolerance.');
assert(numel(sol.levels)==levels,'Requested refinement levels were not built.');
finest=base*2^(levels-1);
result=struct('mode',mode,'matlab',version,'epsilon',epsilon,'alpha',3, ...
    'beta',2,'gamma',1,'kappa',0.01,'boundaryPoints',3000, ...
    'baseCells',base,'finestEquivalent',finest,'levels',levels, ...
    'visibleCells',st.visibleCells,'storedAmrCells',st.totalCells, ...
    'uniformCells',finest^2,'cellFraction',st.visibleCells/finest^2, ...
    'vcycles',st.vcycles,'relativeResidual',st.finalRel, ...
    'balanceDefect',sol.mass.rel,'setupAndSolveSeconds',elapsed);
% Verify that the actual finest patches cover all cell centers in |d|<eps.
% Check at the finest resolution, independently of the tagging/clustering.
x=-2+((1:finest)-0.5)*4/finest; [X,Y]=ndgrid(x,x);
band=abs(g.signedDistance(X,Y))<epsilon;
cover=false(finest); bx=sol.levels(end).box;
for k=1:size(bx,1), cover(bx(k,1):bx(k,2),bx(k,3):bx(k,4))=true; end
result.interfaceBandCells=nnz(band);
result.uncoveredBandCells=nnz(band & ~cover);
assert(result.uncoveredBandCells==0,'Finest grid leaves part of the interface band uncovered.');
result.hOverEpsilon=4/finest/epsilon;
result.patchesPerLevel=[sol.levels.np];
result.cellsPerLevel=[sol.levels.cells];
if strcmp(mode,'full')
    samples=cell(1,2); uniformResidual=zeros(1,2); uniformSeconds=zeros(1,2);
    for j=1:2
        N=256*2^j; timer=tic;
        ref=ebSolve(prob,'BaseCells',N,'MaxLevels',1, ...
            'RelTol',1e-10,'MaxVCycles',100,'Verbose',1);
        uniformSeconds(j)=toc(timer); uniformResidual(j)=ref.stats.stages{end}.finalRel;
        assert(uniformResidual(j)<1e-10,'Uniform reference solve failed.');
        samples{j}=ddToUniform(ref.H,512);
        if j==2
            % Compare at native visible AMR cell centers by averaging the
            % uniform 1024 cell averages onto each composite cell footprint.
            Ufine=ref.H.lev{ref.H.nlev}.Q{1}(2:end-1,2:end-1);
            [result.relativeL2VsUniform1024,result.linfVsUniform1024,difference]=compareComposite(sol.H,Ufine);
        end
        clear ref
    end
    result.uniformGridSizes=[512 1024]; result.uniformRelativeResiduals=uniformResidual;
    result.uniformSetupAndSolveSeconds=uniformSeconds;
    result.uniform512To1024RelativeL2=norm(samples{1}-samples{2},'fro')/norm(samples{2},'fro');
    result.comparison='Uniform 1024 cell averages restricted to visible AMR cells; fixed epsilon.';
    assert(isscalar(result.linfVsUniform1024) && isfinite(result.linfVsUniform1024));
    assert(isfinite(result.relativeL2VsUniform1024) && result.relativeL2VsUniform1024<1e-3, ...
        'Adaptive/uniform difference exceeds the benchmark acceptance threshold.');
    result.acceptanceRelativeL2=1e-3;
    plotStarfish(sol,aux,result,out,X,Y,difference);
else
    plotStarfish(sol,aux,result,out);
end
name=['starfish-' mode];
fid=fopen(fullfile(vr,[name '.json']),'w'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true)); clear cleanup;
history=table((0:numel(st.resHist)-1).',st.resHist(:),st.resHist(:)/st.resHist(1), ...
    'VariableNames',{'vcycle','absoluteResidual','relativeToInitial'});
writetable(history,fullfile(vr,[name '-residual.csv']));
fprintf('STARFISH %s PASSED: %d active cells, %.1f%% of uniform, residual %.3e.\n', ...
    upper(mode),st.visibleCells,100*result.cellFraction,st.finalRel);
disp(result);
end

function [rel,linf,E]=compareComposite(H,U)
N=size(U,1); e2=0; r2=0; linf=0; E=nan(N);
for k=H.rootIdx:H.nlev
    lev=H.lev{k}; ratio=N/lev.n(1);
    assert(ratio==round(ratio));
    for p=1:lev.np
        b=lev.box(p,:); mx=b(2)-b(1)+1; my=b(4)-b(3)+1;
        rows=(b(1)-1)*ratio+1:b(2)*ratio; cols=(b(3)-1)*ratio+1:b(4)*ratio;
        fine=U(rows,cols);
        avg=reshape(mean(mean(reshape(fine,ratio,mx,ratio,my),1),3),mx,my);
        q=lev.Q{p}(2:end-1,2:end-1); vis=~lev.covered{p}; area=prod(lev.h);
        err=q(vis)-avg(vis);
        if isempty(err), continue; end
        e2=e2+sum(err.^2)*area; r2=r2+sum(avg(vis).^2)*area;
        linf=max(linf,max(abs(err),[],'all'));
        expanded=repelem(q-avg,ratio,ratio); mask=repelem(vis,ratio,ratio);
        block=E(rows,cols); block(mask)=expanded(mask); E(rows,cols)=block;
    end
end
rel=sqrt(e2/r2);
assert(all(isfinite(E),'all'),'Composite comparison does not tile the domain.');
end
