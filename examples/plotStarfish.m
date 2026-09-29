function plotStarfish(sol,aux,result,out,Xc,Yc,difference)
%PLOTSTARFISH Export the solved field and the actual adaptive hierarchy.
% All figures are derived from sol.H. No illustrative/synthetic mesh is used.
H=sol.H; V=aux.vertices; nA=numel(sol.levels);
colors=[0.22 0.32 0.53; 0.08 0.53 0.69; 0.08 0.65 0.52; ...
        0.92 0.60 0.16; 0.80 0.24 0.25];
prefix='starfish'; if strcmp(result.mode,'quick'), prefix='starfish-quick'; end
x=linspace(-2,2,361); [X,Y]=ndgrid(x,x); U=ebSample(H,X,Y);

fig=newFigure([1500 660]);
t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t); surf(ax,X,Y,U,'EdgeColor','none','FaceColor','interp'); hold(ax,'on');
z=ebSample(H,V(:,1),V(:,2));
plot3(ax,V(:,1),V(:,2),z+0.035,'Color',[0.15 0.15 0.15],'LineWidth',1.5);
colormap(ax,turbo(256)); cb=colorbar(ax); cb.Label.String='u_\epsilon';
view(ax,[-43 35]); axis(ax,'tight'); grid(ax,'on');
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'u_\epsilon');
title(ax,'Diffuse-domain solution','FontWeight','normal');
ax=nexttile(t); hold(ax,'on');
for L=1:nA, patchLevel(ax,H,L,colors(L,:),0.30); end
plot(ax,V(:,1),V(:,2),'k-','LineWidth',1.8);
axis(ax,'equal'); axis(ax,[-2 2 -2 2]); xlabel(ax,'x'); ylabel(ax,'y');
title(ax,sprintf('%d refinement levels · %s active cells',nA,comma(result.visibleCells)), ...
    'FontWeight','normal');
handles=gobjects(1,nA);
for L=1:nA, handles(L)=plot(ax,nan,nan,'s','MarkerFaceColor',colors(L,:),'Color',colors(L,:)); end
legend(ax,handles,arrayfun(@(L) sprintf('L%d',L),1:nA,'UniformOutput',false), ...
    'Location','southoutside','Orientation','horizontal','Box','off');
title(t,'The starfish transmission problem','FontSize',21,'FontWeight','normal');
subtitle(t,sprintf('\\epsilon = 0.05   |   finest spacing 4/%d   |   %.1f%% of the uniform cell count', ...
    result.finestEquivalent,100*result.cellFraction),'FontSize',12);
saveFigure(fig,out,[prefix '-overview']);

fig=newFigure([1600 860]);
t=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
for L=1:nA
    ax=nexttile(t); hold(ax,'on');
    patchLevel(ax,H,L,colors(L,:),0.22);
    plot(ax,V(:,1),V(:,2),'k-','LineWidth',1.0);
    axis(ax,'equal'); axis(ax,[-2 2 -2 2]); xlabel(ax,'x'); ylabel(ax,'y');
    title(ax,sprintf('Level %d · h = 4/%d',L,sol.levels(L).n(1)),'FontWeight','normal');
    subtitle(ax,sprintf('Patches: %d · stored cells: %s',sol.levels(L).np,comma(sol.levels(L).cells)), ...
        'FontSize',10);
end
ax=nexttile(t,6); hold(ax,'on');
window=[-0.19 0.19 1.50 1.87];
plotCellEdges(ax,H,colors,window);
plot(ax,V(:,1),V(:,2),'k-','LineWidth',2.0);
axis(ax,'equal'); axis(ax,window); xlabel(ax,'x'); ylabel(ax,'y');
title(ax,'A closer look at the upper tip','FontWeight','normal');
subtitle(ax,'Actual cell edges; color denotes refinement level','FontSize',10);
title(t,'Block-structured adaptive grid hierarchy','FontSize',21,'FontWeight','normal');
subtitle(t,'Rectangles show solver patches. Refinement follows the diffuse interface.', 'FontSize',12);
saveFigure(fig,out,[prefix '-hierarchy']);

fig=newFigure([1400 570]);
t=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(t); hist=sol.stats.stages{end}.resHist;
semilogy(ax,0:numel(hist)-1,hist/hist(1),'-o','Color',colors(2,:), ...
    'LineWidth',2,'MarkerSize',5,'MarkerFaceColor','w'); hold(ax,'on');
yline(ax,1e-10,'--','Requested tolerance','Color',[0.4 0.4 0.4], ...
    'LabelHorizontalAlignment','center');
% Fit log(r_k) = a + k*log(rho) over the last five residual samples.
% The fitted rho is the per-V-cycle residual reduction factor; R^2 measures
% goodness of fit in log-residual space over those same five samples.
if numel(hist)>=5
    tail=numel(hist)-4:numel(hist);
    logResidual=log(hist(tail));
    fit=polyfit(tail-1,logResidual,1);
    rho=exp(fit(1));
    sse=sum((logResidual-polyval(fit,tail-1)).^2);
    sst=sum((logResidual-mean(logResidual)).^2);
    rSquared=NaN; % R^2 is undefined for a constant response.
    if sst>0, rSquared=1-sse/sst; end
    text(ax,0.97,0.96,sprintf('\\rho = %.4f\nR^2 = %.5f',rho,rSquared),'Units','normalized', ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'Interpreter','tex','FontSize',16,'BackgroundColor','white', ...
        'EdgeColor',[0.4 0.4 0.4],'Margin',6);
end
grid(ax,'on'); xlabel(ax,'V-cycle'); ylabel(ax,'||r_k||_\infty / ||r_0||_\infty');
title(ax,sprintf('Algebraic convergence · %d V-cycles',result.vcycles),'FontWeight','normal');
ax=nexttile(t);
if nargin>=7
    imagesc(ax,Xc(:,1),Yc(1,:),difference.'); set(ax,'YDir','normal'); hold(ax,'on');
    plot(ax,V(:,1),V(:,2),'k-','LineWidth',1);
    axis(ax,'equal'); axis(ax,[-2 2 -2 2]); colorbar(ax); colormap(ax,parula(256));
    title(ax,'AMR minus uniform 1024^2','FontWeight','normal');
    subtitle(ax,'Uniform reference averaged onto each visible AMR cell','FontSize',10);
else
    imagesc(ax,x,x,aux.phi(X,Y).'); set(ax,'YDir','normal'); hold(ax,'on');
    plot(ax,V(:,1),V(:,2),'k-','LineWidth',1); axis(ax,'equal'); axis(ax,[-2 2 -2 2]);
    colorbar(ax); colormap(ax,parula(256)); title(ax,'Phase field \phi_\epsilon','FontWeight','normal');
end
xlabel(ax,'x'); ylabel(ax,'y');
title(t,'Numerical checks at fixed interface width','FontSize',20,'FontWeight','normal');
saveFigure(fig,out,[prefix '-validation']);
end

function fig=newFigure(sz)
fig=figure('Visible','off','Color','w','Position',[80 80 sz]);
set(fig,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',12, ...
    'DefaultAxesLineWidth',0.8,'DefaultAxesBox','on', ...
    'DefaultAxesXColor',[0.15 0.18 0.22],'DefaultAxesYColor',[0.15 0.18 0.22]);
end

function patchLevel(ax,H,L,color,alpha)
lev=H.lev{H.rootIdx+L-1}; dom=H.prob.domain;
for p=1:lev.np
    b=lev.box(p,:); x0=dom(1)+(b(1)-1)*lev.h(1); x1=dom(1)+b(2)*lev.h(1);
    y0=dom(3)+(b(3)-1)*lev.h(2); y1=dom(3)+b(4)*lev.h(2);
    patch(ax,[x0 x1 x1 x0],[y0 y0 y1 y1],(1-alpha)+alpha*color, ...
        'EdgeColor',color,'LineWidth',0.55);
end
end

function plotCellEdges(ax,H,colors,w)
dom=H.prob.domain;
for k=H.rootIdx:H.nlev
    lev=H.lev{k}; L=k-H.rootIdx+1; hx=lev.h(1); hy=lev.h(2);
    edgesX=[]; edgesY=[];
    for p=1:lev.np
        b=lev.box(p,:);
        [I,J]=ndgrid(b(1):b(2),b(3):b(4));
        x0=dom(1)+(I-1)*hx; y0=dom(3)+(J-1)*hy;
        vis=~lev.covered{p} & x0+hx>=w(1) & x0<=w(2) & y0+hy>=w(3) & y0<=w(4);
        xx=x0(vis).'; yy=y0(vis).'; z=nan(size(xx));
        edgesX=[edgesX,reshape([xx;xx+hx;xx+hx;xx;xx;z],1,[])]; %#ok<AGROW>
        edgesY=[edgesY,reshape([yy;yy;yy+hy;yy+hy;yy;z],1,[])]; %#ok<AGROW>
    end
    plot(ax,edgesX,edgesY,'Color',colors(L,:),'LineWidth',0.45);
end
end

function saveFigure(fig,out,name)
exportgraphics(fig,fullfile(out,[name '.png']),'Resolution',170,'BackgroundColor','white');
% Rasterize the dense surface in its PDF to avoid a huge vector mesh.
kind='vector'; if contains(name,'overview'), kind='image'; end
exportgraphics(fig,fullfile(out,[name '.pdf']),'ContentType',kind, ...
    'Resolution',200,'BackgroundColor','white');
close(fig);
end

function s=comma(n)
s=regexprep(sprintf('%d',n),'\d(?=(\d{3})+$)','$0,');
end
