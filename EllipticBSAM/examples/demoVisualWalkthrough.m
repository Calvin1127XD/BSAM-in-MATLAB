function demoVisualWalkthrough()
%demoVisualWalkthrough A guided, animated tour of how EBSAM works.
%
%   Run it (desktop MATLAB: plays like a slideshow on screen; headless
%   `matlab -batch demoVisualWalkthrough`: renders off-screen).  Either
%   way it writes, into  ../output/walkthrough/ :
%     walkthrough.gif   the slideshow as an animated GIF
%     walkthrough.mp4   the same as an MP4 movie
%     slide_##.png      the individual slides
%
%   The story, told on a gentle sine background with four sharp bumps:
%     1. start from a uniform root grid and solve with FAS V-cycles;
%     2. measure WHERE the discretization is failing with the undivided
%        Laplacian curvature indicator  |D^2_h u|;
%     3. flag cells where |tau| is large, buffer them, and cluster the
%        flags into rectangular patches (Berger-Rigoutsos);
%     4. refine the patches 2x, prolong the solution as the initial
%        guess, and re-solve the COMPOSITE system (QIF ghost coupling +
%        flux-balanced mass corrections at every coarse-fine interface);
%     5. repeat until the requested depth - error drops ~4x per level.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
outdir = fullfile(here, '..', 'output', 'walkthrough');
if ~exist(outdir, 'dir'), mkdir(outdir); end
interactive = usejava('desktop');

% ---------- the problem: gentle sine background + four sharp bumps -----
%   u(x,y) = 0.1 sin(pi*x) sin(2*pi*y)
%            + sum_i A_i exp( -((x-a_i)^2+(y-b_i)^2) / (2 sigma_i^2) )
%   on [0,1]^2.  The sine part vanishes on the boundary EXACTLY
%   (sin(pi*x) = 0 at x = 0,1; sin(2*pi*y) = 0 at y = 0,1), so u satisfies
%   HOMOGENEOUS DIRICHLET data u = 0 on the box.  The four Gaussian bumps
%   sit well inside the domain and have DIFFERENT widths sigma_i, narrow
%   enough that the widest decays to ~3e-9 on the boundary -- negligible
%   against the discretization error.  We therefore solve the homogeneous
%   Dirichlet problem  -div(grad u) + u = f  with f built from u.  This is
%   the same manufactured solution used in the convergence study: a gentle,
%   low-amplitude background carrying a few sharp, well-localized features
%   -- the regime where block-structured AMR pays off.
xi = [0.50 0.50 0.30 0.80];          % bump centres (x)
yi = [0.50 0.55 0.30 0.80];          % bump centres (y)
si = [0.08 0.04 0.04 0.02];          % FOUR different widths
Ai = [1.00 5.00 -1.00 0.50];         % bump amplitudes
g  = @(x,y,a,b,s) exp(-((x-a).^2 + (y-b).^2) / (2*s^2));
lg = @(x,y,a,b,s) g(x,y,a,b,s) .* (((x-a).^2 + (y-b).^2)/s^4 - 2/s^2);
bg  = @(x,y) 0.1 * sin(pi*x) .* sin(2*pi*y);
lbg = @(x,y) -0.5*pi^2 * sin(pi*x) .* sin(2*pi*y);
uex = @(x,y) bg(x,y) ...
   + Ai(1)*g(x,y,xi(1),yi(1),si(1)) + Ai(2)*g(x,y,xi(2),yi(2),si(2)) ...
   + Ai(3)*g(x,y,xi(3),yi(3),si(3)) + Ai(4)*g(x,y,xi(4),yi(4),si(4));
lap = @(x,y) lbg(x,y) ...
   + Ai(1)*lg(x,y,xi(1),yi(1),si(1)) + Ai(2)*lg(x,y,xi(2),yi(2),si(2)) ...
   + Ai(3)*lg(x,y,xi(3),yi(3),si(3)) + Ai(4)*lg(x,y,xi(4),yi(4),si(4));
f = @(x,y) uex(x,y) - lap(x,y);            % -div(grad u) + u = f, D=C=1
prob = ebProblem([0 1 0 1], 1, 1, f, {'dirichlet', 0}, uex);

% Undivided-Laplacian indicator: refine where the discrete curvature
% |D^2_h u| is large.  The per-level threshold scales like h^2 (1/4 per
% level) so a feature stays tagged as it is refined; the sharp bumps far
% exceed it and drive the refinement, while the gently-curved sine
% background is already well resolved on the root grid.
opts = ebOptions('BaseCells', [64 64], 'MaxLevels', 4, ...
  'Indicator', 'ulap', 'TagThreshold', @(L) 0.002 / 4^(L-1), 'TagBuffer', 3, ...
  'Verbose', 0, 'RelTol', 1e-14, 'MaxVCycles', 40, 'MinFillRatio', 0.90, 'MinBoxWidth', 5);
H = ebSetup(prob, opts);
nb = opts.BaseCells(1);          % root cells per side (drives the captions)
nRounds = opts.MaxLevels - 1;    % refinement rounds = levels beyond the root

% S also carries the running convergence accumulators that feed the
% right-hand statistics dashboard (redrawn on every data slide).
S = struct('frames', {{}}, 'dwell', [], 'outdir', outdir, ...
           'interactive', interactive, 'count', 0, ...
           'lvl', [], 'elinf', [], 'el2', [], 'cells', [], ...
           'vc', [], 'fac', []);

% ======================= SLIDE 1: title ================================
[fig, ax] = newslide(interactive, ...
  'EBSAM — a guided tour of block-structured adaptive multigrid', ...
  {'The problem:   -\nabla\cdot(\nabla u) + u = f   on [0,1]^2,  homogeneous Dirichlet data.', ...
   'u = 0.1 sin(\pix)sin(2\piy) + four Gaussian bumps of DIFFERENT widths \sigma\in[0.02,0.08] —', ...
   sprintf('the gentle background is easy, but the tall, sharp bumps are far too steep for a uniform %d^2 grid.', nb)});
[Xs, Ys] = ndgrid(linspace(0, 1, 301));
surf(ax, Xs, Ys, uex(Xs, Ys), 'EdgeColor', 'none');
view(ax, [37 42]); camlight(ax, 'headlight'); lighting(ax, 'gouraud');
shading(ax, 'interp'); colormap(ax, parula(256));
xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'u(x,y)'); grid(ax, 'on');
title(ax, 'the exact solution u(x,y) = 0.1 sin(\pix)sin(2\piy) + 4 bumps');
S = addSlide(S, fig, 3.2, false);

% ======================= SLIDE 2: root grid ============================
[fig, ax] = newslide(interactive, sprintf('Step 1 — the uniform root grid (%d^2)', nb), ...
  {'Everything starts on a coarse, uniform "root" grid.  Below it (not shown) the solver keeps', ...
   'even coarser full-domain grids — the classic multigrid ladder ending in a tiny direct solve.', ...
   sprintf('Cell size h = 1/%d: the gentle background is well resolved, but the sharp bumps still curve too fast.', nb)});
drawGrid(ax, H, H.rootIdx);
drawBumps(ax, xi, yi);
S = addSlide(S, fig, 2.6);

% ======================= SLIDE 3: root solve ===========================
H = ebSolveCycles(H);
st = H.stats.stages{end};
err0 = ebErrorNorms(H);
S = pushStats(S, 1, err0, st);
[fig, ax] = newslide(interactive, 'Step 2 — solve on the root grid', ...
  {sprintf('FAS V(3,3) cycles with red-black Gauss-Seidel: %d cycles, residual reduction %.0e,', ...
     st.vcycles, st.finalRel), ...
   sprintf('a geometric factor of %.3f per cycle.  The solver is happy --- but the solution is not:', st.avgFactor), ...
   'the bumps are smeared.  Algebraic convergence says nothing about DISCRETIZATION error.'});
drawSolution(ax, H, false);
S = addSlide(S, fig, 2.8);

stencilShown = false;
for round = 1:nRounds
  Lcur = H.nlev - H.rootIdx + 1;
  [boxes, info] = ebTagCluster(H);
  if isempty(boxes), break; end

  % ---------------- indicator slide ------------------------------------
  [fig, ax] = newslide(interactive, ...
    sprintf('Step 3 — where is the grid failing?  (curvature, level %d)', Lcur), ...
    {'The undivided Laplacian |D^2_h u| = ((\delta^2_x u)^2 + (\delta^2_y u)^2)^{1/2} measures how', ...
     'sharply the solution curves relative to the cell size.  Where it exceeds the per-level', ...
     'threshold the grid is too coarse (log_{10} scale --- the sharp bumps blaze; the gentle sine glows faintly).'});
  drawIndicator(ax, H, info);
  S = addSlide(S, fig, 2.8);

  % ---------------- tags slide ------------------------------------------
  [fig, ax] = newslide(interactive, ...
    sprintf('Step 4 — flag and buffer (level %d)', Lcur), ...
    {sprintf('Cells with |D^2_h u| above the threshold are flagged, then padded by %d safety cells', H.opts.TagBuffer), ...
     'and clipped to the proper-nesting region (a child patch must sit fully inside its parent', ...
     'level, so every interface stencil stays valid).'});
  drawSolution(ax, H, true);
  drawTags(ax, H, info.T);
  S = addSlide(S, fig, 2.6);

  % ---------------- boxes slide ------------------------------------------
  [fig, ax] = newslide(interactive, ...
    sprintf('Step 5 — cluster the flags (Berger–Rigoutsos, level %d)', Lcur), ...
    {'Signature splitting carves the flagged region into a few well-shaped rectangles:', ...
     'split at holes/inflections, keep a minimum width, merge or drop slivers.', ...
     'These boxes, refined 2x, become the next level''s patches.'});
  drawSolution(ax, H, true);
  drawTags(ax, H, info.T);
  drawBoxes(ax, H, boxes, [0.85 0.1 0.1]);
  S = addSlide(S, fig, 2.6);

  % ---------------- refine + guess ---------------------------------------
  H = ebAddLevel(H, boxes);
  [fig, ax] = newslide(interactive, ...
    sprintf('Step 6 — refine: the %d-level composite mesh', H.nlev - H.rootIdx + 1), ...
    {'Each box is refined 2x (cells stay aligned with their parents).  The coarse solution is', ...
     'prolonged (bilinearly) onto the new patches as the initial guess.  Interfaces are', ...
     '1-conforming: levels always meet with a jump of exactly one.'});
  drawSolution(ax, H, true);
  drawMesh(ax, H);
  S = addSlide(S, fig, 2.6);

  % ---------------- stencil anatomy (once) -------------------------------
  if ~stencilShown
    S = slideStencil(S, interactive, H);
    stencilShown = true;
  end

  % ---------------- composite solve --------------------------------------
  H = ebSolveCycles(H);
  st = H.stats.stages{end};
  err = ebErrorNorms(H);
  S = pushStats(S, H.nlev - H.rootIdx + 1, err, st);
  [fig, ax] = newslide(interactive, ...
    sprintf('Step 7 — solve the %d-level composite', H.nlev - H.rootIdx + 1), ...
    {'One FAS V-cycle now sweeps ALL levels: fine ghosts come from QIF interpolation, and the', ...
     'coarse force is corrected by the flux balance at every interface cell — that is what keeps', ...
     sprintf('mass conserved to machine precision.  %d cycles, factor %.3f.   ||u - u_{exact}||_\\infty = %.2e.', ...
       st.vcycles, st.avgFactor, err.linf)});
  drawSolution(ax, H, false);
  drawMeshFloor(ax, H);
  S = addSlide(S, fig, 2.8);
end

% ======================= error map =====================================
err = ebErrorNorms(H);
[fig, ax] = newslide(interactive, 'Did it work?  The error map', ...
  {'log_{10}|u - u_{exact}| on the composite mesh: the bumps are now resolved on the finest', ...
   'patches and the error is uniformly small — second order in the finest h.', ...
   sprintf('Final:  ||e||_\\infty = %.2e,   ||e||_2 = %.2e.', err.linf, err.l2)});
drawError(ax, H, uex);
drawMesh(ax, H);
S = addSlide(S, fig, 2.8);

% ======================= summary =======================================
md = ebMassCheck(H);
lines = {sprintf('%-46s %s', 'levels (root + refinements):', ...
           sprintf('%d', H.nlev - H.rootIdx + 1))};
tot = 0;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  nc = 0;
  for p = 1:lev.np, nc = nc + numel(lev.FTrue{p}); end
  tot = tot + nc;
  lines{end+1} = sprintf('%-46s %d patch(es), %d cells, h = 1/%d', ...
    sprintf('  level %d:', k - H.rootIdx + 1), lev.np, nc, lev.n(1)); %#ok<AGROW>
end
ueq = H.lev{H.nlev}.n(1);
lines = [lines, { ...
  sprintf('%-46s %d  (uniform %d^2 would need %d)', 'total cells:', tot, ueq, ueq^2), ...
  sprintf('%-46s %.3f', 'V-cycle factor (depth-independent):', H.stats.stages{end}.avgFactor), ...
  sprintf('%-46s %.1e', 'mass-conservation defect:', md.rel), ...
  sprintf('%-46s %.2e', 'final error vs exact (L_inf):', err.linf), ...
  '', ...
  'The loop you just watched:  solve -> curvature -> flag -> cluster -> refine -> solve.', ...
  'Conservative coarse-fine coupling makes every composite solve behave like a', ...
  'single-grid multigrid solve — geometric convergence at any depth.'}];
[fig, axS] = newslide(interactive, 'Summary — what EBSAM just did', {});
set(axS, 'Visible', 'off');
annotation(fig, 'textbox', [0.08 0.10 0.86 0.74], 'String', lines, ...
  'EdgeColor', 'none', 'FontName', 'Menlo', 'FontSize', 12.5, ...
  'VerticalAlignment', 'middle', 'Interpreter', 'none');
S = addSlide(S, fig, 4.0, false);

% ======================= write the movie ===============================
writeMovie(S);
fprintf('\nWalkthrough written to:\n  %s\n  %s\n  (and slide_##.png)\n', ...
  fullfile(outdir, 'walkthrough.gif'), fullfile(outdir, 'walkthrough.mp4'));
end

% ========================================================================
function [fig, ax] = newslide(interactive, titleStr, captionLines)
if interactive, vis = 'on'; else, vis = 'off'; end
fig = figure('Visible', vis, 'Position', [60 60 980 660], 'Color', 'w');
try
  fig.Theme = 'light';   % headless -batch follows the OS dark theme
                         % (gray text, black axes); pin the light one
catch
end
annotation(fig, 'textbox', [0.04 0.905 0.92 0.075], 'String', titleStr, ...
  'EdgeColor', 'none', 'FontSize', 17, 'FontWeight', 'bold', ...
  'HorizontalAlignment', 'left');
if ~isempty(captionLines)
  annotation(fig, 'textbox', [0.045 0.012 0.92 0.14], 'String', captionLines, ...
    'EdgeColor', 'none', 'FontSize', 11.5, 'VerticalAlignment', 'top');
end
% main plotting axes; the right band x in [0.62,1] is reserved for the
% convergence dashboard drawn by drawDashboard (see addSlide).
ax = axes(fig, 'Position', [0.065 0.20 0.525 0.63]);
hold(ax, 'on');
end

% ------------------------------------------------------------------------
function S = addSlide(S, fig, dwell, showStats)
if nargin < 4, showStats = true; end
if showStats, drawDashboard(fig, S); end
S.count = S.count + 1;
tmp = fullfile(S.outdir, sprintf('slide_%02d.png', S.count));
exportgraphics(fig, tmp, 'Resolution', 110);
S.frames{end+1} = imread(tmp);
S.dwell(end+1) = dwell;
if S.interactive
  drawnow;
  pause(min(dwell, 2.0));
end
close(fig);
end

% ------------------------------------------------------------------------
function drawGrid(ax, H, k)
lev = H.lev{k};
dom = H.prob.domain;
n = lev.n;
for i = 0:n(1)
  x = dom(1) + i*lev.h(1);
  plot(ax, [x x], [dom(3) dom(4)], '-', 'Color', [0.75 0.75 0.78], ...
    'LineWidth', 0.4);
end
for j = 0:n(2)
  y = dom(3) + j*lev.h(2);
  plot(ax, [dom(1) dom(2)], [y y], '-', 'Color', [0.75 0.75 0.78], ...
    'LineWidth', 0.4);
end
finishAxes(ax, H);
end

function drawBumps(ax, xi, yi)
plot(ax, xi, yi, 'o', 'MarkerSize', 14, 'LineWidth', 1.6, ...
  'Color', [0.85 0.2 0.2]);
text(ax, xi(1)+0.03, yi(1), 'the bumps live here', 'Color', [0.85 0.2 0.2], ...
  'FontSize', 11);
end

% ------------------------------------------------------------------------
function drawSolution(ax, H, grayed)
%drawSolution  grayed=true : flat 2D top-down backdrop for tag/box/mesh
%   overlays (a map view of the domain).  grayed=false : a 3D surface
%   z = u(x,y) of the COMPOSITE solution, sampled on a uniform grid so it
%   is a single clean sheet (ebSample picks the finest patch per point).
dom = H.prob.domain;
if grayed
  for k = H.rootIdx:H.nlev
    lev = H.lev{k};
    for p = 1:lev.np
      b = lev.box(p, :);
      xf = dom(1) + (b(1)-1:b(2)) * lev.h(1);
      yf = dom(3) + (b(3)-1:b(4)) * lev.h(2);
      [Xg, Yg] = meshgrid(xf, yf);
      C = lev.Q{p}(2:end-1, 2:end-1);
      Cp = nan(size(Xg));
      Cp(1:end-1, 1:end-1) = C.';
      surf(ax, Xg, Yg, zeros(size(Xg)), Cp, 'EdgeColor', 'none', ...
        'FaceColor', 'flat');
    end
  end
  view(ax, 2);
  colormap(ax, gray(256));
  finishAxes(ax, H);
else
  hfin = H.lev{H.nlev}.h(1);
  NS = min(round((dom(2)-dom(1)) / hfin) + 1, 401);
  [Xq, Yq] = ndgrid(linspace(dom(1), dom(2), NS), ...
                    linspace(dom(3), dom(4), NS));
  U = ebSample(H, Xq, Yq);
  surf(ax, Xq, Yq, U, 'EdgeColor', 'none');
  colormap(ax, parula(256));
  finishAxes3(ax, H, 'u(x,y)');
end
end

% ------------------------------------------------------------------------
function drawIndicator(ax, H, info)
dom = H.prob.domain;
k = H.nlev;
lev = H.lev{k};
vmax = -inf;
for p = 1:lev.np
  vmax = max(vmax, max(info.indicator{p}(:)));
end
for p = 1:lev.np
  b = lev.box(p, :);
  xf = dom(1) + (b(1)-1:b(2)) * lev.h(1);
  yf = dom(3) + (b(3)-1:b(4)) * lev.h(2);
  [Xg, Yg] = meshgrid(xf, yf);
  C = log10(info.indicator{p} / vmax + 1e-8);
  Cp = nan(size(Xg));
  Cp(1:end-1, 1:end-1) = C.';
  surf(ax, Xg, Yg, zeros(size(Xg)), Cp, 'EdgeColor', 'none', ...
    'FaceColor', 'flat');
end
view(ax, 2);
colormap(ax, hot(256));
cb = colorbar(ax);
cb.Label.String = 'log_{10} |D^2_h u| (normalized)';
finishAxes(ax, H);
end

% ------------------------------------------------------------------------
function drawTags(ax, H, T)
dom = H.prob.domain;
n = size(T);
red = cat(3, ones(n(2), n(1)), zeros(n(2), n(1)), zeros(n(2), n(1)));
image(ax, 'XData', [dom(1) dom(2)], 'YData', [dom(3) dom(4)], ...
  'CData', red, 'AlphaData', 0.55 * double(T.'));
finishAxes(ax, H);
end

% ------------------------------------------------------------------------
function drawBoxes(ax, H, boxes, col)
lev = H.lev{H.nlev};
dom = H.prob.domain;
for b = 1:size(boxes, 1)
  bb = boxes(b, :);
  x0 = dom(1) + (bb(1)-1)*lev.h(1);
  y0 = dom(3) + (bb(3)-1)*lev.h(2);
  rectangle(ax, 'Position', [x0, y0, (bb(2)-bb(1)+1)*lev.h(1), ...
    (bb(4)-bb(3)+1)*lev.h(2)], 'EdgeColor', col, 'LineWidth', 2.2);
end
end

% ------------------------------------------------------------------------
function drawMesh(ax, H)
nA = H.nlev - H.rootIdx + 1;
cmapL = lines(max(nA, 2));
dom = H.prob.domain;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  L = k - H.rootIdx + 1;
  for p = 1:lev.np
    b = lev.box(p, :);
    x0 = dom(1) + (b(1)-1)*lev.h(1);
    y0 = dom(3) + (b(3)-1)*lev.h(2);
    rectangle(ax, 'Position', [x0, y0, (b(2)-b(1)+1)*lev.h(1), ...
      (b(4)-b(3)+1)*lev.h(2)], 'EdgeColor', cmapL(L,:), 'LineWidth', 1.6);
  end
end
end

% ------------------------------------------------------------------------
function drawError(ax, H, uex)
dom = H.prob.domain;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  for p = 1:lev.np
    b = lev.box(p, :);
    xc = dom(1) + ((b(1):b(2)) - 0.5) * lev.h(1);
    yc = dom(3) + ((b(3):b(4)) - 0.5) * lev.h(2);
    [XC, YC] = ndgrid(xc, yc);
    E = log10(abs(lev.Q{p}(2:end-1, 2:end-1) - uex(XC, YC)) + 1e-14);
    xf = dom(1) + (b(1)-1:b(2)) * lev.h(1);
    yf = dom(3) + (b(3)-1:b(4)) * lev.h(2);
    [Xg, Yg] = meshgrid(xf, yf);
    Cp = nan(size(Xg));
    Cp(1:end-1, 1:end-1) = E.';
    surf(ax, Xg, Yg, zeros(size(Xg)), Cp, 'EdgeColor', 'none', ...
      'FaceColor', 'flat');
  end
end
view(ax, 2);
colormap(ax, parula(256));
cb = colorbar(ax);
cb.Label.String = 'log_{10}|u - u_{exact}|';
finishAxes(ax, H);
end

% ------------------------------------------------------------------------
function finishAxes3(ax, H, zlab)
%finishAxes3  3D finish for a solution surface: domain limits, oblique
%   view, lighting and labels (NO axis-equal, which would flatten z).
dom = H.prob.domain;
xlim(ax, [dom(1) dom(2)]);
ylim(ax, [dom(3) dom(4)]);
view(ax, [-37 42]);
camlight(ax, 'headlight');
lighting(ax, 'gouraud');
shading(ax, 'interp');
grid(ax, 'on');
xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, zlab);
end

% ------------------------------------------------------------------------
function drawMeshFloor(ax, H)
%drawMeshFloor  Project the patch outlines (coloured by level) onto the
%   floor plane z = zmin of a 3D solution surface, so the AMR structure is
%   visible beneath the sharp features that triggered it.
zl = zlim(ax); z0 = zl(1);
nA = H.nlev - H.rootIdx + 1;
cmapL = lines(max(nA, 2));
dom = H.prob.domain;
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  L = k - H.rootIdx + 1;
  for p = 1:lev.np
    b = lev.box(p, :);
    x0 = dom(1) + (b(1)-1)*lev.h(1); x1 = dom(1) + b(2)*lev.h(1);
    y0 = dom(3) + (b(3)-1)*lev.h(2); y1 = dom(3) + b(4)*lev.h(2);
    plot3(ax, [x0 x1 x1 x0 x0], [y0 y0 y1 y1 y0], z0*ones(1,5), '-', ...
      'Color', cmapL(L,:), 'LineWidth', 1.6);
  end
end
end

% ------------------------------------------------------------------------
function S = pushStats(S, level, err, st)
%pushStats  Append one (finest-level) convergence record to the
%   accumulators that feed the right-hand dashboard.
S.lvl(end+1)   = level;
S.elinf(end+1) = err.linf;
S.el2(end+1)   = err.l2;
S.cells(end+1) = st.visibleCells;
S.vc(end+1)    = st.vcycles;
S.fac(end+1)   = st.avgFactor;
end

% ------------------------------------------------------------------------
function drawDashboard(fig, S)
%drawDashboard  Right-hand statistics panel, redrawn on every data slide
%   from the accumulators in S (each slide is a fresh figure).  Shows the
%   latest solver numbers and -- the headline -- the GLOBAL ERROR coming
%   down level by level as the refinement proceeds.
if isempty(S.lvl)
  txt = {'\bf Convergence tracker', '', '(root solve pending)'};
else
  i = numel(S.lvl);
  txt = {'\bf Convergence tracker', '', ...
    sprintf('finest level: %d', S.lvl(i)), ...
    sprintf('visible cells: %d', S.cells(i)), ...
    sprintf('V-cycles: %d  (factor %.3f)', S.vc(i), S.fac(i)), '', ...
    sprintf('L_\\infty error: %.2e', S.elinf(i)), ...
    sprintf('L_2 error: %.2e', S.el2(i))};
  if i >= 2
    txt{end+1} = '';
    txt{end+1} = sprintf('drop this level: %.1f\\times', ...
      S.elinf(i-1) / max(S.elinf(i), realmin));
  end
end
annotation(fig, 'textbox', [0.635 0.555 0.345 0.30], 'String', txt, ...
  'EdgeColor', [0.80 0.80 0.86], 'LineWidth', 0.8, ...
  'BackgroundColor', [0.975 0.975 1.0], 'FontSize', 11, ...
  'FontName', 'Menlo', 'Interpreter', 'tex', ...
  'VerticalAlignment', 'top', 'Margin', 8);

if isempty(S.lvl), return; end
axd = axes(fig, 'Position', [0.705 0.175 0.265 0.275]);
semilogy(axd, S.lvl, S.elinf, 'o-', 'LineWidth', 1.7, ...
  'Color', [0.85 0.20 0.20], 'MarkerFaceColor', [0.85 0.20 0.20]);
hold(axd, 'on');
semilogy(axd, S.lvl, S.el2, 's--', 'LineWidth', 1.3, ...
  'Color', [0.20 0.30 0.85], 'MarkerFaceColor', [0.20 0.30 0.85]);
if numel(S.lvl) >= 2
  ref = S.elinf(1) * 4.^-(S.lvl - S.lvl(1));   % ideal 2nd order: 4x/level
  semilogy(axd, S.lvl, ref, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
  legend(axd, {'\|e\|_\infty', '\|e\|_2', '4\times/level'}, ...
    'Location', 'southwest', 'FontSize', 8, 'Box', 'off');
end
grid(axd, 'on');
xlim(axd, [min(S.lvl)-0.3, max(S.lvl)+0.3]);
set(axd, 'XTick', unique(S.lvl));
xlabel(axd, 'refinement level');
ylabel(axd, 'global error');
title(axd, 'error vs refinement', 'FontWeight', 'normal', 'FontSize', 10);
end

% ------------------------------------------------------------------------
function S = slideStencil(S, interactive, H)
[fig, ax] = newslide(interactive, ...
  'Anatomy — how levels talk: the QIF ghost stencil', ...
  {'A fine patch needs GHOST values just outside its boundary.  Each ghost (square) is the', ...
   'quadratic along the diagonal through one parent-cell center and two interior fine cells:', ...
   'u_g = 8/15 u_C + 2/3 u_{F1} - 1/5 u_{F2}.   In the other direction, the coarse force is', ...
   'corrected so the coarse flux equals the sum of fine fluxes - that is mass conservation.'});
% schematic, in units of the fine cell h = 1
h = 1;
for i = -2:0                       % coarse grid (h_c = 2) on the left
  for j = -1:2
    rectangle(ax, 'Position', [2*i*h, 2*j*h - h, 2*h, 2*h], ...
      'EdgeColor', [0.55 0.55 0.6], 'LineWidth', 1.0);
  end
end
for i = 0:3                        % fine grid on the right
  for j = -2:5
    rectangle(ax, 'Position', [i*h, j*h, h, h], ...
      'EdgeColor', [0.3 0.45 0.8], 'LineWidth', 0.9);
  end
end
plot(ax, [0 0], [-3 7], '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.4);
text(ax, 0.08, 6.45, 'coarse–fine interface', 'FontSize', 11, 'Rotation', 0);
% ghost at fine row r=2 (centers at x=0.5): G=(-0.5, 1.5); C=(-1, 1);
% F1=(0.5, 2.5); F2=(1.5, 3.5)
G = [-0.5 1.5]; C = [-1 1]; F1 = [0.5 2.5]; F2 = [1.5 3.5];
plot(ax, [C(1) F2(1)], [C(2) F2(2)], '-', 'Color', [0.85 0.3 0.1], ...
  'LineWidth', 1.6);
plot(ax, C(1), C(2), 'o', 'MarkerSize', 13, 'MarkerFaceColor', [0.55 0.55 0.6], ...
  'Color', 'k');
plot(ax, G(1), G(2), 's', 'MarkerSize', 13, 'MarkerFaceColor', [0.95 0.8 0.2], ...
  'Color', 'k');
plot(ax, F1(1), F1(2), 'o', 'MarkerSize', 11, 'MarkerFaceColor', [0.3 0.45 0.8], ...
  'Color', 'k');
plot(ax, F2(1), F2(2), 'o', 'MarkerSize', 11, 'MarkerFaceColor', [0.3 0.45 0.8], ...
  'Color', 'k');
text(ax, C(1)-2.6, C(2)-0.05, 'u_C  (8/15)', 'FontSize', 12);
text(ax, G(1)-2.15, G(2)+0.55, 'ghost u_g', 'FontSize', 12, 'Color', [0.75 0.6 0]);
text(ax, F1(1)+0.25, F1(2), 'u_{F1}  (2/3)', 'FontSize', 12);
text(ax, F2(1)+0.25, F2(2), 'u_{F2}  (-1/5)', 'FontSize', 12);
% flux-balance annotation pointing at the highlighted coarse face
yf0 = -1;
plot(ax, [0 0], [yf0 yf0+2], '-', 'Color', [0.1 0.6 0.3], 'LineWidth', 5);
annotation(fig, 'textarrow', [0.60 0.45], [0.27 0.36], ...
  'String', {'flux balance: coarse face flux', '= sum of its fine face fluxes'}, ...
  'FontSize', 10.5);
axis(ax, 'equal');
axis(ax, [-6.2 4.4 -3.2 7.2]);
ax.XTick = []; ax.YTick = [];
S = addSlide(S, fig, 3.4, false);
end

% ------------------------------------------------------------------------
function finishAxes(ax, H)
dom = H.prob.domain;
axis(ax, 'equal');
axis(ax, [dom(1) dom(2) dom(3) dom(4)]);
end

% ------------------------------------------------------------------------
function writeMovie(S)
% normalize frame sizes (pad with white) to multiples of 16: MATLAB's
% MPEG-4 pipeline shears frames whose width is not a macroblock multiple
hs = cellfun(@(f) size(f, 1), S.frames);
ws = cellfun(@(f) size(f, 2), S.frames);
Hh = 16 * ceil(max(hs) / 16);
Ww = 16 * ceil(max(ws) / 16);
for i = 1:numel(S.frames)
  f = S.frames{i};
  P = 255 * ones(Hh, Ww, 3, 'uint8');
  P(1:size(f,1), 1:size(f,2), :) = f;
  S.frames{i} = P;
end

giffile = fullfile(S.outdir, 'walkthrough.gif');
for i = 1:numel(S.frames)
  [A, map] = rgb2ind(S.frames{i}, 256);
  if i == 1
    imwrite(A, map, giffile, 'gif', 'LoopCount', Inf, ...
      'DelayTime', S.dwell(i));
  else
    imwrite(A, map, giffile, 'gif', 'WriteMode', 'append', ...
      'DelayTime', S.dwell(i));
  end
end

vw = VideoWriter(fullfile(S.outdir, 'walkthrough.mp4'), 'MPEG-4');
vw.FrameRate = 30;
vw.Quality = 90;
open(vw);
for i = 1:numel(S.frames)
  nrep = max(1, round(S.dwell(i) * vw.FrameRate));
  for r = 1:nrep
    writeVideo(vw, S.frames{i});
  end
end
close(vw);
end
