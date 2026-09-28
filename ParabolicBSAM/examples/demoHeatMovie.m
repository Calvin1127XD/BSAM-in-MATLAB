function demoHeatMovie()
%DEMOHEATMOVIE A rotating heat source on an adaptive, moving hierarchy:
%   the patches chase the source around the domain.  Writes
%   ../output/heat_movie.gif (+ .mp4) from per-step frames.

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'src', 'ebpAddpaths.m'));
outdir = fullfile(here, '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end
warning('off', 'eb_tag_cluster:dropthin');

% Gaussian heat source of width ss circling the domain center at radius
% Rc once per unit time; the solution starts at zero.
om = 2*pi; Rc = 0.22; ss = 0.03;
cx = @(t) 0.5 + Rc*cos(om*t);
cy = @(t) 0.5 + Rc*sin(om*t);
src = @(x,y,t) 25*exp(-((x-cx(t)).^2 + (y-cy(t)).^2)/(2*ss^2));
pprob = ebpProblem([0 1 0 1], 1, 0, src, {'dirichlet', 0}, @(x,y) 0*x);

% Frame grabber invoked by ebpSolve after every step; keeps every 4th
% step as an in-memory RGB frame (patch outlines + source marker).
frames = {};
every = 4;
  function grab(H, t, n)
    if mod(n, every) ~= 0, return; end
    fig = figure('Visible', 'off', 'Position', [80 80 720 600], 'Color', 'w');
    ax = axes(fig); hold(ax, 'on');
    dom = [0 1 0 1];
    for k = H.rootIdx:H.nlev
      lev = H.lev{k};
      for p = 1:lev.np
        b = lev.box(p, :);
        xf = dom(1) + (b(1)-1:b(2)) * lev.h(1);
        yf = dom(3) + (b(3)-1:b(4)) * lev.h(2);
        [Xg, Yg] = meshgrid(xf, yf);
        C = lev.Q{p}(2:end-1, 2:end-1);
        Cp = nan(size(Xg)); Cp(1:end-1, 1:end-1) = C.';
        surf(ax, Xg, Yg, zeros(size(Xg)), Cp, 'EdgeColor', 'none', ...
          'FaceColor', 'flat');
      end
    end
    nA = H.nlev - H.rootIdx + 1;
    cmapL = lines(max(nA, 2));
    for k = H.rootIdx+1:H.nlev
      lev = H.lev{k};
      for p = 1:lev.np
        b = lev.box(p, :);
        rectangle(ax, 'Position', [dom(1)+(b(1)-1)*lev.h(1), ...
          dom(3)+(b(3)-1)*lev.h(2), (b(2)-b(1)+1)*lev.h(1), ...
          (b(4)-b(3)+1)*lev.h(2)], ...
          'EdgeColor', cmapL(k-H.rootIdx+1, :), 'LineWidth', 1.5);
      end
    end
    plot(ax, cx(t), cy(t), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    view(ax, 2); axis(ax, 'equal'); axis(ax, dom);
    clim(ax, [0 0.06]); colormap(ax, parula); colorbar(ax);
    title(ax, sprintf('t = %.3f   (red x = source center; boxes = patches)', t));
    tmp = fullfile(outdir, 'frame_tmp.png');
    exportgraphics(fig, tmp, 'Resolution', 100);
    frames{end+1} = imread(tmp); %#ok<AGROW>
    close(fig);
  end

sol = ebpSolve(pprob, 'Scheme', 'bdf2', 'Dt', 5e-3, 'TFinal', 1.0, ...
  'RegridEvery', 4, 'StepVerbose', 40, 'Callback', @grab, ...
  'BaseCells', [64 64], 'MaxLevels', 3, ...
  'Indicator', 'ulap', 'TagThreshold', @(L) 2e-3/4^(L-1), ...
  'TagBuffer', 3, 'RelTol', 1e-8);

% normalize and write frames + gif + mp4.  Frame sides are padded up to
% multiples of 16 (H.264 macroblocks): MATLAB's MPEG-4 pipeline decodes
% non-16-multiple widths with a row shear.
hs = cellfun(@(f) size(f,1), frames); ws = cellfun(@(f) size(f,2), frames);
Hh = 16 * ceil(max(hs) / 16); Ww = 16 * ceil(max(ws) / 16);
fdir = fullfile(outdir, 'movie_frames');
if ~exist(fdir, 'dir'), mkdir(fdir); end
for i = 1:numel(frames)
  P = 255*ones(Hh, Ww, 3, 'uint8');
  P(1:size(frames{i},1), 1:size(frames{i},2), :) = frames{i};
  frames{i} = P;
  imwrite(P, fullfile(fdir, sprintf('frame_%03d.png', i)));
end
giff = fullfile(outdir, 'heat_movie.gif');
for i = 1:numel(frames)
  [A, map] = rgb2ind(frames{i}, 256, 'nodither');
  if i == 1
    imwrite(A, map, giff, 'gif', 'LoopCount', Inf, 'DelayTime', 0.12);
  else
    imwrite(A, map, giff, 'gif', 'WriteMode', 'append', 'DelayTime', 0.12);
  end
end
vw = VideoWriter(fullfile(outdir, 'heat_movie.mp4'), 'MPEG-4');
vw.FrameRate = 8; open(vw);
for i = 1:numel(frames), writeVideo(vw, frames{i}); end
close(vw);
delete(fullfile(outdir, 'frame_tmp.png'));

fprintf('\nmean V-cycles/step = %.2f, regrids = %d\n', ...
  mean(sol.cycles), sol.nregrids);
fprintf('movie written to %s (.gif/.mp4)\n', fullfile(outdir, 'heat_movie'));
end
