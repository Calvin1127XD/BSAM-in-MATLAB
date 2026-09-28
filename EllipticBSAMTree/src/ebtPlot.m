function fig = ebtPlot(H, what, fname, titlestr)
%ebtPlot Headless-friendly visualization.
%   ebtPlot(H, 'mesh', fname [, title])     patch layout colored by level
%   ebtPlot(H, 'solution', fname [, title]) composite solution as a 3D
%                                           surface  z = u(x,y)
%   Saves with exportgraphics when fname is nonempty; figure is invisible.
%
%   The 'solution' mode samples the COMPOSITE solution on a uniform grid
%   (ebtSample picks the finest patch covering each point), so the surface
%   is a single clean sheet with no coarse/fine z-fighting.

if nargin < 4, titlestr = ''; end
fig = figure('Visible', 'off', 'Position', [100 100 820 660]);
ax = axes(fig);
hold(ax, 'on');
dom = H.prob.domain;
nA = H.nlev - H.rootIdx + 1;
cmapL = lines(max(nA, 2));

switch lower(what)
  case 'mesh'
    for k = H.rootIdx:H.nlev
      lev = H.lev{k};
      L = k - H.rootIdx + 1;
      for p = 1:lev.np
        b = lev.box(p, :);
        x0 = dom(1) + (b(1)-1)*lev.h(1);
        y0 = dom(3) + (b(3)-1)*lev.h(2);
        w = (b(2)-b(1)+1)*lev.h(1);
        ht = (b(4)-b(3)+1)*lev.h(2);
        rectangle(ax, 'Position', [x0 y0 w ht], ...
          'EdgeColor', cmapL(L,:), 'LineWidth', 1.4, ...
          'FaceColor', [cmapL(L,:) 0.06]);
      end
    end
    leg = arrayfun(@(L) sprintf('level %d', L), 1:nA, 'UniformOutput', false);
    for L = 1:nA
      plot(ax, NaN, NaN, '-', 'Color', cmapL(L,:), 'LineWidth', 2);
    end
    legend(ax, leg, 'Location', 'eastoutside');
    axis(ax, 'equal');
    axis(ax, [dom(1) dom(2) dom(3) dom(4)]);
    xlabel(ax, 'x'); ylabel(ax, 'y');

  case 'solution'
    % Sample at roughly the finest resolution, capped at 401 points per
    % side to keep figure export fast.
    hfin = H.lev{H.nlev}.h(1);
    NS = min(round((dom(2)-dom(1)) / hfin) + 1, 401);
    [Xq, Yq] = ndgrid(linspace(dom(1), dom(2), NS), ...
                      linspace(dom(3), dom(4), NS));
    U = ebtSample(H, Xq, Yq);
    surf(ax, Xq, Yq, U, 'EdgeColor', 'none');
    colormap(ax, parula(256));
    cb = colorbar(ax);
    cb.Label.String = 'u';
    view(ax, [-37 42]);
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    shading(ax, 'interp');
    xlim(ax, [dom(1) dom(2)]); ylim(ax, [dom(3) dom(4)]);
    xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'u(x,y)');
    grid(ax, 'on');

  otherwise
    error('ebt_plot:what', 'what must be mesh or solution');
end

if ~isempty(titlestr), title(ax, titlestr, 'Interpreter', 'none'); end
if nargin >= 3 && ~isempty(fname)
  exportgraphics(fig, fname, 'Resolution', 150);
  close(fig);
  fig = [];
end
end
