%runDdForced Forced-problem study (no manufactured solution anywhere):
%   the data f, g are prescribed directly, the sharp solution on the
%   pentagon has no closed form, and the reference is built by CAUCHY
%   self-convergence on a common projection grid (ddToUniform).
%
%   Problem: pentagon Omega_1 (circumradius 1.2), forcing f = 1 + x/2.
%     F1 (Dirichlet penalty, u = 0 on Gamma):
%        h-Cauchy at fixed eps  -> rate ~ 2 (discretization)
%        eps-Cauchy at fixed h  -> model self-convergence; the full-
%        interior rate is layer-limited (theory: the O(eps) deviation in
%        the O(eps)-wide kink layer contributes eps^1.5 in L2), while the
%        eroded interior (psi < -5 eps_max) shows the clean bulk rate.
%     F2 kink visualization: along the horizontal line y ~ 0 crossing the
%        pentagon's right edge, du/dx jumps from the PDE-determined
%        interior slope to ~0 outside (u ~ g = 0 there), smoothed over an
%        O(eps) layer that sharpens as eps -> 0 (10-90 width ~ 0.73 eps).
%     F3 contrast: the Neumann flux form with g = 0.3 keeps du/dn
%        CONTINUOUS across Gamma at leading order (the jump moves to the
%        second derivative), as predicted by the matched asymptotics.
%
%   Writes ../output/results_forced.mat and figures.
clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', 'src', 'ddAddpaths.m'));
outdir = fullfile(here, '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end
warning('off', 'eb_tag_cluster:dropthin');
R = struct();
PASS = true;
t0all = tic;

box = [-2.5 2.5 -2.5 2.5];
geom = ddGeometry('pentagon', 1.2);
ffun = @(x,y) 1 + 0.5*x;                  % prescribed forcing (no MMS)
dataD = struct('f', ffun, 'g', 0, 'C', 0);    % Dirichlet penalty, u=0 on Gamma
N0 = 256;                                  % common projection grid
[X0, Y0] = ndgrid(box(1) + ((1:N0)-0.5)*(box(2)-box(1))/N0, ...
                  box(3) + ((1:N0)-0.5)*(box(4)-box(3))/N0);
insideM = geom.psi(X0, Y0) < 0;
da0 = ((box(2)-box(1))/N0)^2;

% ======================================================================
fprintf('\n=== F1a: h-Cauchy at fixed eps = 0.2 (Dirichlet penalty) =========\n');
% Bulk (eroded) and full-interior differences are reported separately:
% the layer's L2 discretization constant scales like 1/eps^2, so the
% full-interior norm reaches its asymptotic rate only once h << eps;
% the bulk shows the clean order-2 immediately.
ep = 0.2;
Nh = [256 512 1024 2048];
bulkM = geom.psi(X0, Y0) < -2.5 * ep;        % inside, beyond the layer
Uh = cell(1, numel(Nh));
for i = 1:numel(Nh)
  sol = ddSolve(geom, ep, box, 'dirichlet', dataD, ...
    'BaseCells', (Nh(i)/8)*[1 1], 'MaxLevels', 4, 'TagWidth', 2.5, ...
    'TagBuffer', 2, 'Mu', 1, 'Verbose', 0, 'RelTol', 1e-9, 'MaxVCycles', 40);
  Uh{i} = ddToUniform(sol.H, N0);
  st = sol.stats.stages{end};
  fprintf('  N=%4d: dof=%7d  factor=%.3f  relres=%.1e  mass=%.0e\n', ...
    Nh(i), st.visibleCells, st.avgFactor, st.finalRel, sol.mass.rel);
  % the conservation defect tracks the solve tolerance (RelTol 1e-9)
  % scaled by the eps^-3 penalty reaction; tightening RelTol drives it
  % to roundoff as in the main suite
  if sol.mass.rel > 2e-9, fprintf('  ** F1a FAIL (mass)\n'); PASS = false; end
end
dh = nan(1, numel(Nh)-1); dhb = dh;
for i = 1:numel(Nh)-1
  D = Uh{i} - Uh{i+1};
  dh(i)  = sqrt(sum(D(insideM).^2) * da0);
  dhb(i) = sqrt(sum(D(bulkM).^2) * da0);
end
rh = log2(dh(1:end-1) ./ dh(2:end));
rhb = log2(dhb(1:end-1) ./ dhb(2:end));
fprintf('  ||u_N - u_2N||_L2: inside %s | bulk %s\n', mat2str(dh, 3), mat2str(dhb, 3));
fprintf('  h-rates: inside %s | bulk %s\n', mat2str(rh, 2), mat2str(rhb, 2));
R.F1a = struct('N', Nh, 'd', dh, 'dBulk', dhb, 'rates', rh, 'ratesBulk', rhb);
if rhb(end) < 1.7
  fprintf('  ** F1a FAIL (bulk h-rate < 1.7)\n'); PASS = false;
end
if rh(end) < 1.0
  fprintf('  ** F1a FAIL (full-interior h-rate < 1.0)\n'); PASS = false;
end

% ======================================================================
fprintf('\n=== F1b: eps-Cauchy at fixed h (effective 1024) ===================\n');
epsl = [0.4 0.2 0.1 0.05];
Ue = cell(1, numel(epsl));
for j = 1:numel(epsl)
  sol = ddSolve(geom, epsl(j), box, 'dirichlet', dataD, ...
    'BaseCells', 128*[1 1], 'MaxLevels', 4, 'TagWidth', 2.5, ...
    'TagBuffer', 2, 'Mu', 1, 'Verbose', 0, 'RelTol', 1e-9, 'MaxVCycles', 40);
  Ue{j} = ddToUniform(sol.H, N0);
end
% bulk away from all layers; the pentagon's inradius is only 0.97, so
% erode by 2*eps_max = 0.8 (tanh saturation there: 1-phi ~ 6e-6)
erodedM = geom.psi(X0, Y0) < -2 * max(epsl);
de = nan(1, numel(epsl)-1); dee = de;
for j = 1:numel(epsl)-1
  D = Ue{j} - Ue{j+1};
  de(j)  = sqrt(sum(D(insideM).^2) * da0);
  dee(j) = sqrt(sum(D(erodedM).^2) * da0);
end
re = log2(de(1:end-1) ./ de(2:end));
ree = log2(dee(1:end-1) ./ dee(2:end));
fprintf('  ||u_eps - u_eps/2||: inside %s | eroded bulk %s\n', ...
  mat2str(de, 3), mat2str(dee, 3));
fprintf('  eps-rates: inside %s | eroded bulk %s\n', mat2str(re, 2), mat2str(ree, 2));
R.F1b = struct('eps', epsl, 'd', de, 'dErode', dee, 'rates', re, 'ratesErode', ree);
if any(re < 1.2) || any(ree < 1.4)
  fprintf('  ** F1b FAIL (eps self-convergence too slow)\n'); PASS = false;
end

% ======================================================================
fprintf('\n=== F2: normal-derivative jump at Gamma (Dirichlet penalty) =======\n');
% uniform runs for clean line sampling; right edge crossing at y ~ 0:
% edge normal angle -18 deg, line x_G = apothem / cos(18 deg)
xG = 1.2 * cos(pi/5) / cos(pi/10);
Nu = 1024;                    % keep eps/h >= 5 down to eps = 0.025
h = (box(2)-box(1)) / Nu;
xs = box(1) + ((1:Nu) - 0.5) * h;
j0 = round((0 - box(3)) / h + 0.5);
epk = [0.2 0.1 0.05];         % thin vs the pentagon inradius 0.97
prof = cell(1, numel(epk)); dprof = prof;
width = nan(1, numel(epk)); slopeIn = width;
for j = 1:numel(epk)
  sol = ddSolve(geom, epk(j), box, 'dirichlet', dataD, ...
    'BaseCells', Nu*[1 1], 'MaxLevels', 1, 'Verbose', 0, ...
    'RelTol', 1e-9, 'MaxVCycles', 40);
  Q = sol.H.lev{end}.Q{1};
  u = Q(2:end-1, j0+1).';
  du = nan(size(u));
  du(2:end-1) = (u(3:end) - u(1:end-2)) / (2*h);
  prof{j} = u; dprof{j} = du;
  % Max-slope (shock-thickness) width of the du/dx transition: robust to
  % the smooth bulk variation of du/dx.  slopeIn = the most negative
  % du/dx in the window (the wall shoulder), Lout = saturated exterior
  % value; width = (Lout - slopeIn) / max d(du)/dx in the layer.
  win = find(xs > xG - 3.5*epk(j) & xs < xG + 3.5*epk(j));
  dw = du(win);
  Lout = interp1(xs, du, xG + 2.5*epk(j));
  slopeIn(j) = min(dw);
  ddu = (dw(3:end) - dw(1:end-2)) / (2*h);
  width(j) = (Lout - slopeIn(j)) / max(ddu);
  fprintf('  eps=%5.3f: du/dx wall shoulder=%.3f outside=%.3f  max-slope width=%.4f (width/eps=%.2f)\n', ...
    epk(j), slopeIn(j), Lout, width(j), width(j)/epk(j));
end
fprintf('  width ratios w(eps/2)/w(eps): %s  (layer thickness ~ eps)\n', ...
  mat2str(width(2:end)./width(1:end-1), 2));
fprintf('  wall-shoulder slope converging to the sharp normal derivative: %s\n', ...
  mat2str(slopeIn, 3));
R.F2 = struct('eps', epk, 'width', width, 'slopeIn', slopeIn, 'xG', xG);
if ~all(isfinite(width)) || any(diff(width) >= 0) ...
    || width(end)/width(1) > 0.55 || any(width > 1.5*epk)
  fprintf('  ** F2 FAIL (kink not sharpening like eps)\n'); PASS = false;
end
if abs(slopeIn(end) - slopeIn(end-1)) > abs(slopeIn(end-1) - slopeIn(1))
  fprintf('  note: wall-shoulder slope not yet converging\n');
end

% ======================================================================
fprintf('\n=== F3: Neumann flux form: derivative CONTINUOUS across Gamma =====\n');
dataN = struct('f', ffun, 'g', 0.3, 'C', 1);
profN = cell(1, numel(epk)); dprofN = profN;
for j = 1:numel(epk)
  sol = ddSolve(geom, epk(j), box, 'neumann', dataN, ...
    'BaseCells', Nu*[1 1], 'MaxLevels', 1, 'Verbose', 0, ...
    'RelTol', 1e-9, 'MaxVCycles', 40);
  Q = sol.H.lev{end}.Q{1};
  u = Q(2:end-1, j0+1).';
  du = nan(size(u));
  du(2:end-1) = (u(3:end) - u(1:end-2)) / (2*h);
  profN{j} = u; dprofN{j} = du;
  jumpN = abs(interp1(xs, du, xG + 2*epk(j)) - interp1(xs, du, xG - 2*epk(j)));
  fprintf('  eps=%4.2f: |du/dx(xG+2eps) - du/dx(xG-2eps)| = %.4f (Dirichlet: %.3f)\n', ...
    epk(j), jumpN, abs(interp1(xs, dprof{j}, xG + 2*epk(j)) ...
                     - interp1(xs, dprof{j}, xG - 2*epk(j))));
end
R.F3prof = {xs, prof, dprof, profN, dprofN};

% ======================================================================
% figures
co = lines(numel(epk));
fig = figure('Visible', 'off', 'Position', [80 80 1000 420]);
subplot(1,2,1); hold on; grid on;
for j = 1:numel(epk)
  plot(xs, prof{j}, '-', 'Color', co(j,:), 'LineWidth', 1.4);
end
xline(xG, 'k--');
xlim([0 2]); xlabel('x (y = 0)'); ylabel('u_\epsilon');
title('Dirichlet penalty: solution profile');
legend(arrayfun(@(e) sprintf('\\epsilon = %.2f', e), epk, 'UniformOutput', false), ...
  'Location', 'southwest');
subplot(1,2,2); hold on; grid on;
for j = 1:numel(epk)
  plot(xs, dprof{j}, '-', 'Color', co(j,:), 'LineWidth', 1.4);
end
xline(xG, 'k--');
xlim([0.6 1.5]); xlabel('x (y = 0)'); ylabel('du_\epsilon/dx');
title('du/dx: O(1) jump at \Gamma, smoothed over O(\epsilon)');
exportgraphics(fig, fullfile(outdir, 'f2_kink_dirichlet.png'), 'Resolution', 150);
close(fig);

fig = figure('Visible', 'off', 'Position', [80 80 1000 420]);
subplot(1,2,1); hold on; grid on;
for j = 1:numel(epk)
  plot(xs, dprofN{j}, '-', 'Color', co(j,:), 'LineWidth', 1.4);
end
xline(xG, 'k--');
xlim([0.6 1.5]); xlabel('x (y = 0)'); ylabel('du_\epsilon/dx');
title('Neumann form: du/dx continuous across \Gamma');
legend(arrayfun(@(e) sprintf('\\epsilon = %.2f', e), epk, 'UniformOutput', false), ...
  'Location', 'northeast');
subplot(1,2,2);
loglog(Nh(1:end-1), dh, 'o-', 'LineWidth', 1.4); hold on; grid on;
loglog(Nh(1:end-1), dhb, 's-', 'LineWidth', 1.4);
loglog(Nh(1:end-1), dhb(1)*(Nh(1:end-1)/Nh(1)).^(-2), 'k--', 'LineWidth', 1);
xlabel('N (effective)'); ylabel('||u_N - u_{2N}||_{L^2}');
legend({'inside \Omega_1 (incl. layer)', 'bulk (\psi < -2.5\epsilon)', ...
  'slope -2'}, 'Location', 'southwest');
title(sprintf('h-Cauchy self-convergence, \\epsilon = %.2f', ep));
exportgraphics(fig, fullfile(outdir, 'f3_neumann_and_cauchy.png'), 'Resolution', 150);
close(fig);

% ======================================================================
R.elapsed = toc(t0all);
save(fullfile(outdir, 'results_forced.mat'), 'R');
fprintf('\nTotal forced-study time: %.1f s\n', R.elapsed);
if PASS
  fprintf('ALL FORCED-PROBLEM TESTS PASSED\n');
else
  error('run_dd_forced: at least one test failed');
end

