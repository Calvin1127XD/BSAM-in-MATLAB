%stressTreeInterfaces Adversarial interface-geometry stress suite for
%   EllipticBSAMTree (strict single-parent tree solver).  Mirrors
%   EllipticBSAM/tests/stressInterfaces.m case-for-case and adds the
%   tree-specific verification:
%
%     * ebtTreeCheck is run on EVERY hierarchy (single-parent containment,
%       box/mbounds consistency, alignment, disjointness, children/covered
%       bookkeeping, ring nesting);
%     * the coefficient fill-down duplicate probe: coarse faces shared by
%       abutting patches are stored twice; after ebtCoeffFilldown both
%       copies must be IDENTICAL (the pushFaceX/pushFaceY sibling routing
%       -- the historical T8 conservation-bug locus, only visible with
%       variable D);
%     * cross-solver check: on hand-built layouts both solvers build the
%       IDENTICAL mesh, and every operation is the same arithmetic on the
%       same numbers, so the solutions must agree to machine precision.
%
%   Case families (see the ring suite header for details):
%     S1 combs (interior / boundary / Neumann-side / anisotropic, QIF+LIF)
%     S2 corner-only adjacency        S3 reentrant tag shapes + annulus
%     S4 boundary-hanging interfaces  S5 minimum-width patches
%     S6 adversarial layouts (checkerboards, spiral, corridors, EDGEKISS)
%     S7 deep 4-patch cross point (5 levels)
%     S8 conservation discriminators (defect must scale with RelTol)
%     S9 singular pure-Neumann on exotic meshes (SourceQuad=3)
%     S10 ring-vs-tree bit-agreement on identical hand-built meshes
%
%   Prints PASS/FAIL per case and errors out at the end on any failure.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
addpath(fullfile(here, '..', '..', 'EllipticBSAM', 'src'));   % for S10 only
outdir = fullfile(here, '..', 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end
warning('off', 'ebt_tag_cluster:dropthin');
warning('off', 'eb_tag_cluster:dropthin');
PASS = true;
say = @(varargin) fprintf(varargin{:});
R = struct();
t0 = tic;

% ---------------------------------------------------------------- problems
pv.uex = @(x,y) sin(pi*x).*cos(pi*y);
pv.D   = @(x,y) 1 + 0.5*sin(2*pi*x).*cos(pi*y);
pv.C   = @(x,y) 1 + x.*y;
pvDx = @(x,y) pi*cos(2*pi*x).*cos(pi*y);
pvDy = @(x,y) -0.5*pi*sin(2*pi*x).*sin(pi*y);
pvux = @(x,y) pi*cos(pi*x).*cos(pi*y);
pvuy = @(x,y) -pi*sin(pi*x).*sin(pi*y);
pvf = @(x,y) -(pvDx(x,y).*pvux(x,y) + pvDy(x,y).*pvuy(x,y) ...
             + pv.D(x,y).*(-2*pi^2*pv.uex(x,y))) + pv.C(x,y).*pv.uex(x,y);
probV = ebtProblem([0 1 0 1], pv.D, pv.C, pvf, {'dirichlet', pv.uex}, pv.uex);
probVring = ebProblem([0 1 0 1], pv.D, pv.C, pvf, {'dirichlet', pv.uex}, pv.uex);

pm.uex = @(x,y) exp(x).*sin(pi*y) + 0.5*x.^2;
pm.D   = @(x,y) 1 + 0.5*x.*y;
pm.C   = @(x,y) 1 + x + y;
pmf = @(x,y) -( 0.5*y .* (exp(x).*sin(pi*y) + x) ...
              + (1 + 0.5*x.*y) .* (exp(x).*sin(pi*y) + 1) ...
              + 0.5*x .* (pi*exp(x).*cos(pi*y)) ...
              + (1 + 0.5*x.*y) .* (-pi^2*exp(x).*sin(pi*y)) ) ...
            + pm.C(x,y) .* pm.uex(x,y);
pmgW = @(x,y) -(exp(x).*sin(pi*y) + x);
pmgE = @(x,y)  (exp(x).*sin(pi*y) + x);
pmgS = @(x,y) -(pi*exp(x).*cos(pi*y));
probMD = ebtProblem([0 1 0 1], pm.D, pm.C, pmf, {'dirichlet', pm.uex}, pm.uex);
bcNN = struct('west', {{'neumann', pmgW}}, 'south', {{'neumann', pmgS}}, ...
              'east', {{'dirichlet', pm.uex}}, 'north', {{'dirichlet', pm.uex}});
probMN = ebtProblem([0 1 0 1], pm.D, pm.C, pmf, bcNN, pm.uex);
bcMX = struct('west', {{'dirichlet', pm.uex}}, 'south', {{'neumann', pmgS}}, ...
              'east', {{'neumann', pmgE}}, 'north', {{'dirichlet', pm.uex}});
probMX = ebtProblem([0 1 0 1], pm.D, pm.C, pmf, bcMX, pm.uex);

pa.uex = @(x,y) sin(pi*x).*cos(pi*y);
paf = @(x,y) -0.5*pi*cos(pi*x).*cos(pi*y) ...
           + 2*pi^2*(1 + 0.5*x).*sin(pi*x).*cos(pi*y) + pa.uex(x,y);
probA = ebtProblem([0 2 0 1], @(x,y) 1 + 0.5*x, 1, paf, ...
                   {'dirichlet', pa.uex}, pa.uex);

pnuex = @(x,y) cos(2*pi*x).*cos(2*pi*y);
pnux = @(x,y) -2*pi*sin(2*pi*x).*cos(2*pi*y);
pnuy = @(x,y) -2*pi*cos(2*pi*x).*sin(2*pi*y);
pnf = @(x,y) -(pvDx(x,y).*pnux(x,y) + pvDy(x,y).*pnuy(x,y) ...
             + pv.D(x,y).*(-8*pi^2*pnuex(x,y)));
probN = ebtProblem([0 1 0 1], pv.D, 0, pnf, {'neumann', 0}, pnuex);

ms = [32 64 128];
rateBar = 1.75; facBar = 0.2; massBar = 5e-10;

% ---------------------------------------------------------------- shapes
combfn = @(X,Y) (X>4/32 & X<16/32 & Y>4/32 & Y<28/32) | ...
  (X>16/32 & X<20/32 & ( (Y>4/32 & Y<7/32) | (Y>10/32 & Y<13/32) | ...
                         (Y>16/32 & Y<19/32) | (Y>22/32 & Y<25/32) ));
combb = @(X,Y) (X < 4/32) & (32*Y > 2) & (32*Y < 28) & ...
               (mod(floor((32*Y - 2)/2), 2) == 0);
combA = @(X,Y) (X > 0.5 & X < 1.0 & Y > 1/8 & Y < 7/8) | ...
  (X >= 1.0 & X < 1.25 & ( (Y > 1/8 & Y < 2/8) | (Y > 3/8 & Y < 4/8) | ...
                           (Y > 5/8 & Y < 6/8) ));
u16 = @(X,Y,x1,x2,y1,y2) (X>=x1/16 & X<=x2/16 & Y>=y1/16 & Y<=y2/16);
zig = @(X,Y) (floor(8*X - 8*Y) >= -1 & floor(8*X - 8*Y) <= 0 & ...
              X > 2/16 & X < 14/16 & Y > 2/16 & Y < 14/16);
ann = @(X,Y) max(abs(X-0.5), abs(Y-0.5)) >= 3/16 & ...
             max(abs(X-0.5), abs(Y-0.5)) <= 6/16;

% ======================================================================
say('\n=== S1: comb interfaces (tree) ==================================\n');
[PASS, R.S1a] = tagSweepT(PASS, 'COMB qif', probV, ms, @(X,Y,L) combfn(X,Y), ...
  0.95, 'qif', rateBar, facBar, massBar, 4);
[PASS, R.S1b] = tagSweepT(PASS, 'COMB lif', probV, ms, @(X,Y,L) combfn(X,Y), ...
  0.95, 'lif', rateBar, facBar, massBar, 4);
[PASS, R.S1c] = tagSweepT(PASS, 'COMBB (bdy)', probV, ms, @(X,Y,L) combb(X,Y), ...
  0.85, 'qif', rateBar, facBar, massBar, 5);
[PASS, R.S1d] = tagSweepT(PASS, 'COMBB-NN', probMN, ms, @(X,Y,L) combb(X,Y), ...
  0.85, 'qif', rateBar, facBar, massBar, 5);
errs = nan(1,3); facs = errs; mass = errs;
bases = {[64 16], [128 32], [256 64]};
for i = 1:3
  o = ebtOptions('BaseCells', bases{i}, 'MaxLevels', 2, 'Indicator', 'fun', ...
    'TagFunction', @(X,Y,L) combA(X,Y), 'TagBuffer', 0, 'MinBoxWidth', 2, ...
    'MinFillRatio', 0.95, 'Verbose', 0, 'RelTol', 1e-10, 'MaxVCycles', 60);
  sol = ebtSolve(probA, o);
  ebtTreeCheck(sol.H);
  errs(i) = sol.err.linf; facs(i) = sol.stats.stages{end}.avgFactor;
  mass(i) = sol.mass.rel;
end
r = log2(errs(1:2)./errs(2:3));
say('  %-12s errs=[%.3e %.3e %.3e] rates=[%.2f %.2f] fac<=%.3f mass<=%.1e\n', ...
  'COMB-ANISO', errs, r, max(facs), max(mass));
if any(r < rateBar) || max(facs) > 0.25 || max(mass) > massBar
  say('  ** S1 COMB-ANISO FAIL\n'); PASS = false;
end
R.S1e = struct('errs', errs, 'rates', r, 'facs', facs, 'mass', mass);

% ======================================================================
say('\n=== S2: corner-only adjacency (tree) ============================\n');
S2 = {};
S2{end+1} = {'CORNER2',   {[9 16 9 16; 17 24 17 24]}};
S2{end+1} = {'CORNER2AD', {[17 24 9 16; 9 16 17 24]}};
S2{end+1} = {'CORNERA',   {[5 16 5 16; 17 28 17 24]}};
S2{end+1} = {'CORNER4',   {[9 16 9 16; 17 24 9 16; 9 16 17 24; 17 24 17 24]}};
for c = 1:numel(S2)
  [PASS, R.(sprintf('S2_%d', c))] = boxSweepT(PASS, S2{c}{1}, probV, ms, ...
    S2{c}{2}, rateBar, facBar, massBar);
end

% ======================================================================
say('\n=== S3: reentrant tag shapes (tree) =============================\n');
S3 = {};
S3{end+1} = {'USHAPE', @(X,Y,L) u16(X,Y,4,12,4,6)|u16(X,Y,4,6,6,12)|u16(X,Y,10,12,6,12), 0.95};
S3{end+1} = {'TSHAPE', @(X,Y,L) u16(X,Y,4,12,8,10)|u16(X,Y,7,9,2,8), 0.95};
S3{end+1} = {'PLUS',   @(X,Y,L) u16(X,Y,3,13,7,9)|u16(X,Y,7,9,3,13), 0.95};
S3{end+1} = {'ZIGZAG', @(X,Y,L) zig(X,Y), 0.85};
S3{end+1} = {'ANNULUS',@(X,Y,L) ann(X,Y), 0.85};
for c = 1:numel(S3)
  [PASS, out] = tagSweepT(PASS, S3{c}{1}, probV, ms, S3{c}{2}, S3{c}{3}, ...
    'qif', rateBar, facBar, massBar, 2);
  R.(sprintf('S3_%d', c)) = out;
  if strcmp(S3{c}{1}, 'ANNULUS') && any(out.covCtr)
    say('  ** S3 ANNULUS FAIL (hole covered)\n'); PASS = false;
  end
end

% ======================================================================
say('\n=== S4: hanging nodes at the physical boundary (tree) ===========\n');
S4 = {};
S4{end+1} = {'HANGQ-D',   probMD, {[1 16 1 16]}};
S4{end+1} = {'HANGQ-NN',  probMN, {[1 16 1 16]}};
S4{end+1} = {'HANGQ-MX',  probMX, {[1 16 1 16]}};
S4{end+1} = {'NOTCH-MX',  probMX, {[9 24 1 8]}};
S4{end+1} = {'STRIP-MX',  probMX, {[1 16 1 32]}};
S4{end+1} = {'HANGQ3-MX', probMX, {[1 16 1 16], [1 16 1 16]}};
for c = 1:numel(S4)
  [PASS, R.(sprintf('S4_%d', c))] = boxSweepT(PASS, S4{c}{1}, S4{c}{2}, ms, ...
    S4{c}{3}, rateBar, facBar, massBar);
end

% ======================================================================
say('\n=== S5: minimum-width patches (tree) ============================\n');
chain = zeros(7,4); for i = 1:7, a = 7+2*i; chain(i,:) = [a a+1 a a+1]; end
S5 = {};
S5{end+1} = {'SANDWICH',  {[5 28 15 16; 5 28 17 26; 5 28 27 28]}};
S5{end+1} = {'DOMINOES',  {[9 10 9 24; 11 12 9 24; 13 14 9 24; 15 16 9 24]}};
S5{end+1} = {'CHAIN',     {chain}};
S5{end+1} = {'RINGTOUCH', {[1 4 13 16; 1 4 9 12; 1 4 17 20; 5 8 13 16]}};
S5{end+1} = {'SKINNY3',   {[9 24 9 24], [25 26 19 42]}};
for c = 1:numel(S5)
  [PASS, R.(sprintf('S5_%d', c))] = boxSweepT(PASS, S5{c}{1}, probV, ms, ...
    S5{c}{2}, rateBar, facBar, massBar);
end

% ======================================================================
say('\n=== S6: adversarial hand-built layouts (tree) ===================\n');
[ii, jj] = ndgrid(0:7, 0:7); keep = mod(ii+jj,2)==0;
cb = [4*ii(keep)+1, 4*ii(keep)+4, 4*jj(keep)+1, 4*jj(keep)+4];
[ii, jj] = ndgrid(0:9, 0:9); keep = mod(ii+jj,2)==0;
cb3 = [4*ii(keep)+13, 4*ii(keep)+16, 4*jj(keep)+13, 4*jj(keep)+16];
spiral = [1 32 1 4; 29 32 5 32; 1 28 29 32; 1 4 5 28; 5 28 5 8; ...
          25 28 9 28; 5 24 25 28; 5 8 9 24; 9 24 9 12; 21 24 13 24; ...
          9 20 21 24; 9 12 13 20; 13 20 13 16];
corridors = [5 26 5 6; 27 28 5 26; 7 28 27 28; 5 6 7 28; ...
             9 22 9 10; 23 24 9 22; 11 24 23 24; 9 10 11 22];
S6 = {};
S6{end+1} = {'CHECKER',   {cb}};
S6{end+1} = {'CHECKER3',  {[5 28 5 28], cb3}};
S6{end+1} = {'SPIRAL',    {spiral}};
S6{end+1} = {'CORRIDORS', {corridors}};
S6{end+1} = {'EDGEKISS',  {[5 16 5 28; 17 28 5 28], [25 32 17 40]}};
S6{end+1} = {'EDGEKISS2', {[5 16 5 28; 17 28 5 28], [25 32 17 40; 33 40 29 52]}};
S6{end+1} = {'SURROUND',  {[9 12 9 12; 13 20 9 12; 21 24 9 12; 9 12 13 20; ...
  13 20 13 20; 21 24 13 20; 9 12 21 24; 13 20 21 24; 21 24 21 24]}};
for c = 1:numel(S6)
  [PASS, R.(sprintf('S6_%d', c))] = boxSweepT(PASS, S6{c}{1}, probV, ms, ...
    S6{c}{2}, rateBar, facBar, massBar);
end
% shared-face coefficient duplicate probe (pushFaceX/pushFaceY routing):
% run on the layouts whose children touch parent edges + the checkerboards
for c = [1 2 5 6]
  H = buildBoxesT(probV, 32, S6{c}{2}, 1e-9, 2);
  dup = coeffDupMismatch(H);
  say('  %-10s shared-face coeff duplicate mismatch = %.2e\n', S6{c}{1}, dup);
  if dup > 1e-14
    say('  ** S6 %s FAIL (fill-down left inconsistent face copies)\n', S6{c}{1});
    PASS = false;
  end
end

% ======================================================================
say('\n=== S7: deep 4-patch cross point (tree, 5 levels) ================\n');
dc = cell(1,4);
for lv = 2:5
  n = 32*2^(lv-2); c0 = n/2;
  dc{lv-1} = [c0-7 c0 c0-7 c0; c0+1 c0+8 c0-7 c0; ...
              c0-7 c0 c0+1 c0+8; c0+1 c0+8 c0+1 c0+8];
end
for depth = 3:5
  H = buildBoxesT(probV, 32, dc(1:depth-1), 1e-10, 60);
  st = H.stats.stages{end};
  e = ebtErrorNorms(H); md = ebtMassCheck(H);
  say('  depth=%d: V=%2d fac=%.3f rel=%.1e err=%.3e mass=%.1e\n', ...
    depth, st.vcycles, st.avgFactor, st.finalRel, e.linf, md.rel);
  if st.avgFactor > facBar || md.rel > massBar
    say('  ** S7 FAIL (depth %d degraded)\n', depth); PASS = false;
  end
end
[PASS, R.S7] = boxSweepT(PASS, 'DEEPCROSS5', probV, ms, dc, rateBar, facBar, massBar);

% ======================================================================
say('\n=== S8: conservation discriminators (tree) ======================\n');
disc = {};
disc{end+1} = {'COMB',      [], @(X,Y,L) combfn(X,Y), 0.95};
disc{end+1} = {'ANNULUS',   [], @(X,Y,L) ann(X,Y), 0.85};
disc{end+1} = {'EDGEKISS2', S6{6}{2}, [], []};
disc{end+1} = {'RINGTOUCH', S5{4}{2}, [], []};
disc{end+1} = {'CHECKER3',  S6{2}{2}, [], []};
disc{end+1} = {'DEEPCROSS5',dc,       [], []};
for c = 1:numel(disc)
  md = nan(1,2); tols = [1e-10 1e-12];
  for i = 1:2
    if isempty(disc{c}{2})
      o = ebtOptions('BaseCells', [64 64], 'MaxLevels', 2, 'Indicator', 'fun', ...
        'TagFunction', disc{c}{3}, 'TagBuffer', 0, 'MinBoxWidth', 2, ...
        'MinFillRatio', disc{c}{4}, 'Verbose', 0, 'RelTol', tols(i), ...
        'MaxVCycles', 90);
      sol = ebtSolve(probV, o);
      md(i) = sol.mass.rel;
    else
      lv = cellfun(@(b) scaleBoxes(b, 2), disc{c}{2}, 'UniformOutput', false);
      H = buildBoxesT(probV, 64, lv, tols(i), 90);
      mm = ebtMassCheck(H);
      md(i) = mm.rel;
    end
  end
  say('  %-10s defect %.2e (1e-10) -> %.2e (1e-12)\n', disc{c}{1}, md(1), md(2));
  if md(1) > massBar || md(2) > 1e-12
    say('  ** S8 FAIL (%s: tolerance-independent defect floor)\n', disc{c}{1});
    PASS = false;
  end
end
mdOn = nan; mdOff = nan;
for mc = [true false]
  o = ebtOptions('BaseCells', [64 64], 'MaxLevels', 2, 'Indicator', 'fun', ...
    'TagFunction', @(X,Y,L) combfn(X,Y), 'TagBuffer', 0, 'MinBoxWidth', 2, ...
    'MinFillRatio', 0.95, 'MassCorrection', mc, 'Verbose', 0, ...
    'RelTol', 1e-11, 'MaxVCycles', 90);
  sol = ebtSolve(probV, o);
  if mc, mdOn = sol.mass.rel; else, mdOff = sol.mass.rel; end
end
say('  COMB corrections on/off: %.2e / %.2e\n', mdOn, mdOff);
if ~(mdOff > 1e3 * max(mdOn, 1e-16))
  say('  ** S8 FAIL (corrections-off defect not dramatically larger)\n');
  PASS = false;
end
R.S8 = struct('on', mdOn, 'off', mdOff);

% ======================================================================
say('\n=== S9: singular pure-Neumann on exotic meshes (tree) ===========\n');
errs = nan(1,3); facs = errs; mass = errs;
for i = 1:3
  o = ebtOptions('BaseCells', ms(i)*[1 1], 'MaxLevels', 2, 'Indicator', 'fun', ...
    'TagFunction', @(X,Y,L) ann(X,Y), 'TagBuffer', 0, 'MinBoxWidth', 2, ...
    'MinFillRatio', 0.85, 'SourceQuad', 3, 'Verbose', 0, ...
    'RelTol', 1e-10, 'MaxVCycles', 60);
  sol = ebtSolve(probN, o);
  ebtTreeCheck(sol.H);
  errs(i) = sol.err.linf; facs(i) = sol.stats.stages{end}.avgFactor;
  mass(i) = sol.mass.rel;
end
r = log2(errs(1:2)./errs(2:3));
say('  %-12s errs=[%.3e %.3e %.3e] rates=[%.2f %.2f] fac<=%.3f mass<=%.1e\n', ...
  'ANNULUS-N q3', errs, r, max(facs), max(mass));
if any(r < rateBar) || max(facs) > facBar || max(mass) > 1e-10
  say('  ** S9 FAIL (annulus, singular Neumann, SourceQuad=3)\n'); PASS = false;
end
R.S9 = struct('errs', errs, 'rates', r, 'facs', facs, 'mass', mass);

% ======================================================================
say('\n=== S10: ring-vs-tree bit-agreement on identical meshes ==========\n');
% Identical hand-built boxes -> identical meshes -> identical arithmetic:
% the two implementations must agree to machine precision.  (Interior
% differences would expose divergent window/ring/correction routing.)
XC = {};
XC{end+1} = {'CORNER4',   S2{4}{2}, [32 64]};
XC{end+1} = {'CHAIN',     S5{3}{2}, [32]};
XC{end+1} = {'RINGTOUCH', S5{4}{2}, [32]};
XC{end+1} = {'SKINNY3',   S5{5}{2}, [32]};
XC{end+1} = {'CHECKER',   S6{1}{2}, [32 64]};
XC{end+1} = {'CHECKER3',  S6{2}{2}, [32]};
XC{end+1} = {'SPIRAL',    S6{3}{2}, [32]};
XC{end+1} = {'CORRIDORS', S6{4}{2}, [32]};
XC{end+1} = {'EDGEKISS2', S6{6}{2}, [32 64]};
XC{end+1} = {'DEEPCROSS5', dc,      [32 64]};
XC{end+1} = {'HANGQ3',    {[1 16 1 16], [1 16 1 16]}, [32]};
for c = 1:numel(XC)
  for m = XC{c}{3}
    lv = cellfun(@(b) scaleBoxes(b, m/32), XC{c}{2}, 'UniformOutput', false);
    Ht = buildBoxesT(probV, m, lv, 1e-10, 60);
    Hr = buildBoxesRing(probVring, m, lv, 1e-10, 60);
    d = maxQdiff(Hr, Ht);
    vT = Ht.stats.stages{end}.vcycles; vR = Hr.stats.stages{end}.vcycles;
    say('  %-10s m=%3d  |ring-tree|=%.2e  V(ring/tree)=%d/%d\n', ...
      XC{c}{1}, m, d, vR, vT);
    if d > 1e-11 || vR ~= vT
      say('  ** S10 FAIL (%s m=%d: solvers diverge on the identical mesh)\n', ...
        XC{c}{1}, m);
      PASS = false;
    end
  end
end

% ======================================================================
R.elapsed = toc(t0);
save(fullfile(outdir, 'stressTreeResults.mat'), 'R');
say('\nTotal tree stress-suite time: %.1f s\n', R.elapsed);
if PASS
  say('ALL TREE STRESS-INTERFACE TESTS PASSED\n');
else
  error('stressTreeInterfaces: at least one test failed');
end

% ========================================================================
% helpers
% ========================================================================
function [PASS, out] = tagSweepT(PASS, name, prob, ms, tagfn, fill, ghost, ...
                                 rateBar, facBar, massBar, npMin)
%tagSweepT MMS sweep of a tag-driven 2-level tree mesh (+ ebtTreeCheck).
errs = nan(1,numel(ms)); facs = errs; mass = errs; nps = errs; covc = errs;
for i = 1:numel(ms)
  o = ebtOptions('BaseCells', ms(i)*[1 1], 'MaxLevels', 2, 'Indicator', 'fun', ...
    'TagFunction', tagfn, 'TagBuffer', 0, 'MinBoxWidth', 2, ...
    'MinFillRatio', fill, 'GhostInterp', ghost, 'Verbose', 0, ...
    'RelTol', 1e-10, 'MaxVCycles', 60);
  sol = ebtSolve(prob, o);
  ebtTreeCheck(sol.H);
  errs(i) = sol.err.linf; facs(i) = sol.stats.stages{end}.avgFactor;
  mass(i) = sol.mass.rel; nps(i) = sol.levels(2).np;
  lev = sol.H.lev{sol.H.nlev};
  ctr = false;
  for p = 1:lev.np
    b = lev.box(p, :);
    if b(1) <= lev.n(1)/2 && b(2) >= lev.n(1)/2 && ...
       b(3) <= lev.n(2)/2 && b(4) >= lev.n(2)/2
      ctr = true;
    end
  end
  covc(i) = ctr;
end
r = log2(errs(1:end-1)./errs(2:end));
fprintf('  %-12s np=[%s] errs=[%.3e %.3e %.3e] rates=[%.2f %.2f] fac<=%.3f mass<=%.1e\n', ...
  name, strtrim(sprintf('%d ', nps)), errs, r, max(facs), max(mass));
ok = all(r > rateBar) && max(facs) <= facBar && max(mass) <= massBar ...
     && all(nps >= npMin);
if ~ok
  fprintf('  ** %s FAIL (rates/factor/mass/patch-count)\n', name);
  PASS = false;
end
out = struct('errs', errs, 'rates', r, 'facs', facs, 'mass', mass, ...
             'np', nps, 'covCtr', covc);
end

% ------------------------------------------------------------------------
function [PASS, out] = boxSweepT(PASS, name, prob, ms, levelsBase, ...
                                 rateBar, facBar, massBar)
%boxSweepT MMS sweep of a hand-built tree hierarchy (+ ebtTreeCheck).
errs = nan(1,numel(ms)); facs = errs; mass = errs;
for i = 1:numel(ms)
  s = ms(i)/32;
  lv = cellfun(@(b) scaleBoxes(b, s), levelsBase, 'UniformOutput', false);
  H = buildBoxesT(prob, ms(i), lv, 1e-10, 60);
  e = ebtErrorNorms(H); md = ebtMassCheck(H); st = H.stats.stages{end};
  errs(i) = e.linf; facs(i) = st.avgFactor; mass(i) = md.rel;
end
r = log2(errs(1:end-1)./errs(2:end));
fprintf('  %-12s errs=[%.3e %.3e %.3e] rates=[%.2f %.2f] fac<=%.3f mass<=%.1e\n', ...
  name, errs, r, max(facs), max(mass));
ok = all(r > rateBar) && max(facs) <= facBar && max(mass) <= massBar;
if ~ok
  fprintf('  ** %s FAIL (rates/factor/mass)\n', name);
  PASS = false;
end
out = struct('errs', errs, 'rates', r, 'facs', facs, 'mass', mass);
end

% ------------------------------------------------------------------------
function H = buildBoxesT(prob, m, levels, rtol, maxv)
%buildBoxesT Assemble a tree hierarchy from explicit per-level box lists
%   (boxes in the CURRENT-finest-level index space).  The single parent of
%   each box is resolved by containment; every hierarchy is validated with
%   ebtTreeCheck before it is returned.
opts = ebtOptions('BaseCells', [m m], 'Verbose', 0, 'RelTol', rtol, ...
                  'MaxVCycles', maxv);
H = ebtSetup(prob, opts);
for L = 1:numel(levels)
  bx = levels{L};
  for a = 1:size(bx,1)
    for b = a+1:size(bx,1)
      if bx(a,1) <= bx(b,2) && bx(a,2) >= bx(b,1) && ...
         bx(a,3) <= bx(b,4) && bx(a,4) >= bx(b,3)
        error('stress:overlap', 'level set %d: boxes %d,%d overlap', L, a, b);
      end
    end
  end
  levC = H.lev{H.nlev};
  par = zeros(size(bx,1), 1);
  for i = 1:size(bx,1)
    for q = 1:levC.np
      qb = levC.box(q,:);
      if bx(i,1) >= qb(1) && bx(i,2) <= qb(2) && ...
         bx(i,3) >= qb(3) && bx(i,4) <= qb(4)
        par(i) = q; break;
      end
    end
    if par(i) == 0
      error('stress:parent', 'box %d of level set %d not inside one parent', i, L);
    end
  end
  H = ebtAddLevel(H, [par, bx]);
end
ebtTreeCheck(H);
H = ebtSolveCycles(H);
end

% ------------------------------------------------------------------------
function H = buildBoxesRing(prob, m, levels, rtol, maxv)
%buildBoxesRing Same hierarchy through the EllipticBSAM (ring) solver.
opts = ebOptions('BaseCells', [m m], 'Verbose', 0, 'RelTol', rtol, ...
                 'MaxVCycles', maxv);
H = ebSetup(prob, opts);
for L = 1:numel(levels)
  H = ebAddLevel(H, levels{L});
end
H = ebSolveCycles(H);
end

% ------------------------------------------------------------------------
function b = scaleBoxes(b, s)
b = [(b(:,1)-1)*s+1, b(:,2)*s, (b(:,3)-1)*s+1, b(:,4)*s];
end

% ------------------------------------------------------------------------
function d = maxQdiff(Hr, Ht)
%maxQdiff Max abs difference of the two solvers' patch interiors.
d = 0;
for k = Hr.rootIdx:Hr.nlev
  for p = 1:Hr.lev{k}.np
    A = Hr.lev{k}.Q{p}(2:end-1, 2:end-1);
    B = Ht.lev{k}.Q{p}(2:end-1, 2:end-1);
    d = max(d, max(abs(A(:) - B(:))));
  end
end
end

% ------------------------------------------------------------------------
function dup = coeffDupMismatch(H)
%coeffDupMismatch Max mismatch between the two stored copies of every
%   x/y-face shared by abutting same-level patches, across all levels.
%   Nonzero mismatch after ebtCoeffFilldown means pushFaceX/pushFaceY did
%   not route an averaged boundary face into the abutting sibling's copy
%   (the historical T8 conservation-bug signature; variable D only).
dup = 0;
for k = 1:H.nlev
  lev = H.lev{k};
  for a = 1:lev.np
    ab = lev.box(a, :);
    for b = 1:lev.np
      if a == b, continue; end
      bb = lev.box(b, :);
      if ab(2) + 1 == bb(1)                       % a east | b west
        j0 = max(ab(3), bb(3)); j1 = min(ab(4), bb(4));
        if j0 <= j1
          va = lev.DeW{a}(ab(2)-ab(1)+2, j0-ab(3)+1 : j1-ab(3)+1);
          vb = lev.DeW{b}(1, j0-bb(3)+1 : j1-bb(3)+1);
          dup = max(dup, max(abs(va - vb)));
        end
      end
      if ab(4) + 1 == bb(3)                       % a north | b south
        i0 = max(ab(1), bb(1)); i1 = min(ab(2), bb(2));
        if i0 <= i1
          va = lev.DnS{a}(i0-ab(1)+1 : i1-ab(1)+1, ab(4)-ab(3)+2);
          vb = lev.DnS{b}(i0-bb(1)+1 : i1-bb(1)+1, 1);
          dup = max(dup, max(abs(va - vb)));
        end
      end
    end
  end
end
end
