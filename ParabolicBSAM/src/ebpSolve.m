function sol = ebpSolve(pprob, varargin)
%EBPSOLVE Time-dependent block-structured adaptive multigrid solve of
%     du/dt - div(D grad u) + C u = f      (ParabolicBSAM / EBSAM-P).
%
%   sol = ebpSolve(pprob [, opts | 'Name', value, ...])
%
%   Implicit stepping (BE / BDF2 / CN) where every step is one
%   mass-conservative AFAS solve of the EBSAM elliptic kernel with
%   C_eff = C + alpha/Dt; the hierarchy follows the solution via periodic
%   conservative regridding.  Warm starts keep the per-step cost at a few
%   V-cycles.
%
%   Returns sol with: H (final hierarchy), t, mass(t) (composite integral
%   of u), cycles(n), relres(n), levels/cells history at regrids, err
%   (vs pprob.exact at TFinal, if given), stats.

if numel(varargin) == 1 && isstruct(varargin{1})
  opts = varargin{1};
else
  opts = ebpOptions(varargin{:});
end

nsteps = max(1, round(opts.TFinal / opts.Dt));
Dt = opts.TFinal / nsteps;
alphaOf = @(s) 1*strcmp(s,'be') + 1.5*strcmp(s,'bdf2') + 2*strcmp(s,'cn');
scheme1 = opts.Scheme;
if strcmp(scheme1, 'bdf2'), scheme1 = 'be'; end   % bootstrap

H = ebpInit(pprob, opts, alphaOf(scheme1)/Dt);

mass = zeros(1, nsteps+1);
mass(1) = compMass(H);
cycles = zeros(1, nsteps);
relres = zeros(1, nsteps);
nregrids = 0;

for n = 1:nsteps
  t1 = n * Dt;
  if n == 1
    scheme = scheme1;
  else
    scheme = opts.Scheme;
  end

  % Assemble CN's explicit term with the OLD coefficients and boundary
  % data, before refreshing the implicit operator to the new time.
  H = ebpRhs(H, pprob, t1, Dt, scheme);
  H = ebpSetTime(H, pprob, t1, alphaOf(scheme)/Dt);

  H = ebSolveCycles(H);
  st = H.stats.stages{end};
  cycles(n) = st.vcycles;
  relres(n) = st.finalRel;

  % advance the history (post fill-down state: composite-consistent)
  H.par.Uoo = H.par.Uo;
  H.par.Uo = interiors(H);
  mass(n+1) = compMass(H);

  if ~isempty(opts.Callback)
    opts.Callback(H, t1, n);
  end

  if opts.StepVerbose > 0 && (mod(n, opts.StepVerbose) == 0 || n == nsteps)
    nc = 0;
    for k = H.rootIdx:H.nlev
      for p = 1:H.lev{k}.np, nc = nc + numel(H.lev{k}.FTrue{p}); end
    end
    fprintf('[ebsam-p] step %4d/%d  t=%.4f  V=%d  rel=%.1e  levels=%d  cells=%d  mass=%.6e\n', ...
      n, nsteps, t1, cycles(n), relres(n), H.nlev - H.rootIdx + 1, nc, mass(n+1));
  end

  if opts.RegridEvery > 0 && mod(n, opts.RegridEvery) == 0 ...
      && n < nsteps && opts.eb.MaxLevels > 1
    H = ebpRegrid(H, pprob, t1, alphaOf(opts.Scheme)/Dt, opts);
    nregrids = nregrids + 1;
  end
end

sol.H = H;
sol.t = (0:nsteps) * Dt;
sol.mass = mass;
sol.cycles = cycles;
sol.relres = relres;
sol.nregrids = nregrids;
sol.stats = H.stats;
if ~isempty(pprob.exact)
  T = opts.TFinal;
  sol.err = ebErrorNorms(H, @(x, y) pprob.exact(x, y, T));
else
  sol.err = [];
end
for k = H.rootIdx:H.nlev
  lev = H.lev{k};
  ncell = 0;
  for p = 1:lev.np, ncell = ncell + numel(lev.FTrue{p}); end
  sol.levels(k - H.rootIdx + 1) = struct('n', lev.n, 'np', lev.np, ...
    'cells', ncell);
end
end

% ------------------------------------------------------------------------
function m = compMass(H)
[mbar, vol] = ebCompositeMean(H);
m = mbar * vol;
end

function U = interiors(H)
U = cell(1, H.nlev);
for k = 1:H.nlev
  for p = 1:H.lev{k}.np
    U{k}{p} = H.lev{k}.Q{p}(2:end-1, 2:end-1);
  end
end
end
