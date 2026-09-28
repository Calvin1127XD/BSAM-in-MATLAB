function opts = ebtOptions(varargin)
%ebtOptions Solver options for the block-structured adaptive FAS multigrid.
%
%   opts = ebtOptions('Name', value, ...)
%
%   Grid / hierarchy
%     BaseCells      [nx ny] root grid (default [64 64]); each must be
%                    even enough to coarsen down to CoarsestCells.
%     MaxLevels      number of AMR levels including the root (default 1).
%     CoarsestCells  stop sub-root coarsening when min(n) <= this (16).
%
%   Discretization / interface
%     GhostInterp    'qif' (quadratic, default) | 'lif' (linear).
%     MassCorrection true (default): flux-balanced corrections to the
%                    coarse-level force (the conservative scheme).  Set
%                    false to recover the non-conservative "QI" scheme.
%
%   Multigrid
%     NPre, NPost    smoothing sweeps per level (default 3, 3).
%     Omega          relaxation factor for red-black GS (default 1.0).
%     MaxVCycles     per solve stage (default 60).
%     RelTol         relative composite-residual reduction target (1e-9).
%
%   Adaptivity
%     Indicator      'ulap' undivided Laplacian of u (paper eq. 4.32),
%                    or 'fun' to use TagFunction.
%     TagFunction    @(X,Y) -> logical, cells to refine (used when
%                    Indicator='fun'); X,Y are cell centers.
%     TagThreshold   threshold C_k for 'ulap' (scalar or @(level)) (1e-3).
%     TagBuffer      safety buffer in parent cells around tags (2).
%     MinBoxWidth    minimum cluster box width in parent cells (4).
%     MinFillRatio   Berger-Rigoutsos efficiency target (0.70).
%
%   Misc
%     Verbose        0 silent, 1 per-stage, 2 per-V-cycle (default 1).
%     InitialGuess   @(X,Y) initial iterate (default 0).
%     MeanZero       fix the additive constant for pure-Neumann problems
%                    ('auto' default: on iff all-Neumann and C==0).
%     SourceQuad     source quadrature points per direction (default 1 =
%                    midpoint; nq > 1 uses the tensor Gauss-Legendre cell
%                    average, making the composite source integral exact
%                    to order 2*nq -- see srcquad in ebMakeLevel).

p = inputParser;
p.addParameter('BaseCells', [64 64]);
p.addParameter('MaxLevels', 1);
p.addParameter('CoarsestCells', 16);
p.addParameter('GhostInterp', 'qif');
p.addParameter('MassCorrection', true);
p.addParameter('NPre', 3);
p.addParameter('NPost', 3);
p.addParameter('Omega', 1.0);
p.addParameter('MaxVCycles', 60);
p.addParameter('RelTol', 1e-9);
p.addParameter('Indicator', 'ulap');
p.addParameter('TagFunction', []);
p.addParameter('TagThreshold', 1e-3);
p.addParameter('TagBuffer', 2);
p.addParameter('MinBoxWidth', 4);
p.addParameter('MinFillRatio', 0.70);
p.addParameter('Verbose', 1);
p.addParameter('InitialGuess', []);
p.addParameter('MeanZero', 'auto');
p.addParameter('SourceQuad', 1);
p.parse(varargin{:});
opts = p.Results;

opts.BaseCells = round(opts.BaseCells(:).');
if numel(opts.BaseCells) == 1, opts.BaseCells = [opts.BaseCells opts.BaseCells]; end
opts.GhostInterp = lower(opts.GhostInterp);
if ~ismember(opts.GhostInterp, {'qif', 'lif'})
  error('ebt_options:ghost', 'GhostInterp must be qif or lif');
end
end
