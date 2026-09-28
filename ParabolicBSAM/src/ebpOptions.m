function opts = ebpOptions(varargin)
%EBPOPTIONS Options for the parabolic solver.  Time-stepping fields:
%
%   Scheme       'bdf2' (default; BE bootstrap on step 1) | 'be' | 'cn'
%   Dt           requested positive time step (required)
%   TFinal       positive final time (required). A constant step is chosen
%                as TFinal/max(1,round(TFinal/Dt)) to land exactly at TFinal.
%   RegridEvery  regrid the moving hierarchy every k steps (8; 0 = never)
%   StepVerbose  print a one-line summary every k steps (10; 0 = silent)
%   Callback     @(H, t, n) called after every accepted step (e.g. frames)
%   MaxVCyclesStep  V-cycle cap per step (20)
%
%   All remaining name-value pairs are forwarded to ebOptions (BaseCells,
%   MaxLevels, Indicator/TagThreshold/TagBuffer, RelTol, GhostInterp,
%   MassCorrection, NPre/NPost, ...).  The elliptic Verbose is forced to
%   0; use StepVerbose for time-loop output.

timeDefaults = struct('Scheme', 'bdf2', 'Dt', [], 'TFinal', [], ...
  'RegridEvery', 8, 'StepVerbose', 10, 'Callback', [], ...
  'MaxVCyclesStep', 20);

tf = fieldnames(timeDefaults);
rest = {};
opts = timeDefaults;
if mod(numel(varargin), 2) ~= 0
  error('ebp_options:pairs', 'Options must be name-value pairs.');
end
i = 1;
while i <= numel(varargin)
  nm = varargin{i};
  hit = find(strcmpi(nm, tf), 1);
  if ~isempty(hit)
    opts.(tf{hit}) = varargin{i+1};
  else
    rest(end+1:end+2) = varargin(i:i+1);
  end
  i = i + 2;
end

opts.Scheme = lower(opts.Scheme);
if ~ismember(opts.Scheme, {'be', 'bdf2', 'cn'})
  error('ebp_options:scheme', 'Scheme must be be, bdf2 or cn.');
end
if isempty(opts.Dt) || isempty(opts.TFinal)
  error('ebp_options:time', 'Dt and TFinal are required.');
end
validateattributes(opts.Dt, {'numeric'}, {'real','scalar','finite','positive'}, mfilename, 'Dt');
validateattributes(opts.TFinal, {'numeric'}, {'real','scalar','finite','positive'}, mfilename, 'TFinal');
validateattributes(opts.RegridEvery, {'numeric'}, {'real','scalar','finite','integer','nonnegative'}, mfilename, 'RegridEvery');
validateattributes(opts.StepVerbose, {'numeric'}, {'real','scalar','finite','integer','nonnegative'}, mfilename, 'StepVerbose');
validateattributes(opts.MaxVCyclesStep, {'numeric'}, {'real','scalar','finite','integer','positive'}, mfilename, 'MaxVCyclesStep');
if ~isempty(opts.Callback) && ~isa(opts.Callback, 'function_handle')
  error('ebp_options:callback', 'Callback must be empty or a function handle.');
end

opts.eb = ebOptions(rest{:});
opts.eb.Verbose = 0;
opts.eb.MaxVCycles = opts.MaxVCyclesStep;
end
