function sol = ddSolve(geom, epsi, box, bctype, data, varargin)
%ddSolve Solve an elliptic problem on an irregular domain by the diffuse
%   domain method with block-structured adaptive multigrid refinement of
%   the interfacial layer.
%
%   sol = ddSolve(geom, eps, box, bctype, data, 'Name', value, ...)
%
%   Accepts all ddProblem options plus all ebOptions options, plus:
%     'TagWidth'  (2.5)  refine where |psi| < TagWidth * eps
%   Defaults tuned for the layer: Indicator='fun' with the geometric band
%   tag, MaxLevels as given.  Returns the ebSolve solution augmented with
%     sol.aux        ddProblem auxiliaries (phi, psi, inside, ...)
%     sol.errInside  error vs data.exact measured on Omega_1 (psi < 0)

% split varargin into ddProblem's and ours/eb's
ddNames = {'PhiFloor', 'Mu', 'OuterBC'};
ddArgs = {}; rest = {};
tagWidth = 2.5;
i = 1;
if mod(numel(varargin), 2) ~= 0
  error('dd_solve:pairs', 'Options must be name-value pairs.');
end
while i <= numel(varargin)
  nm = varargin{i};
  if any(strcmpi(nm, ddNames))
    ddArgs(end+1:end+2) = varargin(i:i+1);
  elseif strcmpi(nm, 'TagWidth')
    tagWidth = varargin{i+1};
  else
    rest(end+1:end+2) = varargin(i:i+1);
  end
  i = i + 2;
end
validateattributes(tagWidth, {'numeric'}, {'real','scalar','finite','positive'}, mfilename, 'TagWidth');

[prob, aux] = ddProblem(geom, epsi, box, bctype, data, ddArgs{:});

opts = ebOptions('Indicator', 'fun', 'TagFunction', aux.tagfun(tagWidth), ...
                  rest{:});

sol = ebSolve(prob, opts);
sol.aux = aux;
if isfield(data, 'exact') && ~isempty(data.exact)
  sol.errInside = ddErrorInside(sol.H, data.exact, aux);
else
  sol.errInside = [];
end
end
