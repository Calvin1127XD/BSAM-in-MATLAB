function [boxes, info] = ebTagCluster(H)
%ebTagCluster Tag cells of the current finest level and cluster them
%   into rectangular boxes (Berger-Rigoutsos signature splitting).
%
%   Tagging: 'ulap'  - undivided Laplacian of u exceeds TagThreshold
%                      (paper eq. 4.32),
%            'rte'   - relative truncation error |L_2h(Ru) - R(L_h u)|
%                      exceeds TagThreshold (paper sec. 4.4; ebRte), or
%            'fun'   - user TagFunction(X,Y) -> logical at cell centers.
%
%   Optional second output (for visualization/diagnostics):
%     info.indicator{p}  per-patch indicator field (or logical for 'fun')
%     info.Traw          tag mask before buffering (level index space)
%     info.T             tags after buffer + proper-nesting clip
%     info.allowed       the proper-nesting region
%   Tags are buffered by TagBuffer cells, then clipped to the proper-
%   nesting region (cells whose full 8-neighborhood is covered by the
%   current finest level), so child patches always have valid coarse-fine
%   interpolation stencils and one-level interface jumps (1-conforming).
%   Boxes are returned in level-k index space; the child patches are their
%   2x refinement, hence automatically aligned to parent cells, with even
%   dimensions and at least MinBoxWidth parent cells across.

opts = H.opts;
k = H.nlev;
lev = H.lev{k};
n = lev.n;

% ---------------- tags
T = false(n);
dom = H.prob.domain;
wantInfo = nargout > 1;
info = struct('indicator', {{}}, 'Traw', [], 'T', [], 'allowed', []);
for p = 1:lev.np
  b = lev.box(p, :);
  switch lower(opts.Indicator)
    case 'fun'
      if isempty(opts.TagFunction)
        error('eb_tag_cluster:tagfun', 'Indicator=fun needs TagFunction.');
      end
      xc = dom(1) + ((b(1):b(2)) - 0.5) * lev.h(1);
      yc = dom(3) + ((b(3):b(4)) - 0.5) * lev.h(2);
      [XC, YC] = ndgrid(xc, yc);
      if nargin(opts.TagFunction) >= 3
        % L = current finest AMR level (1-based); refining L creates L+1
        tp = logical(opts.TagFunction(XC, YC, k - H.rootIdx + 1));
      else
        tp = logical(opts.TagFunction(XC, YC));
      end
    case 'ulap'
      Q = lev.Q{p};
      mx = b(2) - b(1) + 1; my = b(4) - b(3) + 1;
      sdx = Q(3:mx+2, 2:my+1) - 2*Q(2:mx+1, 2:my+1) + Q(1:mx, 2:my+1);
      sdy = Q(2:mx+1, 3:my+2) - 2*Q(2:mx+1, 2:my+1) + Q(2:mx+1, 1:my);
      s = sqrt(sdx.^2 + sdy.^2);
      thr = opts.TagThreshold;
      if isa(thr, 'function_handle'), thr = thr(k - H.rootIdx + 1); end
      tp = s > thr;
      if wantInfo, info.indicator{p} = s; end
    case 'rte'
      s = ebRte(H, k, p);
      thr = opts.TagThreshold;
      if isa(thr, 'function_handle'), thr = thr(k - H.rootIdx + 1); end
      tp = s > thr;
      if wantInfo, info.indicator{p} = s; end
    otherwise
      error('eb_tag_cluster:ind', 'unknown Indicator %s', opts.Indicator);
  end
  if wantInfo && strcmpi(opts.Indicator, 'fun'), info.indicator{p} = tp; end
  T(b(1):b(2), b(3):b(4)) = T(b(1):b(2), b(3):b(4)) | tp;
end
if wantInfo, info.Traw = T; end

% ---------------- buffer (dilation by TagBuffer)
T = dilate(T, opts.TagBuffer);

% ---------------- proper-nesting clip
cover = false(n);
for p = 1:lev.np
  b = lev.box(p, :);
  cover(b(1):b(2), b(3):b(4)) = true;
end
allowed = erode8(cover);
T = T & allowed;
if wantInfo
  info.T = T;
  info.allowed = allowed;
end

if ~any(T(:))
  boxes = zeros(0, 4);
  return;
end

% ---------------- Berger-Rigoutsos clustering
boxes = brCluster(T, opts.MinBoxWidth, opts.MinFillRatio);

% clip boxes to the allowed region (bisect until clean), drop tagless
boxes = clipAllowed(boxes, T, allowed);

% grow boxes below the minimum width (within allowed, no overlaps)
boxes = growMin(boxes, allowed, opts.MinBoxWidth, n);

% boxes still thinner than 2 cells (ungrowable after the nesting clip):
% merge into a neighbor when the bounding union stays allowed and
% non-overlapping, otherwise drop them (locally less refinement is
% preferred over fragile slivers - this is what made earlier attempts
% diverge).
boxes = mergeOrDropThin(boxes, allowed);
end

% ========================================================================
function M = dilate(M, b)
% Morphological dilation by a (2b+1)^2 square: OR of all shifted copies,
% computed on a zero-padded array so tags never wrap around edges.
if b <= 0, return; end
[n1, n2] = size(M);
P = false(n1 + 2*b, n2 + 2*b);
P(b+1:b+n1, b+1:b+n2) = M;
acc = false(n1, n2);
for dx = -b:b
  for dy = -b:b
    acc = acc | P(b+1+dx : b+n1+dx, b+1+dy : b+n2+dy);
  end
end
M = acc;
end

% ------------------------------------------------------------------------
function A = erode8(C)
%ERODE8 allowed(c) = cell and all in-domain 8-neighbors covered.
[n1, n2] = size(C);
P = true(n1+2, n2+2);            % outside the domain counts as covered
P(2:n1+1, 2:n2+1) = C;
A = true(n1, n2);
for dx = -1:1
  for dy = -1:1
    A = A & P(2+dx : n1+1+dx, 2+dy : n2+1+dy);
  end
end
end

% ------------------------------------------------------------------------
function out = brCluster(T, minW, eff)
% Berger-Rigoutsos: recursively shrink each candidate box to the tag
% bounding box, accept it once the tag fill ratio reaches eff (or it
% cannot shrink further), otherwise split at the best signature cut and
% recurse on both halves (explicit stack instead of recursion).
out = zeros(0, 4);
stack = {tightBbox(T, [1 size(T,1) 1 size(T,2)])};
while ~isempty(stack)
  bb = stack{end}; stack(end) = [];
  if isempty(bb), continue; end
  bb = tightBbox(T, bb);
  if isempty(bb), continue; end
  w = [bb(2)-bb(1)+1, bb(4)-bb(3)+1];
  sub = T(bb(1):bb(2), bb(3):bb(4));
  fill = nnz(sub) / numel(sub);
  if fill >= eff || all(w <= minW)
    out(end+1, :) = bb; %#ok<AGROW>
    continue;
  end
  % Signatures = tag counts per row / per column; the split search looks
  % for holes or curvature jumps in these 1D profiles.
  sigx = sum(sub, 2);
  sigy = sum(sub, 1).';
  [sx, qx] = pickSplit(sigx, minW);
  [sy, qy] = pickSplit(sigy, minW);
  if qx == 0 && qy == 0
    % no admissible split: halve the longest dimension if possible
    if w(1) >= 2*minW || w(2) >= 2*minW
      if w(1) >= w(2)
        s = floor(w(1)/2);
        stack{end+1} = [bb(1), bb(1)+s-1, bb(3), bb(4)];
        stack{end+1} = [bb(1)+s, bb(2), bb(3), bb(4)];
      else
        s = floor(w(2)/2);
        stack{end+1} = [bb(1), bb(2), bb(3), bb(3)+s-1];
        stack{end+1} = [bb(1), bb(2), bb(3)+s, bb(4)];
      end
    else
      out(end+1, :) = bb; %#ok<AGROW>
    end
    continue;
  end
  % prefer holes (q=2) over inflections (q=1); tie -> longer dimension
  usex = qx > qy || (qx == qy && w(1) >= w(2));
  if qx == 0, usex = false; end
  if qy == 0, usex = true; end
  if usex
    stack{end+1} = [bb(1), bb(1)+sx-1, bb(3), bb(4)];
    stack{end+1} = [bb(1)+sx, bb(2), bb(3), bb(4)];
  else
    stack{end+1} = [bb(1), bb(2), bb(3), bb(3)+sy-1];
    stack{end+1} = [bb(1), bb(2), bb(3)+sy, bb(4)];
  end
end
end

% ------------------------------------------------------------------------
function bb = tightBbox(T, bb)
sub = T(bb(1):bb(2), bb(3):bb(4));
rows = any(sub, 2); cols = any(sub, 1);
if ~any(rows), bb = []; return; end
i1 = find(rows, 1, 'first'); i2 = find(rows, 1, 'last');
j1 = find(cols, 1, 'first'); j2 = find(cols, 1, 'last');
bb = [bb(1)+i1-1, bb(1)+i2-1, bb(3)+j1-1, bb(3)+j2-1];
end

% ------------------------------------------------------------------------
function [s, quality] = pickSplit(sig, minW)
%pickSplit Split index s (left = 1..s, right = s+1..end).
%   quality: 2 = at a hole (zero signature), 1 = at the strongest
%   Laplacian-of-signature jump, 0 = none admissible.
len = numel(sig);
s = 0; quality = 0;
% holes
holes = find(sig == 0);
holes = holes(holes > minW & holes <= len - minW);
if ~isempty(holes)
  [~, ii] = min(abs(holes - len/2));
  s = holes(ii); quality = 2;
  return;
end
if len < 2*minW, return; end
d2 = sig(3:len) - 2*sig(2:len-1) + sig(1:len-2);  % at 2..len-1
best = -1;
for i = 2:len-2                                    % jump between i, i+1
  j = abs(d2(i) - d2(i-1));
  cut = i;                                         % split after cell i
  if cut >= minW && (len - cut) >= minW && j > best
    best = j; s = cut;
  end
end
if s > 0, quality = 1; end
end

% ------------------------------------------------------------------------
function out = clipAllowed(boxes, T, allowed)
% Keep bisecting any box that pokes outside the proper-nesting region
% until every surviving piece is fully allowed; pieces without tags drop.
out = zeros(0, 4);
stack = num2cell(boxes, 2).';
while ~isempty(stack)
  bb = stack{end}; stack(end) = [];
  sub = T(bb(1):bb(2), bb(3):bb(4));
  if ~any(sub(:)), continue; end
  ok = allowed(bb(1):bb(2), bb(3):bb(4));
  if all(ok(:))
    out(end+1, :) = bb; %#ok<AGROW>
    continue;
  end
  w = [bb(2)-bb(1)+1, bb(4)-bb(3)+1];
  if w(1) >= w(2)
    s = floor(w(1)/2);
    stack{end+1} = [bb(1), bb(1)+s-1, bb(3), bb(4)];
    stack{end+1} = [bb(1)+s, bb(2), bb(3), bb(4)];
  else
    s = floor(w(2)/2);
    stack{end+1} = [bb(1), bb(2), bb(3), bb(3)+s-1];
    stack{end+1} = [bb(1), bb(2), bb(3)+s, bb(4)];
  end
end
end

% ------------------------------------------------------------------------
function boxes = growMin(boxes, allowed, minW, n)
% Widen undersized boxes one strip at a time, alternating sides, as long
% as the new strip stays inside the allowed region and off other boxes.
nb = size(boxes, 1);
for b = 1:nb
  for dim = 1:2
    lo = boxes(b, 2*dim-1); hi = boxes(b, 2*dim);
    nd = n(dim);
    side = 0;
    while hi - lo + 1 < minW
      grown = false;
      for tryside = [side, 1-side]
        if tryside == 0 && lo > 1
          cand = [lo-1, hi];
          if stripOk(boxes, allowed, b, dim, lo-1, n)
            lo = lo - 1; grown = true; break;
          end
        elseif tryside == 1 && hi < nd
          cand = [lo, hi+1]; %#ok<NASGU>
          if stripOk(boxes, allowed, b, dim, hi+1, n)
            hi = hi + 1; grown = true; break;
          end
        end
      end
      if ~grown, break; end
      side = 1 - side;
    end
    boxes(b, 2*dim-1) = lo; boxes(b, 2*dim) = hi;
  end
end
end

function boxes = mergeOrDropThin(boxes, allowed)
thin = @(bb) (bb(2)-bb(1)+1 < 2) || (bb(4)-bb(3)+1 < 2);
b = 1;
while b <= size(boxes, 1)
  if ~thin(boxes(b, :))
    b = b + 1;
    continue;
  end
  merged = false;
  % candidate partners ordered by union area; a merge is accepted only if
  % the bounding union adds little area beyond the two boxes (<= 50%
  % slack), so distant band pieces are never bridged across untagged
  % regions (which would swallow e.g. the interior of an annular band)
  nb = size(boxes, 1);
  areab = (boxes(b,2)-boxes(b,1)+1) * (boxes(b,4)-boxes(b,3)+1);
  ua = inf(nb, 1);
  for o = 1:nb
    if o == b, continue; end
    u = [min(boxes(b,1), boxes(o,1)), max(boxes(b,2), boxes(o,2)), ...
         min(boxes(b,3), boxes(o,3)), max(boxes(b,4), boxes(o,4))];
    ua(o) = (u(2)-u(1)+1) * (u(4)-u(3)+1);
  end
  [~, order] = sort(ua);
  for o = order.'
    if o == b || ~isfinite(ua(o)), continue; end
    % only merge with an ADJACENT box (sliver expanded by 1 touches it)
    sb = boxes(b, :); ob = boxes(o, :);
    if ~(sb(1)-1 <= ob(2) && sb(2)+1 >= ob(1) && ...
         sb(3)-1 <= ob(4) && sb(4)+1 >= ob(3))
      continue;
    end
    u = [min(boxes(b,1), boxes(o,1)), max(boxes(b,2), boxes(o,2)), ...
         min(boxes(b,3), boxes(o,3)), max(boxes(b,4), boxes(o,4))];
    areao = (boxes(o,2)-boxes(o,1)+1) * (boxes(o,4)-boxes(o,3)+1);
    if ua(o) > 1.5 * (areab + areao), continue; end
    A = allowed(u(1):u(2), u(3):u(4));
    if ~all(A(:)), continue; end
    clash = false;
    for t = 1:nb
      if t == b || t == o, continue; end
      tb = boxes(t, :);
      if u(1) <= tb(2) && u(2) >= tb(1) && u(3) <= tb(4) && u(4) >= tb(3)
        clash = true; break;
      end
    end
    if clash, continue; end
    boxes(o, :) = u;
    boxes(b, :) = [];
    merged = true;
    break;
  end
  if ~merged
    warning('eb_tag_cluster:dropthin', ...
      'dropping an ungrowable sliver cluster box (less local refinement); disable with warning(''off'',''eb_tag_cluster:dropthin'')');
    boxes(b, :) = [];
  end
end
end

% ------------------------------------------------------------------------
function ok = stripOk(boxes, allowed, b, dim, pos, n) %#ok<INUSD>
%stripOk Strip at index pos along dim of box b is allowed and outside
%   every other box.
bb = boxes(b, :);
if dim == 1
  strip = allowed(pos, bb(3):bb(4));
  s = [pos, pos, bb(3), bb(4)];
else
  strip = allowed(bb(1):bb(2), pos);
  s = [bb(1), bb(2), pos, pos];
end
ok = all(strip(:));
if ~ok, return; end
for o = 1:size(boxes, 1)
  if o == b, continue; end
  ob = boxes(o, :);
  if s(1) <= ob(2) && s(2) >= ob(1) && s(3) <= ob(4) && s(4) >= ob(3)
    ok = false;
    return;
  end
end
end
