function H = ebCoarsest(H, mode)
%ebCoarsest Direct solver for level 1 (single full-domain patch).
%   ebCoarsest(H,'assemble') builds and factors the ghost-eliminated
%   5-point system (quadratic Dirichlet ghosts / Neumann fluxes folded in).
%   For the singular pure-Neumann, C==0 case an augmented KKT system pins
%   the mean.  ebCoarsest(H,'solve') solves L u = F (the FAS coarse
%   equation) and refreshes the level-1 ghosts.

lev = H.lev{1};
n1 = lev.n(1); n2 = lev.n(2);
N = n1 * n2;

switch mode
  case 'assemble'
    hx = lev.h(1); hy = lev.h(2);
    hx2 = hx*hx; hy2 = hy*hy;
    DeW = lev.DeW{1}; DnS = lev.DnS{1}; CC = lev.CC{1};
    bcd = lev.bc{1};

    % Assemble in triplet form; start with the diagonal reaction term CC.
    idx = reshape(1:N, n1, n2);
    rows = idx(:); colsd = idx(:); vals = CC(:);
    bcRhs = zeros(N, 1);

    % internal x-faces
    % Each interior face contributes +c to both diagonals and -c to the
    % two off-diagonal couplings (conservative flux difference).
    c = DeW(2:n1, :) / hx2;
    iL = idx(1:n1-1, :); iR = idx(2:n1, :);
    rows = [rows; iL(:); iR(:); iL(:); iR(:)];
    colsd = [colsd; iL(:); iR(:); iR(:); iL(:)];
    vals = [vals; c(:); c(:); -c(:); -c(:)];
    % internal y-faces
    c = DnS(:, 2:n2) / hy2;
    iL = idx(:, 1:n2-1); iR = idx(:, 2:n2);
    rows = [rows; iL(:); iR(:); iL(:); iR(:)];
    colsd = [colsd; iL(:); iR(:); iR(:); iL(:)];
    vals = [vals; c(:); c(:); -c(:); -c(:)];

    % west boundary (cells (1,:), face DeW(1,:))
    d = DeW(1, :).' ; i1 = idx(1, :).'; i2 = idx(2, :).'; g = bcd(1).g(:);
    [rows, colsd, vals, bcRhs] = bface(rows, colsd, vals, bcRhs, ...
      bcd(1).code, d, g, i1, i2, hx, hx2);
    % east
    d = DeW(n1+1, :).'; i1 = idx(n1, :).'; i2 = idx(n1-1, :).'; g = bcd(2).g(:);
    [rows, colsd, vals, bcRhs] = bface(rows, colsd, vals, bcRhs, ...
      bcd(2).code, d, g, i1, i2, hx, hx2);
    % south
    d = DnS(:, 1); i1 = idx(:, 1); i2 = idx(:, 2); g = bcd(3).g(:);
    [rows, colsd, vals, bcRhs] = bface(rows, colsd, vals, bcRhs, ...
      bcd(3).code, d, g, i1, i2, hy, hy2);
    % north
    d = DnS(:, n2+1); i1 = idx(:, n2); i2 = idx(:, n2-1); g = bcd(4).g(:);
    [rows, colsd, vals, bcRhs] = bface(rows, colsd, vals, bcRhs, ...
      bcd(4).code, d, g, i1, i2, hy, hy2);

    A = sparse(rows, colsd, vals, N, N);

    % singular iff all-Neumann and C == 0 across the hierarchy
    % (then u is defined only up to a constant; pin the mean with a
    % one-row/one-column KKT augmentation instead of regularizing A).
    allneu = all(arrayfun(@(s) s.code ~= 1, bcd));
    cmax = 0;
    for k = 1:H.nlev
      for p = 1:H.lev{k}.np
        cmax = max(cmax, max(abs(H.lev{k}.CC{p}(:))));
      end
    end
    singular = allneu && cmax == 0;
    if singular
      e = ones(N, 1);
      H.coarsest.dec = decomposition([A, e; e.', 0]);
    else
      H.coarsest.dec = decomposition(A);
    end
    H.coarsest.singular = singular;
    H.coarsest.bcRhs = bcRhs;

  case 'solve'
    rhs = H.lev{1}.F{1}(:) + H.coarsest.bcRhs;
    if H.coarsest.singular
      % Augmented system: last unknown is the Lagrange multiplier that
      % enforces zero mean; discard it after the solve.
      z = H.coarsest.dec \ [rhs; 0];
      u = z(1:N);
    else
      u = H.coarsest.dec \ rhs;
    end
    H.lev{1}.Q{1}(2:end-1, 2:end-1) = reshape(u, n1, n2);
    H = ebFillGhosts(H, 1);
end
end

% ------------------------------------------------------------------------
function [rows, cols, vals, bcRhs] = bface(rows, cols, vals, bcRhs, ...
                                           code, d, g, i1, i2, h, h2)
%BFACE Fold one physical side's boundary faces into the matrix/RHS.
%   Quadratic Dirichlet ghost u_g = (8/3)g - 2u_1 + (1/3)u_2 gives
%     diag += 3 D/h^2,  A(i1,i2) -= D/(3 h^2),  rhs += (8/3) D g/h^2.
%   Neumann (du/dn = g): rhs += D g / h.
if code == 1
  rows = [rows; i1; i1];
  cols = [cols; i1; i2];
  vals = [vals; 3*d/h2; -d/(3*h2)];
  bcRhs(i1) = bcRhs(i1) + (8/3) * d .* g / h2;
elseif code == 2
  bcRhs(i1) = bcRhs(i1) + d .* g / h;
end
end
