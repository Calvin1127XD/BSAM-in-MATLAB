function H = ebVcycle(H, k)
%ebVcycle Recursive AFAS V-cycle on the composite hierarchy (paper
%   algorithm a0-a9 with the recursive coarse solve of Remark 3.7):
%
%     relax(k)                                  a1
%     restrict u_k -> QC, parent                a3
%     fill parent ghosts; gather QC ring
%     coarse load:  R(f - L u) + L_c(R u)       a2,a4
%                   + flux-balance corrections  a4 (mass conservation)
%     recurse k-1                               a5
%     correct: u_k += P(parent - QC)            a6,a7
%     relax(k)                                  a8
%
%   Level 1 is solved directly (prefactored sparse system).

% Bottom of the cycle: direct sparse solve on the global coarsest grid.
if k == 1
  H = ebCoarsest(H, 'solve');
  return;
end

% Pre-smooth, then form the coarse problem: restrict the solution, refresh
% the parent ghosts, gather the QC ring, and build the FAS coarse load
% (with flux-balance corrections at the CF interface).
H = ebRelaxLevel(H, k, H.opts.NPre);
H = ebTransfer(H, k, 'restrict');
H = ebFillGhosts(H, k-1);
H = ebTransfer(H, k, 'ring');
H = ebCoarseLoad(H, k);
% Solve the coarse problem recursively, then apply the prolongated
% correction and post-smooth.
H = ebVcycle(H, k-1);
H = ebTransfer(H, k, 'correct');
H = ebFillGhosts(H, k);
H = ebRelaxLevel(H, k, H.opts.NPost);
end
