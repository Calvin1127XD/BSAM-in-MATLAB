function bsamSetup
%BSAMSETUP Add the four solver packages and portfolio examples to the path.
% Test folders are deliberately separate: ring and tree suites each contain
% their own mmsCases.m. Run each suite from its own tests directory.
root = fileparts(mfilename('fullpath'));
for folder = {'EllipticBSAM','EllipticBSAMTree','ParabolicBSAM','DiffuseDomainBSAM'}
    addpath(fullfile(root, folder{1}, 'src'));
end
addpath(fullfile(root, 'examples'));
fprintf('BSAM in MATLAB: elliptic (QC-ring/tree), parabolic, diffuse domain.\n');
end
