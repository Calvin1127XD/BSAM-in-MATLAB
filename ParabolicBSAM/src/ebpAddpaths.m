function ebpAddpaths()
%EBPADDPATHS Put ParabolicBSAM and the EBSAM elliptic kernel on the path.
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', 'EllipticBSAM', 'src'));
end
