function ddAddpaths()
%ddAddpaths Put the DiffuseDomainBSAM and EllipticBSAM sources on the path.
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', 'EllipticBSAM', 'src'));
end
