function [currVer, installType] = checkInstallation(developMode, checkBinaries)
% checkInstallation  DEPRECATED. Use checkRaven instead.
%
%   Forwards all arguments and outputs to checkRaven.
%
%   This function runs before RAVEN is on the MATLAB path. It therefore
%   calls no other RAVEN function until checkRaven has added RAVEN to the
%   path, and it does not use the shared deprecationWarning function in
%   the utils folder.
%
% Usage: [currVer, installType] = checkInstallation(developMode, checkBinaries)

warning('RAVEN:deprecated', ...
    'checkInstallation is deprecated and will be removed in the next major release. Use checkRaven instead.')
if nargin < 1
    developMode = false;
end
if nargin < 2
    checkBinaries = true;
end
if nargout > 0
    [currVer, installType] = checkRaven(developMode, checkBinaries);
else
    checkRaven(developMode, checkBinaries);
end
end
