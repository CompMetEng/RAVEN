function [ravenPath, prevDir] = findRAVENroot()
% findRAVENroot
%   Finds the root of the RAVEN directory by searching for the path to
%   version.txt. Can also record the current directory, for functions that use
%   the ravenPath to navigate to a precise folder and return to the previous
%   directory afterwards. See, for example, optimizeProb calling glpk.
%
%   The RAVEN.ravenPath preference is used if it points to a directory with a
%   version.txt file. Otherwise the search starts at the location of this
%   function file and moves up the directory tree.
%
%   The resolved path is cached for the rest of the MATLAB session, so the
%   preference is read once although the function is called on every solver
%   invocation. `clear findRAVENroot` drops the cached path after a change of
%   the preference (checkRaven does this) without a restart of MATLAB.

persistent cachedPath
prevDir = pwd();
if ~isempty(cachedPath)
    ravenPath = cachedPath;
    return;
end

% A stored preference is only trusted if it still points at a real RAVEN
% install; otherwise fall through to walking up from the currently
% executing copy of this file. Without this check, a stale preference
% left over from a different RAVEN checkout on the same machine silently
% resolves to that other copy's data, not the one actually running.
ravenPath = '';
if ispref('RAVEN','ravenPath')
    prefPath = getpref('RAVEN','ravenPath');
    if isfile(fullfile(prefPath,'version.txt'))
        ravenPath = prefPath;
    end
end
if isempty(ravenPath)
    ST=dbstack('-completenames');
    ravenPath = ST(strcmp({ST.name},'findRAVENroot')).file;
    rootFound = 0;
    while rootFound == 0
        isRoot = isfile(fullfile(ravenPath,'version.txt'));
        if isRoot
            rootFound = 1;
        else
            ravenPathOld = ravenPath;
            ravenPath = fileparts(ravenPath);
            if strcmp(ravenPathOld,ravenPath)
                error('Cannot find the RAVEN root directory. Make sure you have not removed the version.txt file from your RAVEN installation.')
            end
        end
    end
end
cachedPath = ravenPath;
end
