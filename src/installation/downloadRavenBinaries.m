function downloadRavenBinaries(tools)
% downloadRavenBinaries
%   Download RAVEN's external command-line binaries on demand from the
%   toolbox-assets repository into the src/external directory of RAVEN. Only
%   the binary matching the current platform is fetched, and only when it is
%   missing.
%
%   Input:
%   tools   cell array, subset of {'blast+','diamond','hmmer'} to
%           download. Defaults to all of them. Tools whose target executable is
%           already present are skipped.
%
%   The download requires an internet connection. An offline installation
%   uses the offline bundle (RAVEN-<tag>-full.zip) attached to each RAVEN
%   release, which contains the binaries of all platforms in src/external.
%
%   The binaries are hosted as release assets in
%   https://github.com/CompMetEng/toolbox-assets.
%
%   Usage: downloadRavenBinaries({'blast+','diamond'})
%
% getBlast, getDiamond and getKEGGModelForOrganism call this function
% when a binary is missing. checkRaven calls it for all three tools when its
% checkBinaries input is true.

if nargin < 1 || isempty(tools)
    tools = {'blast+','diamond','hmmer'};
elseif ischar(tools)
    tools = {tools};
end

ravenDir = findRAVENroot();

% toolbox-assets platform key for the per-platform binary ZIPs.
if ispc
    plat = 'windows-x86_64';
elseif ismac
    plat = 'macos-arm64';
else
    plat = 'linux-x86_64';
end

for i = 1:numel(tools)
    tool = tools{i};
    switch tool
        case 'blast+'
            tag = 'blast-2.17.0';
            asset = ['blast-2.17.0-' plat '.zip'];
            destDir = fullfile(ravenDir,'src','external','blast+');
            sentinel = fullfile(destDir,'blastp');           % bare name on every OS (.exe on Windows)
            execs = {'blastp','makeblastdb'};
        case 'diamond'
            tag = 'diamond-2.1.17';
            asset = ['diamond-2.1.17-' plat '.zip'];
            destDir = fullfile(ravenDir,'src','external','diamond');
            sentinel = fullfile(destDir,'diamond');
            execs = {'diamond'};
        case 'hmmer'
            if ispc
                tag = 'hmmer-3.3.2';                          % native Windows build (Cygwin); reads HMMER3/f
                asset = ['hmmer-3.3.2-' plat '.zip'];
            else
                tag = 'hmmer-3.4.0';
                asset = ['hmmer-3.4.0-' plat '.zip'];
            end
            destDir = fullfile(ravenDir,'src','external','hmmer');
            sentinel = fullfile(destDir,'hmmsearch');
            execs = {'hmmsearch'};
        otherwise
            error('downloadRavenBinaries: unknown tool ''%s''.',tool);
    end

    % macOS uses a .mac suffix and Windows a .exe suffix (Linux: bare name).
    if ispc
        sentinel = [sentinel '.exe']; %#ok<AGROW>
    elseif ismac
        sentinel = [sentinel '.mac']; %#ok<AGROW>
    end
    if exist(sentinel,'file')
        continue;   % already provisioned
    end

    zipPath = fetchToolboxAsset(fullfile(ravenDir,'src','external'),tag,asset);
    unzip(zipPath,destDir);
    delete(zipPath);

    % macOS uses the .mac suffix throughout RAVEN, but the toolbox-assets macOS ZIP
    % ships bare names, so rename the extracted executables to match.
    if ismac
        for k = 1:numel(execs)
            bare = fullfile(destDir,execs{k});
            if isfile(bare)
                movefile(bare, fullfile(destDir,[execs{k} '.mac']));
            end
        end
    end

    % MATLAB's unzip does not preserve Unix execute permissions; restore them.
    if isunix
        if ismac; sfx = '.mac'; else; sfx = ''; end
        for k = 1:numel(execs)
            system(['chmod +x "' fullfile(destDir,[execs{k} sfx]) '"']);
        end
        % Downloaded binaries are quarantined by macOS Gatekeeper and will not
        % execute until the com.apple.quarantine attribute is cleared.
        if ismac
            system(['xattr -dr com.apple.quarantine "' destDir '"']);
        end
    end
end
end
