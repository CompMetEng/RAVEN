function localPath=fetchToolboxAsset(destDir,releaseTag,assetName)
% fetchToolboxAsset  Download (if needed) a toolbox-assets release asset.
%
% Downloads assetName from the given release tag of
% https://github.com/CompMetEng/toolbox-assets into destDir
% and returns its local path. An already-downloaded copy is reused as-is.
% Generic across toolbox-assets' release families (KEGG artefacts, HMM
% libraries, binaries, ...): nothing here is specific to any one of
% them, so callers supply the release tag and asset name rather than
% this function hardcoding either.
%
% A download is checked against the SHA-256 checksum that the
% toolbox-assets manifest (release manifest-v3) lists for it. A file whose
% checksum differs is deleted and an error is raised. An asset that the
% manifest does not list, or a manifest that cannot be read, gives a
% warning and the file is kept unchecked.
%
% Parameters
% ----------
% destDir : char
%     directory to download into (created if it does not exist yet).
% releaseTag : char
%     the toolbox-assets release tag, e.g. 'kegg118'.
% assetName : char
%     the release asset file name, e.g. 'kegg118_core.tar.gz'.
%
% Returns
% -------
% localPath : char
%     path to the downloaded file.

if ~isfolder(destDir)
    mkdir(destDir);
end
localPath=fullfile(destDir,assetName);
if isfile(localPath)
    return;
end
url=['https://github.com/CompMetEng/toolbox-assets/releases/download/' releaseTag '/' assetName];
partPath=[localPath '.part'];
fprintf(['Downloading ' assetName '... ']);
try
    websave(partPath,url);
catch ME
    if isfile(partPath)
        delete(partPath);
    end
    if strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError')
        error('Failed to download %s, the server returned a 404 error, try again later. If the problem persists please report it on the RAVEN GitHub Issues page: https://github.com/CompMetEng/RAVEN/issues',assetName)
    else
        rethrow(ME);
    end
end
expected=manifestChecksum(url);
if isempty(expected)
    warning('RAVEN:fetchToolboxAsset:unverified',['The toolbox-assets manifest gives no ' ...
        'checksum for %s; the download is not verified.'],assetName);
else
    actual=sha256File(partPath);
    if ~strcmpi(actual,expected)
        delete(partPath);
        error('RAVEN:fetchToolboxAsset:checksum',['The download of %s is not the published ' ...
            'file (SHA-256 %s, expected %s) and was deleted. Try again; if it persists, ' ...
            'report it on the RAVEN GitHub Issues page: https://github.com/CompMetEng/RAVEN/issues'], ...
            assetName,actual,expected)
    end
end
movefile(partPath,localPath);
fprintf('COMPLETE\n');
end

function sha=manifestChecksum(url)
% The sha256 the toolbox-assets manifest lists for an asset URL, or '' when
% it lists none or cannot be read. The manifest is read once per session.
persistent manifest
sha='';
if isempty(manifest)
    try
        manifest=webread('https://github.com/CompMetEng/toolbox-assets/releases/download/manifest-v3/manifest.json', ...
            weboptions('ContentType','json','Timeout',30));
    catch ME
        warning('RAVEN:fetchToolboxAsset:manifest','Cannot read the toolbox-assets manifest: %s',ME.message);
        return
    end
end
sha=findUrl(manifest,url);
end

function sha=findUrl(node,url)
% Depth-first search for a struct with this url and a sha256 field. The
% manifest's file names are not valid MATLAB field names, so the search
% goes by value.
sha='';
if isstruct(node)
    for k=1:numel(node)
        if isfield(node(k),'url') && isfield(node(k),'sha256') && strcmp(node(k).url,url)
            sha=node(k).sha256;
            return
        end
        f=fieldnames(node(k));
        for i=1:numel(f)
            sha=findUrl(node(k).(f{i}),url);
            if ~isempty(sha)
                return
            end
        end
    end
elseif iscell(node)
    for i=1:numel(node)
        sha=findUrl(node{i},url);
        if ~isempty(sha)
            return
        end
    end
end
end

function sha=sha256File(path)
% SHA-256 of a file as lower-case hex, read in chunks.
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'r');
cleanup=onCleanup(@() fclose(fid));
while true
    chunk=fread(fid,2^24,'*uint8');
    if isempty(chunk)
        break
    end
    md.update(typecast(chunk,'int8'));
end
digest=typecast(md.digest(),'uint8');
sha=lower(reshape(dec2hex(digest,2)',1,[]));
end
