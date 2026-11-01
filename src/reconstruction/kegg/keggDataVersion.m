function version=keggDataVersion()
% keggDataVersion  The toolbox-assets release tag for KEGG artefacts.
%
% Reads the first line of src/reconstruction/kegg/KEGG_VERSION.md. The file is
% the single place where the version is recorded, so that every function that
% downloads a KEGG artefact (getModelFromKEGG, getKEGGModelForOrganism) uses
% the same release.
%
% Returns
% -------
% version : char
%     the toolbox-assets release tag, for example 'kegg118'.

ravenPath=findRAVENroot();
versionFile=fullfile(ravenPath,'src','reconstruction','kegg','KEGG_VERSION.md');
fid=fopen(versionFile,'r');
if fid==-1
    error('keggDataVersion:fileNotFound','Cannot read %s.',strrep(versionFile,'\','/'))
end
version=strtrim(fgetl(fid));
fclose(fid);
end
