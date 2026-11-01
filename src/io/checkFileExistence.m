function files=checkFileExistence(files,varargin)
% checkFileExistence  Check whether files exist.
%
% A file name without a full path refers to a file in the current folder,
% which is prepended to the file name.
%
% Parameters
% ----------
% files : char or cell
%     string or cell array of strings with the path or file name of one or
%     more files.
%
% Name-Value Arguments
% --------------------
% fullOrTemp : double
%     controls path handling (default 0):
%
%     - 0 : do not change path to file(s)
%     - 1 : return full path to file(s)
%     - 2 : copy file(s) to system default temporary folder and return
%       full path
% allowSpace : logical
%     true if the space character is allowed in the path (default true).
% checkExist : logical
%     true if file existence is checked. Setting it to false returns the full
%     path to a file that does not exist yet (default true). Can only be set
%     to false if fullOrTemp is set to 1.
%
% Returns
% -------
% files : char or cell
%     string or cell array of strings with updated paths if fullOrTemp is 1
%     or 2, otherwise the original paths.
%
% Examples
% --------
%     files = checkFileExistence(files, fullOrTemp, allowSpace, checkExist);

p=parseRAVENargs(varargin, {'fullOrTemp',0; 'allowSpace',true; 'checkExist',true});
fullOrTemp=p.fullOrTemp; allowSpace=p.allowSpace; checkExist=p.checkExist;
files=convertCharArray(files);
if numel(files)==1
    oneFile=true;
else
    oneFile=false;
end
filesOriginal = files;

%Make all full paths before check of file existence
if ispc % full path starts like "C:\" or "C:/", or is a UNC path "\\server\share"
    inCurrDir = cellfun(@isempty,regexpi(files,'^([a-z]\:[\\\/]|\\\\)'));
else %isunix full path starts like "/"
    inCurrDir = cellfun(@isempty,regexpi(files,'^\/'));
end
files(inCurrDir) = fullfile(cd,files(inCurrDir));

%Check existence
if checkExist
    for i=1:numel(files)
        if ~isfile(files{i})
            error('File "%s" cannot be found\n',files{i});
        elseif allowSpace == false & strfind(files{i},' ')
            error('File "%s" has an invalid space in the filename or path, please remove this before running this function\n',files{i});
        end
    end
end

switch fullOrTemp
    case 0
        files = filesOriginal;
    case 1
        % files already contains full path
    case 2
        for i=1:numel(files)
            tmpFile=tempname;
            copyfile(files{i},tmpFile);
            files{i}=tmpFile;
        end
end

if oneFile == true
    files = files{1};
end