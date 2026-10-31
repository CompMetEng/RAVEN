function legacyMethodNotice(oldFunction, preferred)
% legacyMethodNotice  Warn once that a function has a recommended replacement.
%
% For functions that remain supported and maintained but have a recommended
% replacement for new work. Unlike deprecationWarning, the notice announces no
% removal.
%
% The warning RAVEN:legacyMethod is issued at the first call of each
% oldFunction in a MATLAB session; later calls with the same oldFunction issue
% nothing. "clear functions" resets this. warning('off','RAVEN:legacyMethod')
% silences the notice.
%
% Parameters
% ----------
% oldFunction : char
%     name of the function being called. The notice is issued once per name.
% preferred : char
%     name of the recommended replacement for new work, written as in the
%     warning text.
%
% Examples
% --------
%     legacyMethodNotice('getINITModel','ftINIT');
%
% See also
% --------
% deprecationWarning

persistent alreadyNoticed

if isempty(alreadyNoticed)
    alreadyNoticed=containers.Map('KeyType','char','ValueType','logical');
end

oldFunction=char(oldFunction);
if isKey(alreadyNoticed,oldFunction)
    return
end
alreadyNoticed(oldFunction)=true;

warning('RAVEN:legacyMethod', ...
    ['%s is the legacy implementation and remains supported, but %s is the ' ...
     'recommended method for new work. Silence this with ' ...
     'warning(''off'',''RAVEN:legacyMethod'').'], oldFunction, char(preferred));
end
