function setRavenProgress(mode)
% setRavenProgress  Set the RAVEN progress-reporting backend.
%
% Sets the preference that progressReport uses to decide how progress is
% shown. The setting is stored as the RAVEN preference progressBar and
% persists across MATLAB sessions.
%
% Parameters
% ----------
% mode : char
%     one of:
%
%     - 'auto' : choose automatically (default): 'gui' on the MATLAB desktop,
%       'cli' in a terminal, 'log' on a -batch or deployed session, and no
%       output during unit tests.
%     - 'gui'  : always use a MATLAB waitbar.
%     - 'cli'  : always use the animated command-line bar.
%     - 'log'  : always print plain milestone lines, which are readable in
%       redirected logs and on headless servers.
%     - 'none' or 'off' : never show progress.
%
% Examples
% --------
%     setRavenProgress('log');   % for example on a compute cluster
%     setRavenProgress('off');   % silence all progress reporting
%     setRavenProgress('auto');  % restore the default behaviour
%
% See also
% --------
% progressReport

if nargin < 1 || isempty(mode)
    mode = 'auto';
end
mode = lower(char(mode));
if strcmp(mode,'off')
    mode = 'none';
end
if ~ismember(mode,{'auto','gui','cli','log','none'})
    error('Unknown progress mode "%s". Use auto, gui, cli, log or off.', mode);
end
setpref('RAVEN','progressBar',mode);
end
