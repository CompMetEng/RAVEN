classdef progressReport < handle
% progressReport  Universal progress reporting for RAVEN functions.
%
% Progress-reporting class used by RAVEN functions. The rendering backend is
% selected automatically. The same calls work on an interactive desktop, in a
% command-line session and on a headless server. During unit tests no output
% is produced. The backends are:
%
%   - 'gui'  : MATLAB waitbar (interactive desktop).
%   - 'cli'  : animated command-line bar (interactive terminal).
%   - 'log'  : plain milestone lines, for example "Running BLAST: 30%", printed
%              every 10%. The lines contain no control characters, so they
%              are readable in redirected logs and on headless and -batch
%              servers.
%   - 'none' : no output (during tests, or when disabled).
%
% setRavenProgress sets the backend: setRavenProgress('log') selects plain
% lines, for example on a compute cluster, and setRavenProgress('off')
% disables all progress reporting.
%
% Usage
% -----
% To register each iteration, call count once per pass:
%
%     PB = progressReport(N, 'Task name');   % before the loop
%     for i = 1:N
%         ...                                % loop body
%         PB.count;                          % or count(PB), also in parfor
%     end
%     PB.done;                               % optional; finalises the bar on early exit
%
% To avoid the per-iteration overhead in long loops of short iterations, set
% the absolute position at intervals:
%
%     PB = progressReport(N, 'Task name');
%     for i = 1:N
%         ...
%         if rem(i-1, 1000) == 0
%             PB.update(i);                  % or update(PB, i)
%         end
%     end
%     PB.done;
%
% Parameters
% ----------
% N : double
%     number of iterations the loop performs (a positive integer). For any other
%     value, no progress is shown.
% taskName : char, optional
%     label shown alongside the progress (default '').
%
% See also
% --------
% setRavenProgress

    properties (SetAccess = private)
        mode char            % resolved backend: 'gui', 'cli', 'log' or 'none'
        N    double
        taskName char
    end
    properties (Access = private, Transient)
        counter   = 0
        lastPct   = -1       % last percentage rendered
        logStep   = 10       % percent between 'log' lines
        barWidth  = 30
        prevLen   = 0        % characters on the current 'cli' line
        guiHandle = []
        Queue     = []       % DataQueue for parfor counting
        Listener  = []
        noPCT     = true
        finished  = false
    end

    methods
        function obj = progressReport(N, taskName)
            if nargin < 2 || isempty(taskName)
                taskName = '';
            end
            obj.N = double(N);
            obj.taskName = truncateName(char(taskName), 40);
            obj.mode = resolveProgressMode();

            % Nothing sensible to show for an empty or unknown-length loop.
            if ~(isscalar(obj.N) && isfinite(obj.N) && obj.N >= 1)
                obj.mode = 'none';
            end
            if strcmp(obj.mode,'none')
                return
            end

            obj.noPCT = isempty(ver('parallel'));
            if ~obj.noPCT
                obj.Queue = parallel.pool.DataQueue;
                obj.Listener = afterEach(obj.Queue, @(v) localUpdate(obj, v));
            end

            % Initial render.
            switch obj.mode
                case 'gui'
                    name = obj.taskName;
                    if isempty(name); name = 'Progress'; end
                    try
                        obj.guiHandle = waitbar(0, '0%', 'Name', name);
                    catch
                        obj.mode = 'log';   % no display available
                        emitLog(obj, 0);
                    end
                    obj.lastPct = 0;
                case 'cli'
                    drawCli(obj, 0);
                    obj.lastPct = 0;
                case 'log'
                    emitLog(obj, 0);
                    obj.lastPct = 0;
            end
        end

        function count(obj)
            % Register one completed iteration.
            if strcmp(obj.mode,'none')
                return
            end
            if obj.noPCT || isempty(getCurrentJob)
                localUpdate(obj, -1);       % serial execution
            else
                send(obj.Queue, -1);        % parfor worker -> client
            end
        end

        function update(obj, n)
            % Set the absolute number of completed iterations to n.
            if strcmp(obj.mode,'none')
                return
            end
            if obj.noPCT || isempty(getCurrentJob)
                localUpdate(obj, n);        % serial execution
            else
                send(obj.Queue, n);         % parfor worker -> client
            end
        end

        function done(obj)
            % Finalise the bar (call when a loop may exit before reaching N).
            finishUp(obj);
        end
    end

    methods (Hidden = true)
        function delete(obj)
            if ~isempty(obj.Queue)
                delete(obj.Queue);
            end
            finishUp(obj);
        end
    end

    methods (Access = private)
        function localUpdate(obj, v)
            if v < 0
                obj.counter = obj.counter + 1;
            else
                obj.counter = max(obj.counter, v);
            end
            render(obj);
        end

        function render(obj)
            if obj.finished
                return
            end
            pct = min(100, floor(100 * obj.counter / obj.N));
            switch obj.mode
                case 'gui'
                    if pct ~= obj.lastPct
                        obj.lastPct = pct;
                        if ~isempty(obj.guiHandle) && isvalid(obj.guiHandle)
                            waitbar(pct/100, obj.guiHandle, sprintf('%d%%', pct));
                        end
                    end
                    if obj.counter >= obj.N; finishUp(obj); end
                case 'cli'
                    if pct ~= obj.lastPct
                        obj.lastPct = pct;
                        drawCli(obj, pct);
                    end
                    if obj.counter >= obj.N; finishUp(obj); end
                case 'log'
                    if pct >= obj.lastPct + obj.logStep ...
                            || (obj.counter >= obj.N && obj.lastPct < 100)
                        if obj.counter >= obj.N; pct = 100; end
                        obj.lastPct = pct;
                        emitLog(obj, pct);
                    end
            end
        end

        function drawCli(obj, pct)
            nfill = round(obj.barWidth * pct / 100);
            bar = [repmat('#',1,nfill), repmat('-',1,obj.barWidth-nfill)];
            if isempty(obj.taskName)
                s = sprintf('[%s] %3d%%', bar, pct);
            else
                s = sprintf('%s [%s] %3d%%', obj.taskName, bar, pct);
            end
            fprintf('%s%s', repmat(sprintf('\b'),1,obj.prevLen), s);
            obj.prevLen = numel(s);
        end

        function emitLog(obj, pct)
            if isempty(obj.taskName)
                fprintf('Progress: %d%%\n', pct);
            else
                fprintf('%s: %d%%\n', obj.taskName, pct);
            end
        end

        function finishUp(obj)
            if obj.finished
                return
            end
            switch obj.mode
                case 'gui'
                    if ~isempty(obj.guiHandle) && isvalid(obj.guiHandle)
                        delete(obj.guiHandle);
                    end
                case 'cli'
                    if obj.lastPct < 100
                        drawCli(obj, 100);  % show completion
                    end
                    if obj.prevLen > 0
                        fprintf('\n');      % move off the progress line
                    end
                case 'log'
                    if obj.lastPct < 100
                        emitLog(obj, 100);
                    end
            end
            obj.finished = true;
        end
    end
end

% ---------------------------- Local functions ----------------------------
function mode = resolveProgressMode()
% Decide the backend from the RAVEN preference, the test context and the
% runtime environment.
mode = 'auto';
try
    mode = char(getpref('RAVEN','progressBar','auto'));
catch
end
if strcmpi(mode,'off')
    mode = 'none';
end
if ~strcmpi(mode,'auto')
    if ~ismember(lower(mode),{'gui','cli','log','none'})
        mode = 'cli';
    end
    mode = lower(mode);
    return
end

% 'auto': silence inside the test framework first.
if inTestContext()
    mode = 'none';
    return
end

% Headless / deployed / -batch -> clean log lines.
isBatch = false;
try
    isBatch = batchStartupOptionUsed;
catch
end
if isBatch || isdeployed
    mode = 'log';
    return
end

% Desktop available -> GUI, otherwise an interactive terminal -> CLI.
hasDesktop = false;
try
    hasDesktop = usejava('desktop');
catch
end
if hasDesktop
    mode = 'gui';
else
    mode = 'cli';
end
end

function tf = inTestContext()
% Best-effort detection of the matlab.unittest framework on the call stack.
% RavenTestCase also forces the preference to 'none', so the maintained
% suite is silenced even if this heuristic misses a frame.
tf = false;
st = dbstack('-completenames');
for i = 1:numel(st)
    if contains(st(i).name,'matlab.unittest','IgnoreCase',true) ...
            || contains(st(i).name,'TestRunner','IgnoreCase',true) ...
            || contains(st(i).name,'runtests','IgnoreCase',true)
        tf = true;
        return
    end
end
end

function name = truncateName(name, maxLen)
if numel(name) > maxLen
    name = name(1:maxLen);
end
end
