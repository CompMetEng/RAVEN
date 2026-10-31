function args = parseRAVENargs(rawArgs, spec)
% parseRAVENargs  Resolve optional arguments passed positionally or by name.
%
% Helper for RAVEN functions that accept optional arguments in positional
% order, as name-value pairs, or as leading positional values followed by
% name-value pairs. The required leading arguments are explicit in the function
% signature. The optional arguments are collected in varargin and passed to
% this function together with a specification of the optional parameters.
%
% Parameters
% ----------
% rawArgs : cell
%     the varargin of the function, holding the arguments supplied after the
%     required positional arguments.
% spec : cell
%     an N-by-2 or N-by-3 cell array describing the optional parameters, one
%     per row:
%
%     - column 1 : parameter name (char).
%     - column 2 : default value, used when the parameter is not supplied.
%     - column 3 : optional validator, a function handle that is called with
%       the resolved value (the default when the parameter was not supplied)
%       and should error on an invalid value, for example @mustBeNumeric. [] means no
%       validation.
%
% Returns
% -------
% args : struct
%     a struct with one field per parameter name, holding the resolved value.
%
% Notes
% -----
% The calling convention is resolved in these steps:
%
%   1. If rawArgs is empty, all parameters take their defaults.
%   2. The first element in rawArgs that is a string matching a parameter
%      name (case-insensitive) marks position k. The elements before k are
%      positional and are assigned in the order of spec. The elements from k
%      onward must form valid name-value pairs: an even count in which every
%      name-position element is a parameter name.
%   3. If no element is a matching string, or the elements from k do not form
%      valid name-value pairs, all of rawArgs is positional.
%
% Three styles are accepted:
%
%     f(model, v1, v2)                   % purely positional
%     f(model, "p1", v1, "p2", v2)       % purely named
%     f(model, v1, "p2", v2, "p3", v3)   % hybrid: v1 positional, rest named
%
% A positional value that is a string equal to a parameter name is read as a
% parameter name only when the arguments from there on form valid name-value
% pairs; otherwise all arguments are positional. The name-value form passes
% such a value unambiguously.
%
% Examples
% --------
%     % inside exportToExcelFormat(model, varargin):
%     p = parseRAVENargs(varargin, { ...
%         "fileName", []; ...
%         "sortIds",  false});
%     fileName = p.fileName;
%     sortIds  = p.sortIds;
%
%     % the caller can use any of these:
%     exportToExcelFormat(model, "model.xlsx", true)              % positional
%     exportToExcelFormat(model, "sortIds", true)                 % named
%     exportToExcelFormat(model, "model.xlsx", "sortIds", true)   % hybrid

if nargin < 2 || isempty(spec)
    error('parseRAVENargs:noSpec', 'A parameter specification is required.');
end

names = spec(:, 1);
defaults = spec(:, 2);
hasValidators = size(spec, 2) >= 3;

% Seed every parameter with its default.
args = struct();
for i = 1:numel(names)
    args.(names{i}) = defaults{i};
end

if ~isempty(rawArgs)
    % Find which elements are strings matching a known parameter name.
    isParamName = cellfun(@(x) (ischar(x) || isstring(x)) && isscalar(string(x)) && ...
                          any(strcmpi(x, names)), rawArgs);

    % Position of the first parameter-name string — where name-value begins.
    nvStart = find(isParamName, 1);

    nvOK = false;
    if ~isempty(nvStart)
        nv = rawArgs(nvStart:end);
        % Valid name-value section: even count, all name-position elements
        % are known parameter names.
        nvOK = mod(numel(nv), 2) == 0 && ...
               all(cellfun(@(x) (ischar(x) || isstring(x)) && isscalar(string(x)) && ...
                   any(strcmpi(x, names)), nv(1:2:end)));
    end

    if ~isempty(nvStart) && nvOK
        % Assign leading positional arguments (may be zero when nvStart == 1).
        nPos = nvStart - 1;
        if nPos > numel(names)
            error('parseRAVENargs:tooManyArgs', ...
                'Too many positional arguments: expected at most %d, got %d.', ...
                numel(names), nPos);
        end
        for k = 1:nPos
            args.(names{k}) = rawArgs{k};
        end
        % Assign name-value pairs.
        nv = rawArgs(nvStart:end);
        for k = 1:2:numel(nv)
            args.(names{strcmpi(nv{k}, names)}) = nv{k + 1};
        end
    else
        % Purely positional.
        if numel(rawArgs) > numel(names)
            error('parseRAVENargs:tooManyArgs', ...
                'Too many positional arguments: expected at most %d, got %d.', ...
                numel(names), numel(rawArgs));
        end
        for k = 1:numel(rawArgs)
            args.(names{k}) = rawArgs{k};
        end
    end
end

% Apply optional validators.
if hasValidators
    for i = 1:numel(names)
        validator = spec{i, 3};
        if ~isempty(validator)
            validator(args.(names{i}));
        end
    end
end
end
