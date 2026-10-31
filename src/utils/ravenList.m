function msg = ravenList(header, items, trim)
% ravenList  Format a header and a list of items for warning() or error().
%
% Parameters
% ----------
% header : char
%     message header string.
% items : cell
%     cell array of items, listed one per line after the header.
% trim : logical
%     true to limit the list to 10 lines (default true). When the list has more
%     than 10 items, the tenth line reads '...and N more' and replaces items
%     10 to the end, with N the number of replaced items.
%
% Returns
% -------
% msg : char
%     header, a newline, and the items separated by newlines with a tab before
%     each item. If items is empty, msg is header. Passed to warning() or
%     error() as an argument after the format '%s'.
%
% Examples
% --------
%     warning("RAVEN:warning", "%s", ravenList("Missing metabolites:", badMets));
%     error("RAVEN:badInput", "%s", ravenList("Invalid reactions:", badRxns, false));

if nargin < 3; trim = true; end
items = convertCharArray(items);
if isempty(items)
    msg = header;
    return
end
if trim && numel(items) > 10
    items{10} = sprintf('...and %d more', numel(items)-9);
    items(11:end) = [];
end
% {char(9)} rather than char(9): strcat strips trailing whitespace from a
% plain char argument, and a lone tab is entirely trailing whitespace, so
% every item's indent was silently dropped. Wrapped in a cell, strcat
% treats it as a cell-array input instead, which it does not trim.
msg = [header, newline, strjoin(strcat({char(9)}, items), newline)];
end
