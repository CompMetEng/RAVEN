function fluxData = loadFluxData(fluxDataFile)
% loadFluxData  Read measured growth rates and exchange fluxes of several conditions from a file.
%
% The file is tab-delimited with one row per condition. The first three columns are the condition
% name, the total protein content Ptot (g/gDCW) and the growth rate grRate (1/h). Every further
% column is an exchange flux (mmol/gDCW/h, negative for uptake) headed '<metabolite name>
% (<reaction id>)':
%
%     Condition   Ptot   grRate   glucose (EX_glc)   ethanol (EX_etoh)   weight   source
%     glucose     0.40   0.40     -10.0              2.5                 1        study A
%     ethanol     0.35   0.12     NaN                -8.0                1        study B
%
% An empty or non-numeric cell is NaN (not measured). A flux of -1000 or 1000 means unconstrained
% uptake or excretion. The optional columns weight and source may stand in any position. weight
% is the weight of the condition in the distance; 0 leaves the condition out of the simulation
% and the distance (see simulateConditions). source is free text.
%
% The header is split at the parenthesis matching the final closing one, so reaction ids may
% contain parentheses ('X (EX_a(e))' gives 'X' and 'EX_a(e)'). A header without that form is used
% as both metabolite name and reaction id, with a warning (loadFluxData:headerForm).
%
% Parameters
% ----------
% fluxDataFile : char
%     path to the tab-delimited file.
%
% Returns
% -------
% fluxData : struct
%     structure with one row per condition and the fields:
%
%     - conds : cell column with the condition names. The carbon source of a condition is the
%       exchange metabolite with that name.
%     - Ptot : total protein content (g/gDCW), NaN if not measured.
%     - grRate : growth rate (1/h).
%     - exchFluxes : fluxes (mmol/gDCW/h), conditions in rows and exchange reactions in columns.
%     - exchMets : metabolite names of the columns of exchFluxes.
%     - exchRxnIDs : reaction ids of the columns of exchFluxes.
%     - weight : weight per condition; only if the file has the column.
%     - source : origin of each condition; only if the file has the column.
%
% Notes
% -----
% Errors: loadFluxData:fileNotFound, loadFluxData:emptyFile, loadFluxData:noHeader (fewer than
% three columns), loadFluxData:raggedRow (a line with another number of cells than the header) and
% loadFluxData:noExchangeColumns.
%
% Examples
% --------
%     fluxData = loadFluxData('data/fluxData.tsv');
%
% See also
% --------
% loadExperimentData, simulateConditions, applyFluxDataConstraints

if ~isfile(fluxDataFile)
    error('loadFluxData:fileNotFound', 'The flux data file %s does not exist.', fluxDataFile);
end
lines = splitlines(string(fileread(fluxDataFile)));
lines = erase(lines, char(13));
lines = lines(strlength(strtrim(lines)) > 0);
if isempty(lines)
    error('loadFluxData:emptyFile', 'The flux data file %s is empty.', fluxDataFile);
end

raw = cell(numel(lines), 1);
for i = 1:numel(lines)
    raw{i} = cellstr(strsplit(lines(i), char(9), 'CollapseDelimiters', false));
end
header = raw{1};
if numel(header) < 3
    error('loadFluxData:noHeader', ...
        '%s: the header has %d column(s); expected at least condition, Ptot and grRate.', ...
        fluxDataFile, numel(header));
end
for i = 2:numel(raw)
    if numel(raw{i}) ~= numel(header)
        error('loadFluxData:raggedRow', ...
            'Line %d of %s has %d cells, the header has %d.', ...
            i, fluxDataFile, numel(raw{i}), numel(header));
    end
end
cells = cell(numel(raw), numel(header));
for i = 1:numel(raw)
    cells(i, :) = raw{i};
end
header = strtrim(header);

fluxData = struct();
extra = false(1, numel(header));
weightIdx = find(strcmp(header, 'weight'), 1);
sourceIdx = find(strcmp(header, 'source'), 1);
if ~isempty(weightIdx)
    fluxData.weight = str2double(cells(2:end, weightIdx));
    extra(weightIdx) = true;
end
if ~isempty(sourceIdx)
    fluxData.source = cells(2:end, sourceIdx);
    extra(sourceIdx) = true;
end

kept = find(~extra);
if numel(kept) < 4
    error('loadFluxData:noExchangeColumns', ...
        '%s has no exchange flux columns after the columns condition, Ptot and grRate.', ...
        fluxDataFile);
end
exchHeaders = header(kept(4:end));
exchMets = exchHeaders;
exchRxns = exchHeaders;
noForm = {};
for i = 1:numel(exchHeaders)
    [exchMets{i}, exchRxns{i}, ok] = splitExchHeader(exchHeaders{i});
    if ~ok
        noForm{end+1} = exchHeaders{i}; %#ok<AGROW>
    end
end
if ~isempty(noForm)
    warning('loadFluxData:headerForm', ...
        ['%d exchange column header(s) lack the ''<met> (<rxn>)'' form; their ' ...
        'reaction id is the full header text: %s'], ...
        numel(noForm), strjoin(noForm(1:min(5, end)), ', '));
end

fluxData.conds      = cells(2:end, kept(1));
fluxData.Ptot       = str2double(cells(2:end, kept(2)));
fluxData.grRate     = str2double(cells(2:end, kept(3)));
fluxData.exchFluxes = str2double(cells(2:end, kept(4:end)));
fluxData.exchMets   = exchMets(:)';
fluxData.exchRxnIDs = exchRxns(:)';
% Order the fields as documented.
order = {'conds','Ptot','grRate','exchFluxes','exchMets','exchRxnIDs','weight','source'};
order = order(isfield(fluxData, order));
fluxData = orderfields(fluxData, order);
end

function [met, rxn, ok] = splitExchHeader(header)
% Split '<metabolite> (<reaction id>)' at the parenthesis that matches the
% final closing one. A header of another form gives the whole header.
header = strtrim(header);
met = header;
rxn = header;
ok = false;
if isempty(header) || header(end) ~= ')'
    return
end
depth = 0;
for k = numel(header):-1:1
    if header(k) == ')'
        depth = depth + 1;
    elseif header(k) == '('
        depth = depth - 1;
        if depth == 0
            candidate = strtrim(header(k+1:end-1));
            if ~isempty(candidate)
                met = strtrim(header(1:k-1));
                rxn = candidate;
                ok = true;
            end
            return
        end
    end
end
end
