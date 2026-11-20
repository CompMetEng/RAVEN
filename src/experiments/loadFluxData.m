function fluxData = loadFluxData(fluxDataFile)
% loadFluxData  Read the measured growth rates and exchange fluxes of several conditions from a file.
%
% The file is tab-delimited with one header row and one row per condition. The first three
% columns are the name of the condition, the total protein content Ptot (g/gDCW) and the growth
% rate grRate (1/h); every further column is an exchange flux (mmol/gDCW/h, negative for uptake)
% with a header of the form '<metabolite name> (<reaction id>)':
%
%     Condition   Ptot   grRate   glucose (EX_glc)   ethanol (EX_etoh)   weight   source
%     glucose     0.40   0.40     -10.0              2.5                 1        study A
%     ethanol     0.35   0.12     NaN                -8.0                1        study B
%
% An empty or non-numeric cell is read as NaN, which means that the flux was not measured. A flux
% of -1000 or 1000 means unconstrained uptake or excretion.
%
% The columns weight and source are optional and may stand in any position. weight is the weight
% of the condition in the distance to the measurements, and 0 leaves the condition out of the
% simulation and of the distance (see simulateConditions). source is free text. Both are taken
% out of the flux matrix.
%
% The header of an exchange column is split at the parenthesis that matches the final closing
% one, so a reaction id may contain parentheses ('X (EX_a(e))' gives the metabolite 'X' and the
% reaction 'EX_a(e)'). A header without that form is used as both metabolite name and reaction
% id, with a warning.
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
%     - conds : cell column with the names of the conditions. The carbon source of a condition
%       is the exchange metabolite (exchMets) that has the name of the condition.
%     - Ptot : total protein content per condition (g/gDCW). NaN if not measured.
%     - grRate : measured growth rate per condition (1/h).
%     - exchFluxes : measured exchange fluxes, conditions in rows and exchange reactions in
%       columns (mmol/gDCW/h). NaN means that the flux was not measured; -1000 and 1000 mean
%       unconstrained uptake and excretion.
%     - exchMets : names of the exchanged metabolites, matching the columns of exchFluxes.
%     - exchRxnIDs : exchange reaction ids, matching the columns of exchFluxes.
%     - weight : only if the file has the column weight. Weight of each condition in the
%       distance to the measurements; 0 leaves the condition out of the simulation and of the
%       distance.
%     - source : only if the file has the column source. Where the measurements of each
%       condition come from.
%
% Notes
% -----
% The function raises loadFluxData:fileNotFound when the file does not exist,
% loadFluxData:emptyFile when it has no content, loadFluxData:noHeader when the header has
% fewer than three columns, loadFluxData:raggedRow when a line has a different number of cells
% than the header, loadFluxData:noExchangeColumns when the file has no exchange columns, and
% loadFluxData:renamedColumn when the file has the column evotuneRMSEweight (the column is
% called weight).
%
% Examples
% --------
%     fluxData = loadFluxData('data/fluxData.tsv');
%
% See also
% --------
% loadExperimentData : reads this file and the other files of a folder.
% simulateConditions : simulates the conditions of the file in a model.
% applyFluxDataConstraints : constrains a model to one condition of the file.

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

if any(strcmp(header, 'evotuneRMSEweight'))
    error('loadFluxData:renamedColumn', ...
        '%s: the column ''evotuneRMSEweight'' is now called ''weight''.', fluxDataFile);
end

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
