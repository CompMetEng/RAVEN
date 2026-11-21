function carbonWeights = computeCarbonWeights(model, rxnIds, bioRxnId)
% computeCarbonWeights  Weight of each reaction that puts its flux on a common carbon scale.
%
% Multiplying an exchange flux by the number of carbon atoms of the exchanged metabolite gives
% Cmmol/gDCW/h, so that fluxes of different metabolites can be added in one distance. The growth
% reaction gets 41 (assumed biomass carbon content in Cmmol/gDCW).
%
% Parameters
% ----------
% model : struct
%     model with metabolite formulas (metFormulas).
% rxnIds : cell
%     reactions to weight: the exchange reactions of the data sets, the zero-flux reactions and
%     the growth reaction.
% bioRxnId : char
%     id of the growth reaction.
%
% Returns
% -------
% carbonWeights : containers.Map
%     weight per reaction id. The growth reaction gets 41. Every other reaction gets the number
%     of carbon atoms in the formula of its single exchanged metabolite. A reaction with several
%     metabolites, a metabolite without a readable formula and a metabolite without carbon get 1,
%     so that mismatches in their fluxes still count.
%
%     Formulas are read by parseFormulas: element symbols with an optional integer or decimal
%     count (C6H12O6, C2.5H5); charge signs and * are ignored (C6H12O6*H2O). Parentheses
%     (C6H12O6(H2O)2, (C5H8)n), empty formulas and unknown elements are not readable.
%
% Notes
% -----
% A model without formulas gives weight 1 to every exchange and the warning
% computeCarbonWeights:noFormulas. A reaction that is not in the model raises
% computeCarbonWeights:missingReaction.
%
% Examples
% --------
%     carbonWeights = computeCarbonWeights(model, {'EX_glc','EX_o2','growth'}, 'growth');
%     carbonWeights('EX_glc')   % 6 for glucose (C6H12O6)
%
% See also
% --------
% datasetRmse, experimentDistance

rxnIds = cellstr(rxnIds);
biomassCarbon = 41;  % Cmmol/gDCW, assumed carbon content of biomass

if ~isfield(model, 'metFormulas') || all(cellfun(@isempty, model.metFormulas))
    warning('computeCarbonWeights:noFormulas', ...
        'The model has no metabolite formulas; every exchange reaction gets the weight 1.');
    formulas = repmat({''}, numel(model.mets), 1);
else
    formulas = model.metFormulas;
end

[found, rxnIdx] = ismember(rxnIds, model.rxns);
if ~all(found)
    error('computeCarbonWeights:missingReaction', ...
        'The reactions %s are not in the model.', strjoin(rxnIds(~found), ', '));
end

% The exchanged metabolite of every reaction that has exactly one, and the carbon atoms in the
% formulas of these metabolites (NaN when a formula cannot be read), parsed in one call.
metIdx = zeros(numel(rxnIds), 1);
for i = 1:numel(rxnIds)
    mets = find(model.S(:, rxnIdx(i)));
    if ~strcmp(rxnIds{i}, bioRxnId) && numel(mets) == 1
        metIdx(i) = mets;
    end
end
carbon = carbonAtoms(formulas(metIdx(metIdx > 0)));
carbonOfRxn = NaN(numel(rxnIds), 1);
carbonOfRxn(metIdx > 0) = carbon;

carbonWeights = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:numel(rxnIds)
    if strcmp(rxnIds{i}, bioRxnId)
        carbonWeights(rxnIds{i}) = biomassCarbon;
    elseif ~isnan(carbonOfRxn(i)) && carbonOfRxn(i) ~= 0
        carbonWeights(rxnIds{i}) = carbonOfRxn(i);
    else
        carbonWeights(rxnIds{i}) = 1;
    end
end
end

function carbon = carbonAtoms(formulas)
% Number of carbon atoms in each formula, NaN where the formula is empty or cannot be read.
% A * between parts of a formula is dropped, as cobra does, before the formulas are parsed.
carbon = NaN(numel(formulas), 1);
formulas = strtrim(strrep(formulas(:), '*', ''));
readable = ~cellfun(@isempty, formulas);
if ~any(readable)
    return
end
[elements, useMat, exitFlag] = parseFormulas(formulas(readable));
carbonCol = strcmp(elements.abbrevs, 'C');
count = NaN(nnz(readable), 1);
if any(carbonCol)
    count = useMat(:, carbonCol);
else
    count(:) = 0;
end
count(exitFlag ~= 1) = NaN;
carbon(readable) = count;
end
