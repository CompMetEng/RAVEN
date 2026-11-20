function carbonWeights = computeCarbonWeights(model, rxnIds, bioRxnId)
% computeCarbonWeights  Weight of each reaction that puts its flux on a common carbon scale.
%
% Multiplying an exchange flux by the number of carbon atoms of the metabolite that it exchanges
% gives the flux in Cmmol/gDCW/h, so that fluxes of different metabolites can be compared and
% added in one distance. The biomass reaction gets 41, the assumed carbon content of biomass in
% Cmmol/gDCW (a carbon mass fraction of about 49%), which converts the growth rate to the same
% scale.
%
% Parameters
% ----------
% model : struct
%     model that provides the metabolite formulas (metFormulas).
% rxnIds : cell
%     ids of the reactions to compute a weight for: the exchange reactions of the data sets, the
%     reactions of the zero-flux list and the biomass reaction.
% bioRxnId : char
%     id of the biomass reaction.
%
% Returns
% -------
% carbonWeights : containers.Map
%     weight per reaction id. The biomass reaction gets 41. Every other reaction gets the number
%     of carbon atoms in the formula of its single exchanged metabolite. A reaction with more than
%     one metabolite, a metabolite without a formula, a formula that is not a plain list of
%     elements and counts, and a metabolite with no carbon atoms (O2, water) get the weight 1, so
%     that their flux mismatches still count in the distance.
%
% Notes
% -----
% A model without formulas gives every exchange the weight 1 and a warning
% (computeCarbonWeights:noFormulas). A reaction that is not in the model raises
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

carbonWeights = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:numel(rxnIds)
    if strcmp(rxnIds{i}, bioRxnId)
        carbonWeights(rxnIds{i}) = biomassCarbon;
        continue
    end
    mets = find(model.S(:, rxnIdx(i)));
    value = 1;
    if numel(mets) == 1
        carbon = carbonAtoms(formulas{mets});
        if ~isnan(carbon) && carbon ~= 0
            value = carbon;
        end
    end
    carbonWeights(rxnIds{i}) = value;
end
end

function carbon = carbonAtoms(formula)
% Number of carbon atoms in a formula, NaN when the formula has no carbon or
% cannot be read as a plain list of elements and counts.
carbon = NaN;
if isempty(formula)
    return
end
tokens = regexp(formula, '([A-Z][a-z]?)(\d*\.?\d*)', 'tokens');
if isempty(tokens) || ~strcmp(strjoin(cellfun(@(t) [t{1} t{2}], tokens, 'UniformOutput', false), ''), formula)
    return
end
for i = 1:numel(tokens)
    if strcmp(tokens{i}{1}, 'C')
        count = str2double(tokens{i}{2});
        if isnan(count)
            count = 1;
        end
        if isnan(carbon)
            carbon = 0;
        end
        carbon = carbon + count;
    end
end
end
