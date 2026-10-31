function fractions = getBiomassFractions(model, biomassConfig)
% getBiomassFractions  Compute mass fraction per biomass component.
%
% Compute the mass fraction (g/gDW) per biomass component plus the total.
% The MATLAB counterpart of yeast-GEM's legacy sumBioMass.
%
% The biomassConfig struct describes the per-organism biomass layout (see
% Parameters). A component contributes 0 when its pseudoreaction is missing
% from the model or has no substrates. For a component whose pseudoreaction
% exists and has substrates, and whose mass_strategy is not 'grams', an error
% is raised when the molecular weight of a substrate is 0 (an empty
% metFormula field, or a formula that consists only of R), when its formula
% contains an element outside the built-in element table, or when
% mass_strategy is unknown. A component with an unknown mass_strategy whose
% pseudoreaction is missing or has no substrates contributes 0 without an
% error.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% biomassConfig : struct
%     Struct describing the biomass layout, with fields:
%
%     - biomass_rxn : rxn id of the top-level biomass pseudoreaction.
%     - proton_met : met id of cytosolic H+ (used by
%       scaleBiomassPseudoreaction, not by getBiomassFractions).
%     - components : cell array of component structs, each with fields:
%
%       - name : component name (for example 'protein').
%       - pseudoreaction_name : model.rxnNames entry to identify the
%         pseudoreaction.
%       - mass_strategy : 'mw', 'mw_minus_2h', 'mw_minus_water' or
%         'grams' (see Notes).
%
% Returns
% -------
% fractions : struct
%     Struct with one field per component name and a field 'total'
%     (for example fractions.protein and fractions.RNA). All values are in
%     g/gDW.
%
% Examples
% --------
%     fractions = getBiomassFractions(model, biomassConfig);
%
% Notes
% -----
% mass_strategy values:
%
%     'mw'              MW from chemical formula
%     'mw_minus_2h'     MW - 2.016 g/mol (two protons released per charged
%                       tRNA; protein pseudoreaction substrates)
%     'mw_minus_water'  MW - 18.015 g/mol (water released per
%                       polymerisation step; RNA and DNA)
%     'grams'           stoichiometry is already in g/gDW (lipid backbone)

fractions = struct();
total = 0;
for i = 1:numel(biomassConfig.components)
    comp = biomassConfig.components{i};
    f = computeComponentFraction(model, comp);
    fractions.(comp.name) = f;
    total = total + f;
end
fractions.total = total;
end

function f = computeComponentFraction(model, comp)
rxnPos = strcmp(model.rxnNames, comp.pseudoreaction_name);
if ~any(rxnPos)
    f = 0;
    return;
end
S_col = model.S(:, rxnPos);
isSub = find(S_col < 0);
if isempty(isSub)
    f = 0;
    return;
end
if strcmp(comp.mass_strategy, 'grams')
    f = full(-sum(S_col(isSub)));
    return;
end
offset = mwOffset(comp.mass_strategy);
formulas = model.metFormulas(isSub);
MWs = zeros(numel(formulas), 1);
for i = 1:numel(formulas)
    MWs(i) = computeFormulaMW(formulas{i});
end
zeroMW = MWs == 0;
if any(zeroMW)
    error('getBiomassFractions:emptyFormula', ...
        'Biomass metabolite %s has an empty metFormula field.', ...
        model.mets{isSub(find(zeroMW, 1))});
end
MWs = MWs + offset;
f = full(-sum(S_col(isSub) .* MWs) / 1000);
end

function offset = mwOffset(strategy)
switch strategy
    case 'mw'
        offset = 0;
    case 'mw_minus_2h'
        offset = -2.016;
    case 'mw_minus_water'
        offset = -18.015;
    otherwise
        error('getBiomassFractions:unknownStrategy', ...
            'Unknown mass_strategy: %s', strategy);
end
end

function mw = computeFormulaMW(formula)
% Molecular weight in g/mol from a Hill-style chemical formula. Same
% element table as the legacy yeast-GEM sumBioMass.
tokens = regexp(formula, '([A-Z][a-z]*)(\d*)', 'tokens');
if isempty(tokens)
    mw = 0;
    return;
end
tokensMatrix = vertcat(tokens{:});
tokensMatrix(cellfun(@isempty, tokensMatrix(:,2)), 2) = {'1'};
elements = tokensMatrix(:, 1);
counts = str2double(tokensMatrix(:, 2));
elem = {'C', 12.01; 'H', 1.008; 'N', 14.007; 'O', 15.999; ...
        'P', 30.974; 'S', 32.06; 'R', 0; ...
        'Fe', 55.845; 'K', 39.098; 'Na', 22.99; 'Cl', 35.45; ...
        'Mn', 54.938; 'Zn', 65.38; 'Ca', 40.078; 'Mg', 24.305; 'Cu', 63.546};
[~, elemMatch] = ismember(elements, elem(:,1));
if any(elemMatch == 0)
    error('getBiomassFractions:unknownElement', ...
        'Unknown element in formula %s', formula);
end
mw = sum(counts .* transpose([elem{elemMatch, 2}]), 'all');
end
