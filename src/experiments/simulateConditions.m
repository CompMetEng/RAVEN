function sims = simulateConditions(model, fluxData, constrain, zeroFluxRxns, bioRxnId, varargin)
% simulateConditions  Simulate every condition of a dataset in a model.
%
% Performs the simulation half of the distance to measurements: for each
% condition of a flux dataset the uptake of the condition's carbon source is
% set, the model is solved, and the growth rate and the fluxes that are
% compared with the measurements are recorded. datasetRmse and
% experimentDistance turn the result into a distance.
%
% The steps for each condition are:
%
% 1. The lower bound of every exchange reaction that is the carbon source of
%    some condition is set to 0, so that the carbon source of one condition
%    cannot carry flux in another.
% 2. The lower bound of the carbon source of this condition is set to the
%    measured uptake rate if constrain is true, or to -1000 if it is false
%    (the best possible growth on the carbon source).
% 3. If the oxygen flux of the condition is exactly 0 and makeAnaerobic is
%    given, the model is made anaerobic. If Ptot of the condition is given
%    and not 0 and changeProteinBiomass is given, the protein content of the
%    model is changed.
% 4. The model is solved. If constrain is true, the simulated flux of every
%    measured exchange and of every reaction in zeroFluxRxns (except the
%    carbon source of the condition) is recorded.
%
% Parameters
% ----------
% model : struct
%     model with all other constraints already applied.
% fluxData : struct
%     dataset to simulate, as returned by loadFluxData. The carbon source of
%     a condition is the exchange metabolite with the name of the condition.
% constrain : logical
%     true for flux data, which fixes the uptake at the measured rate; false
%     for maximum growth data, which opens the uptake fully.
% zeroFluxRxns : cell
%     reaction ids that carry no flux in any condition.
% bioRxnId : char
%     id of the biomass reaction.
%
% Name-Value Arguments
% --------------------
% makeAnaerobic : function_handle
%     function that takes a model and returns the anaerobic model. These
%     steps depend on the organism, and there is no default; without the
%     function the anaerobic condition is not applied (default []).
% changeProteinBiomass : function_handle
%     function that takes a model and Ptot and returns a model with that
%     protein content (default []).
% optimize : function_handle
%     function that solves a model and returns a structure with the fields
%     x, f and stat, as solveLP does (default @solveLP).
%
% Returns
% -------
% sims : struct
%     structure array with one element per condition of fluxData and the
%     fields:
%
%     - feasible : false if the model is infeasible or the solution not
%       optimal. The other fields have no meaning then.
%     - growth : flux of the biomass reaction.
%     - exchRxnIDs : ids of the measured exchange reactions (empty if
%       constrain is false).
%     - exchFluxes : simulated fluxes of exchRxnIDs.
%     - blockRxnIDs : ids of the zero-flux reactions that were recorded.
%     - blockFluxes : simulated fluxes of blockRxnIDs.
%
% Notes
% -----
% A condition that cannot be matched by name to an exchange metabolite of
% fluxData raises simulateConditions:unmatchedCondition, and a reaction that
% is not in the model raises simulateConditions:missingReaction.
%
% Examples
% --------
%     sims = simulateConditions(model, experimentData.fluxData, true, ...
%         experimentData.zeroFlux, 'growth');
%
% See also
% --------
% loadExperimentData, datasetRmse, experimentDistance

p = parseRAVENargs(varargin, { ...
    'makeAnaerobic', []; ...
    'changeProteinBiomass', []; ...
    'optimize', []});
makeAnaerobic = p.makeAnaerobic;
changeProteinBiomass = p.changeProteinBiomass;
optimize = p.optimize;
if isempty(optimize)
    optimize = @(m) solveLP(m);
end
zeroFluxRxns = cellstr(zeroFluxRxns);

[condToExch, ~] = matchConditions(fluxData);
conditionRxns = unique(fluxData.exchRxnIDs(unique(condToExch)));
allRxns = [conditionRxns(:); fluxData.exchRxnIDs(:); zeroFluxRxns(:); {bioRxnId}];
missing = unique(allRxns(~ismember(allRxns, model.rxns)));
if ~isempty(missing)
    error('simulateConditions:missingReaction', ...
        'The reactions %s are not in the model.', strjoin(missing(:)', ', '));
end
o2Col = find(strcmp(fluxData.exchMets, 'oxygen'), 1);
bioIdx = find(strcmp(model.rxns, bioRxnId));
conditionIdx = find(ismember(model.rxns, conditionRxns));

nConds = numel(fluxData.conds);
template = struct('feasible', false, 'growth', NaN, 'exchRxnIDs', {{}}, ...
    'exchFluxes', zeros(0, 1), 'blockRxnIDs', {{}}, 'blockFluxes', zeros(0, 1));
sims = repmat(template, nConds, 1);

for i = 1:nConds
    carbonRxn = fluxData.exchRxnIDs{condToExch(i)};
    carbonIdx = find(strcmp(model.rxns, carbonRxn));
    m = model;
    m.lb(conditionIdx) = 0;
    if constrain
        m.lb(carbonIdx) = fluxData.exchFluxes(i, condToExch(i));
    else
        m.lb(carbonIdx) = -1000;
    end

    if ~isempty(o2Col) && fluxData.exchFluxes(i, o2Col) == 0 && ~isempty(makeAnaerobic)
        m = makeAnaerobic(m);
    end
    ptot = fluxData.Ptot(i);
    if ~(isnan(ptot) || ptot == 0) && ~isempty(changeProteinBiomass)
        m = changeProteinBiomass(m, ptot);
    end

    sol = optimize(m);
    [~, isOptimal] = checkSolution(sol);
    if ~isOptimal || isempty(sol.f) || isnan(sol.f)
        continue
    end

    sims(i).feasible = true;
    sims(i).growth = sol.x(bioIdx);
    if constrain
        measured = find(~isnan(fluxData.exchFluxes(i, :)));
        sims(i).exchRxnIDs = fluxData.exchRxnIDs(measured);
        [~, exchIdx] = ismember(sims(i).exchRxnIDs, m.rxns);
        sims(i).exchFluxes = sol.x(exchIdx);
        blockRxns = zeroFluxRxns(~strcmp(zeroFluxRxns, carbonRxn));
        [~, blockIdx] = ismember(blockRxns, m.rxns);
        sims(i).blockRxnIDs = blockRxns;
        sims(i).blockFluxes = sol.x(blockIdx);
    end
end
end

function [condToExch, missing] = matchConditions(fluxData)
% Column of exchMets for the name of each condition.
[found, condToExch] = ismember(fluxData.conds, fluxData.exchMets);
missing = unique(fluxData.conds(~found));
if ~isempty(missing)
    error('simulateConditions:unmatchedCondition', ...
        ['Carbon source(s) "%s" in the provided fluxData or maxGrowth cannot be ' ...
        'matched by name with an exchange reaction.'], strjoin(missing(:)', '; '));
end
end
