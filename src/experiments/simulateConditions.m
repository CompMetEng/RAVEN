function sims = simulateConditions(model, fluxData, constrain, zeroFluxRxns, bioRxnId, varargin)
% simulateConditions  Simulate every condition of a data set in a model.
%
% For each condition the carbon-source uptake is set, the model is solved, and the growth rate
% and the fluxes covered by the measurements are recorded. datasetRmse and experimentDistance
% turn the result into a distance. The model is not changed.
%
% The model objective must be the maximisation of growth. The function does not set it: each
% condition is solved with model.c, and the growth rate is the flux of bioRxnId in that solution.
% Set it first with setParam(model, 'obj', bioRxnId, 1). Otherwise the growth in a
% maximum-growth condition is not the maximal growth rate. A warning
% (simulateConditions:objectiveNotGrowth) is issued once per call if the objective is not
% bioRxnId alone with a positive coefficient, also when optimize is given, since the check reads
% the model.
%
% fluxData is either the flux conditions of experimentData.fluxData (constrain true) or the
% maximum-growth conditions of experimentData.maxGrowth (constrain false). The carbon source of a
% condition is the exchange reaction whose metabolite name (exchMets) is the condition name. For
% each condition:
%
% 1. The lower bound of every carbon-source exchange of a condition with weight other than 0 is
%    set to 0.
% 2. The lower bound of this condition's carbon source is set to the measured flux (constrain
%    true) or to -1000 (constrain false).
% 3. If the measured oxygen flux is exactly 0 and makeAnaerobic is given, the model is made
%    anaerobic. If Ptot is given and not 0 and changeProteinBiomass is given, the protein content
%    is changed.
% 4. The model is solved. With constrain true, the fluxes of the measured exchanges and of
%    zeroFluxRxns (except this carbon source) are recorded.
%
% A condition with weight 0 is skipped: it is not matched to an exchange reaction or simulated,
% and its element of sims has skipped true and feasible false.
%
% Parameters
% ----------
% model : struct
%     model with all other constraints already applied.
% fluxData : struct
%     data set to simulate, as returned by loadFluxData.
% constrain : logical
%     true for flux conditions (uptake limited to the measured rate, fluxes recorded), false for
%     maximum-growth conditions (uptake unlimited, only growth recorded).
% zeroFluxRxns : cell
%     ids of reactions that carry no flux in any condition (experimentData.zeroFlux).
% bioRxnId : char
%     id of the growth reaction, whose flux is the growth rate.
%
% Name-Value Arguments
% --------------------
% makeAnaerobic : function_handle
%     takes a model and returns the anaerobic model. No default, since this depends on the
%     organism; without it the oxygen flux of the data is ignored (default []).
% changeProteinBiomass : function_handle
%     takes a model and Ptot (g/gDCW) and returns the model with that protein content. Without it
%     Ptot is ignored (default []).
% optimize : function_handle
%     solves a model and returns a structure with x, f and stat, as solveLP does (default
%     @solveLP).
%
% Returns
% -------
% sims : struct
%     structure array with one element per condition and the fields:
%
%     - feasible : false if the solution is infeasible or not optimal, or the condition was
%       skipped. The other fields have no meaning then.
%     - skipped : true for a condition with weight 0.
%     - growth : flux of bioRxnId (1/h).
%     - exchRxnIDs : measured exchange reactions (empty if constrain is false).
%     - exchFluxes : their simulated fluxes (mmol/gDCW/h).
%     - blockRxnIDs : recorded zero-flux reactions.
%     - blockFluxes : their simulated fluxes.
%
% Notes
% -----
% Errors: simulateConditions:unmatchedCondition (a condition with weight other than 0 whose name
% matches no exchange metabolite) and simulateConditions:missingReaction (a reaction of
% fluxData, zeroFluxRxns or bioRxnId that is not in the model).
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

active = activeConditions(fluxData);
condToExch = matchConditions(fluxData, active);
conditionRxns = unique(fluxData.exchRxnIDs(unique(condToExch(active))));
allRxns = [conditionRxns(:); fluxData.exchRxnIDs(:); zeroFluxRxns(:); {bioRxnId}];
missing = unique(allRxns(~ismember(allRxns, model.rxns)));
if ~isempty(missing)
    error('simulateConditions:missingReaction', ...
        'The reactions %s are not in the model.', strjoin(missing(:)', ', '));
end
o2Col = find(strcmp(fluxData.exchMets, 'oxygen'), 1);
bioIdx = find(strcmp(model.rxns, bioRxnId));
if ~isfield(model, 'c') || nnz(model.c) ~= 1 || model.c(bioIdx) <= 0
    warning('simulateConditions:objectiveNotGrowth', ...
        ['The objective of the model is not the maximisation of the growth reaction %s. ' ...
        'The growth rate of a condition is the flux of %s in a solution of the objective ' ...
        'that the model has, and that is not the maximal growth rate. Set the objective with ' ...
        'setParam(model, ''obj'', ''%s'', 1).'], bioRxnId, bioRxnId, bioRxnId);
end
conditionIdx = find(ismember(model.rxns, conditionRxns));

nConds = numel(fluxData.conds);
template = struct('feasible', false, 'skipped', false, 'growth', NaN, ...
    'exchRxnIDs', {{}}, 'exchFluxes', zeros(0, 1), 'blockRxnIDs', {{}}, ...
    'blockFluxes', zeros(0, 1));
sims = repmat(template, nConds, 1);

for i = 1:nConds
    if ~active(i)
        sims(i).skipped = true;
        continue
    end
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

function active = activeConditions(fluxData)
% Conditions with a weight other than 0.
n = numel(fluxData.conds);
if isfield(fluxData, 'weight') && ~isempty(fluxData.weight)
    active = fluxData.weight(:) ~= 0;
else
    active = true(n, 1);
end
end

function condToExch = matchConditions(fluxData, active)
% Column of exchMets for the name of each active condition (0 for the others).
[found, condToExch] = ismember(fluxData.conds(:), fluxData.exchMets);
missing = unique(fluxData.conds(~found & active));
if ~isempty(missing)
    error('simulateConditions:unmatchedCondition', ...
        ['Carbon source(s) "%s" in the provided fluxData or maxGrowth cannot be ' ...
        'matched by name with an exchange reaction.'], strjoin(missing(:)', '; '));
end
end
