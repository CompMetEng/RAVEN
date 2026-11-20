function sims = simulateConditions(model, fluxData, constrain, zeroFluxRxns, bioRxnId, varargin)
% simulateConditions  Simulate every condition of a data set in a model.
%
% This is the simulation half of the distance between a model and measurements. For each condition
% the model is set to the uptake of the condition's carbon source, solved, and the growth rate and
% the fluxes that the measurements cover are recorded. datasetRmse and experimentDistance turn the
% result into a distance. The model is not changed.
%
% IMPORTANT: the model objective must be the maximisation of the growth. The function does not set
% the objective: each condition is solved with the objective that the model has (model.c, which
% solveLP maximises), and the growth rate is the flux of the biomass reaction bioRxnId in that
% solution. Set the objective with setParam(model, 'obj', bioRxnId, 1) before the simulation. With
% another objective the biomass flux is the value in one of the optimal solutions of that
% objective, not the maximal growth rate, and in a maximum-growth condition the distance then
% compares a number that is not the maximal growth rate with the measurement. A warning
% (simulateConditions:objectiveNotGrowth) is issued once per call if the objective of the model
% is not bioRxnId alone with a positive coefficient (another reaction, several reactions or no
% objective). The warning reads the model, so it is also issued when a function is given as
% optimize, which decides what is optimised: it solves the objective of the model unless it sets
% its own.
%
% A data set is either the flux conditions of experimentData.fluxData (constrain true) or the
% maximum-growth conditions of experimentData.maxGrowth (constrain false); see
% loadExperimentData.
%
% The carbon source of a condition is the exchange reaction whose metabolite name (exchMets) is the
% name of the condition.
%
% For each condition:
%
% 1. The lower bound of every exchange reaction that is the carbon source of some condition with
%    a weight other than 0 is set to 0, so that the carbon source of one condition cannot take up
%    substrate in another.
% 2. The lower bound of the carbon source of this condition is set to the measured flux if
%    constrain is true, so that the model can take up at most the measured rate, and to -1000 if
%    it is false (unlimited uptake, which gives the best possible growth on the carbon source).
% 3. If the measured flux of the exchange metabolite named 'oxygen' is exactly 0 and makeAnaerobic
%    is given, the model is made anaerobic. If Ptot of the condition is given and not 0 and
%    changeProteinBiomass is given, the protein content of the model is changed.
% 4. The model is solved. If constrain is true, the simulated flux of every measured exchange and
%    of every reaction in zeroFluxRxns (except the carbon source of the condition) is recorded.
%
% A condition with weight 0 (the field weight of fluxData) is left out: it is not matched to an
% exchange reaction and not simulated, and its element of sims has skipped true and feasible
% false. Its name need not match an exchange metabolite.
%
% Parameters
% ----------
% model : struct
%     model with all other constraints already applied.
% fluxData : struct
%     data set to simulate, as returned by loadFluxData.
% constrain : logical
%     true for flux conditions, which limit the uptake to the measured rate and record the
%     fluxes; false for maximum-growth conditions, which leave the uptake unlimited and record
%     only the growth rate.
% zeroFluxRxns : cell
%     ids of reactions that carry no flux in any condition (experimentData.zeroFlux). A flux
%     through them is counted as an error in the distance.
% bioRxnId : char
%     id of the biomass reaction, whose flux in the solution is the growth rate. It should be the
%     objective of the model (see above).
%
% Name-Value Arguments
% --------------------
% makeAnaerobic : function_handle
%     function that takes a model and returns the anaerobic model. How a model is made anaerobic
%     depends on the organism, so there is no default; without the function the oxygen flux of
%     the data is ignored (default []).
% changeProteinBiomass : function_handle
%     function that takes a model and Ptot (g/gDCW) and returns a model with that protein content.
%     Without the function Ptot is ignored (default []).
% optimize : function_handle
%     function that solves a model and returns a structure with the fields x, f and stat, as
%     solveLP does (default @solveLP).
%
% Returns
% -------
% sims : struct
%     structure array with one element per condition of fluxData and the fields:
%
%     - feasible : false if the model is infeasible or the solution not optimal, or if the
%       condition was skipped. The other fields have no meaning then.
%     - skipped : true for a condition with weight 0, which is not simulated.
%     - growth : flux of the biomass reaction (1/h).
%     - exchRxnIDs : ids of the exchange reactions that were measured in the condition (empty if
%       constrain is false).
%     - exchFluxes : simulated fluxes of exchRxnIDs (mmol/gDCW/h).
%     - blockRxnIDs : ids of the zero-flux reactions that were recorded.
%     - blockFluxes : simulated fluxes of blockRxnIDs.
%
% Notes
% -----
% A condition with a weight other than 0 that cannot be matched by name to an exchange metabolite
% of fluxData raises simulateConditions:unmatchedCondition. A reaction of fluxData, zeroFluxRxns
% or bioRxnId that is not in the model raises simulateConditions:missingReaction. An objective
% that is not the maximisation of bioRxnId issues simulateConditions:objectiveNotGrowth.
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
