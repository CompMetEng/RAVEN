function [rmse, detail] = experimentDistance(experimentData, fluxSims, maxGrowthSims, carbonWeights, varargin)
% experimentDistance  Distance between a model and all measured conditions, as one number.
%
% Combines the distance of the flux conditions (experimentData.fluxData) and of the
% maximum-growth conditions (experimentData.maxGrowth), each from datasetRmse, as
% (dFlux + w * dMaxGrowth) / (1 + w) with w = maxGrowthWeight. A smaller distance is a better fit.
% Use it to compare model versions on the same data or as the objective when fitting parameters.
% Fluxes are compared in Cmmol/gDCW/h (computeCarbonWeights). A data set that is absent from
% experimentData, or in which every condition has weight 0, does not contribute.
%
% The maximum-growth comparison is the maximal growth rate only if the simulated model has growth
% as objective; simulateConditions warns if it does not.
%
% Parameters
% ----------
% experimentData : struct
%     measured conditions, as returned by loadExperimentData.
% fluxSims : struct
%     simulateConditions result for experimentData.fluxData (constrain true); [] if there is no
%     fluxData.
% maxGrowthSims : struct
%     simulateConditions result for experimentData.maxGrowth (constrain false); [] if there is no
%     maxGrowth.
% carbonWeights : containers.Map
%     weights from computeCarbonWeights.
%
% Name-Value Arguments
% --------------------
% penalty : double
%     error of a condition that could not be solved (default 99).
% maxGrowthWeight : double
%     weight of the maximum-growth distance against the flux distance (default 1, the average).
%
% Returns
% -------
% rmse : double
%     combined distance (Cmmol/gDCW/h). NaN if no condition has a weight other than 0.
% detail : struct
%     per-condition errors (rmseList of datasetRmse) in the fields fluxData and maxGrowth, for
%     the data sets that were scored.
%
% Notes
% -----
% experimentDistance:missingSims is raised if a data set has no simulation result.
%
% Examples
% --------
% One flux condition with measured growth 0.9 1/h and glucose uptake 10 mmol/gDCW/h, and a model
% that grows at 1.0 1/h at the same uptake:
%
%     fluxData.conds = {'glc'}; fluxData.Ptot = NaN; fluxData.grRate = 0.9;
%     fluxData.exchFluxes = -10; fluxData.exchMets = {'glc'}; fluxData.exchRxnIDs = {'EX_glc'};
%     experimentData = struct('fluxData', fluxData, 'maxGrowth', [], 'zeroFlux', {{}});
%     sims = struct('feasible', true, 'skipped', false, 'growth', 1.0, ...
%         'exchRxnIDs', {{'EX_glc'}}, 'exchFluxes', -10, 'blockRxnIDs', {{}}, ...
%         'blockFluxes', zeros(0, 1));
%     carbonWeights = containers.Map({'EX_glc', 'growth'}, [6 41]);
%     rmse = experimentDistance(experimentData, sims, [], carbonWeights)   % 2.899
%
% See also
% --------
% datasetRmse, simulateConditions, loadExperimentData

p = parseRAVENargs(varargin, {'penalty', 99; 'maxGrowthWeight', 1});
penalty = p.penalty;
maxGrowthWeight = p.maxGrowthWeight;

parts = [];
weights = [];
detail = struct();
if ~isempty(experimentData.fluxData)
    if isempty(fluxSims)
        error('experimentDistance:missingSims', ...
            'experimentData.fluxData is set but fluxSims is empty.');
    end
    [r, perCond] = datasetRmse(experimentData.fluxData, fluxSims, true, ...
        carbonWeights, 'penalty', penalty);
    if ~isnan(r)
        parts(end+1) = r;
        weights(end+1) = 1;
    end
    detail.fluxData = perCond;
end
if ~isempty(experimentData.maxGrowth)
    if isempty(maxGrowthSims)
        error('experimentDistance:missingSims', ...
            'experimentData.maxGrowth is set but maxGrowthSims is empty.');
    end
    [r, perCond] = datasetRmse(experimentData.maxGrowth, maxGrowthSims, false, ...
        carbonWeights, 'penalty', penalty);
    if ~isnan(r)
        parts(end+1) = r;
        weights(end+1) = maxGrowthWeight;
    end
    detail.maxGrowth = perCond;
end
if isempty(parts)
    rmse = NaN;
else
    rmse = sum(parts .* weights) / sum(weights);
end
end
