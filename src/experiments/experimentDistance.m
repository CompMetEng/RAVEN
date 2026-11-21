function [rmse, detail] = experimentDistance(experimentData, fluxSims, maxGrowthSims, carbonWeights, varargin)
% experimentDistance  Distance between a model and all measured conditions, as one number.
%
% Combines the distance of the flux conditions (experimentData.fluxData) and of the
% maximum-growth conditions (experimentData.maxGrowth), each from datasetRmse, as
% (dFlux + w * dMaxGrowth) / (1 + w), where w is maxGrowthWeight. Use it to compare model versions
% with the same data, or as the objective when parameters of a model are adjusted to the data; a
% smaller distance is a better fit. Fluxes are compared in Cmmol/gDCW/h: each flux is multiplied by
% the number of carbon atoms of the metabolite it exchanges (computeCarbonWeights).
%
% In a flux condition the uptake of the carbon source is limited to the measured rate, and the
% growth rate and the measured exchange fluxes of the model are compared with the measurements. In
% a maximum-growth condition the uptake of the carbon source is unlimited, and only the maximal
% growth rate of the model is compared with the measured growth rate; this is the maximal growth
% rate if the objective of the simulated model is the growth; simulateConditions warns if it is
% not. See
% loadExperimentData for the files and the workflow.
%
% A data set that is not in experimentData does not contribute, and neither does a data set in
% which every condition has weight 0.
%
% Parameters
% ----------
% experimentData : struct
%     the measured conditions, as returned by loadExperimentData.
% fluxSims : struct
%     simulation result of experimentData.fluxData from simulateConditions with constrain true.
%     [] if experimentData has no fluxData.
% maxGrowthSims : struct
%     simulation result of experimentData.maxGrowth from simulateConditions with constrain false.
%     [] if experimentData has no maxGrowth.
% carbonWeights : containers.Map
%     weights from computeCarbonWeights.
%
% Name-Value Arguments
% --------------------
% penalty : double
%     error of a condition that could not be solved (default 99).
% maxGrowthWeight : double
%     weight of the maximum-growth conditions against the flux conditions in the combination
%     (default 1, which averages the two distances). With 2 the maximum-growth conditions count
%     twice as much as the flux conditions.
%
% Returns
% -------
% rmse : double
%     the combined distance in Cmmol/gDCW/h. NaN if experimentData has neither data set or no
%     condition of its data sets has a weight other than 0.
% detail : struct
%     error of each condition (see rmseList of datasetRmse), in the fields fluxData and maxGrowth
%     for the data sets that were scored.
%
% Notes
% -----
% A data set in experimentData without simulation result raises experimentDistance:missingSims.
%
% Examples
% --------
% One flux condition in which the organism grew at 0.9 1/h and took up glucose at 10 mmol/gDCW/h,
% and the model grew at 1.0 1/h at the same uptake:
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
% datasetRmse : the distance of one data set.
% simulateConditions : produces the simulations.
% loadExperimentData : produces experimentData.

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
