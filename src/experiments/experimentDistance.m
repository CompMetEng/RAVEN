function [rmse, detail] = experimentDistance(experimentData, fluxSims, maxGrowthSims, carbonWeights, bioRxnId, varargin)
% experimentDistance  Combine the distances of both datasets into one value.
%
% The distance of the flux dataset and of the maximum growth dataset (see
% datasetRmse) are combined as (rmseFlux + w * rmseMaxGrowth) / (1 + w),
% where w is maxGrowthWeight. A dataset that is not in experimentData does
% not contribute, and neither does a dataset in which every condition has
% weight 0.
%
% Parameters
% ----------
% experimentData : struct
%     measured conditions, as returned by loadExperimentData.
% fluxSims : struct
%     simulation result of experimentData.fluxData from simulateConditions
%     with constrain true. [] if experimentData has no fluxData.
% maxGrowthSims : struct
%     simulation result of experimentData.maxGrowth from simulateConditions
%     with constrain false. [] if experimentData has no maxGrowth.
% carbonWeights : containers.Map
%     weights from computeCarbonWeights.
% bioRxnId : char
%     id of the biomass reaction.
%
% Name-Value Arguments
% --------------------
% penalty : double
%     error of an infeasible condition (default 99).
% maxGrowthWeight : double
%     relative weight of the maximum growth dataset against the flux dataset
%     (default 1). With 2, the maximum growth dataset counts twice as much as
%     the flux dataset in the combination.
%
% Returns
% -------
% rmse : double
%     combined distance, NaN if experimentData has neither dataset or no
%     condition of its datasets has a weight other than 0.
% detail : struct
%     error of each condition, in the fields fluxData and maxGrowth for the
%     datasets that were scored.
%
% Notes
% -----
% A dataset in experimentData without simulation result raises
% experimentDistance:missingSims.
%
% Examples
% --------
%     [rmse, detail] = experimentDistance(experimentData, fluxSims, ...
%         maxGrowthSims, carbonWeights, 'growth', 'maxGrowthWeight', 2);
%
% See also
% --------
% datasetRmse, simulateConditions, computeCarbonWeights

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
        carbonWeights, bioRxnId, 'penalty', penalty);
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
        carbonWeights, bioRxnId, 'penalty', penalty);
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
