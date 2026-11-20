function [rmse, rmseList] = datasetRmse(fluxData, sims, constrain, carbonWeights, varargin)
% datasetRmse  Distance between the simulation of one data set and its measurements.
%
% The error of a condition is the root mean square of the differences between the carbon-weighted
% simulated and measured values. For a flux condition (constrain true) the values are the growth
% rate, every measured exchange flux, and every reaction of the zero-flux list (measured as 0); a
% zero-flux reaction that is simulated at exactly 0 is left out so that correct predictions do not
% dilute the error. For a maximum-growth condition (constrain false) the error is the absolute
% difference of the carbon-weighted growth rates. A condition that the model could not solve gets
% a fixed penalty. The distance of the data set is the mean of the errors of its conditions, each
% multiplied by the weight of the condition if the data has weights (the mean is over the
% conditions, not divided by the sum of the weights). A condition with weight 0 is left out: its
% error is NaN and it is not part of the mean.
%
% The values are in Cmmol/gDCW/h; see computeCarbonWeights.
%
% Parameters
% ----------
% fluxData : struct
%     the data set that was simulated: experimentData.fluxData or experimentData.maxGrowth, as
%     returned by loadFluxData.
% sims : struct
%     simulation result of the data set from simulateConditions. The elements of the conditions
%     with weight 0 are not read.
% constrain : logical
%     must be the value that was given to simulateConditions: true compares growth and exchange
%     fluxes (flux conditions), false compares growth only (maximum-growth conditions).
% carbonWeights : containers.Map
%     weights from computeCarbonWeights. A reaction that is not in the map gets the weight 1.
%
% Name-Value Arguments
% --------------------
% penalty : double
%     error of a condition that could not be solved (default 99). It is much larger than the
%     error of a condition that is solved, so that a model that cannot grow in a condition has a
%     worse distance than a model that grows with a large error.
%
% Returns
% -------
% rmse : double
%     weighted mean of the errors of the conditions with a weight other than 0, in Cmmol/gDCW/h.
%     NaN if every condition has weight 0.
% rmseList : double
%     error of each condition, with the penalty for unsolved conditions and before the weighting;
%     NaN for a condition with weight 0.
%
% Notes
% -----
% A sims with another number of elements than fluxData has conditions raises
% datasetRmse:wrongSimCount.
%
% Examples
% --------
% One maximum-growth condition in which the organism grew at 0.4 1/h and the model at 0.5 1/h; the
% error is 0.1 * 41:
%
%     fluxData.conds = {'glc'}; fluxData.grRate = 0.4;
%     sims = struct('feasible', true, 'growth', 0.5);
%     [rmse, rmseList] = datasetRmse(fluxData, sims, false, containers.Map())   % 4.1
%
% See also
% --------
% simulateConditions : produces sims.
% computeCarbonWeights : produces carbonWeights.
% experimentDistance : combines the distances of the flux and the maximum-growth conditions.

p = parseRAVENargs(varargin, {'penalty', 99});
penalty = p.penalty;
biomassCarbon = 41;  % Cmmol/gDCW, assumed carbon content of biomass

n = numel(fluxData.conds);
if numel(sims) ~= n
    error('datasetRmse:wrongSimCount', ...
        'sims has %d elements; expected %d to match fluxData.conds.', numel(sims), n);
end

if isfield(fluxData, 'weight') && ~isempty(fluxData.weight)
    active = fluxData.weight(:) ~= 0;
else
    active = true(n, 1);
end

rmseList = zeros(n, 1);
for i = 1:n
    if ~active(i)
        rmseList(i) = NaN;
        continue
    end
    sim = sims(i);
    if ~sim.feasible
        rmseList(i) = penalty;
        continue
    end
    bioMeas = fluxData.grRate(i) * biomassCarbon;
    bioSim = sim.growth * biomassCarbon;
    if ~constrain
        rmseList(i) = abs(bioMeas - bioSim);
        continue
    end

    measured = bioMeas;
    simulated = bioSim;
    for k = 1:numel(sim.exchRxnIDs)
        w = weightOf(carbonWeights, sim.exchRxnIDs{k});
        col = find(strcmp(fluxData.exchRxnIDs, sim.exchRxnIDs{k}), 1);
        measured(end+1) = w * fluxData.exchFluxes(i, col); %#ok<AGROW>
        simulated(end+1) = w * sim.exchFluxes(k); %#ok<AGROW>
    end
    for k = 1:numel(sim.blockRxnIDs)
        weightedSim = weightOf(carbonWeights, sim.blockRxnIDs{k}) * sim.blockFluxes(k);
        if weightedSim == 0
            continue
        end
        measured(end+1) = 0; %#ok<AGROW>
        simulated(end+1) = weightedSim; %#ok<AGROW>
    end
    rmseList(i) = sqrt(mean((measured - simulated).^2));
end

if ~any(active)
    rmse = NaN;
    return
end
weighted = rmseList;
if isfield(fluxData, 'weight') && ~isempty(fluxData.weight)
    weighted = rmseList .* fluxData.weight(:);
end
rmse = mean(weighted(active));
end

function w = weightOf(carbonWeights, rxnId)
if isKey(carbonWeights, rxnId)
    w = carbonWeights(rxnId);
else
    w = 1;
end
end
