function [rmse, rmseList] = datasetRmse(fluxData, sims, constrain, carbonWeights, varargin)
% datasetRmse  Distance between the simulation of one data set and its measurements.
%
% The error of a flux condition (constrain true) is the root mean square of the differences
% between carbon-weighted simulated and measured values: the growth rate, every measured exchange
% flux and every zero-flux reaction (measured as 0). Zero-flux reactions simulated at exactly 0
% are left out, so that correct predictions do not dilute the error. The error of a maximum-growth
% condition (constrain false) is the absolute difference of the carbon-weighted growth rates. A
% condition that could not be solved gets the penalty. The distance is the mean of the condition
% errors, each multiplied by the condition weight if the data has weights (the mean is over the
% conditions, not divided by the sum of the weights). A condition with weight 0 has error NaN and
% is not part of the mean. Values are in Cmmol/gDCW/h (see computeCarbonWeights).
%
% Parameters
% ----------
% fluxData : struct
%     the data set that was simulated, as returned by loadFluxData.
% sims : struct
%     result of simulateConditions for that data set.
% constrain : logical
%     the value given to simulateConditions: true compares growth and exchange fluxes, false
%     growth only.
% carbonWeights : containers.Map
%     weights from computeCarbonWeights. A reaction that is not in the map gets weight 1.
%
% Name-Value Arguments
% --------------------
% penalty : double
%     error of a condition that could not be solved (default 99), large enough that a model that
%     cannot grow scores worse than one that grows with a large error.
%
% Returns
% -------
% rmse : double
%     weighted mean of the errors of the conditions with weight other than 0 (Cmmol/gDCW/h). NaN
%     if every condition has weight 0.
% rmseList : double
%     error of each condition before weighting, with the penalty for unsolved conditions; NaN for
%     weight 0.
%
% Notes
% -----
% datasetRmse:wrongSimCount is raised if sims and fluxData have different numbers of conditions.
%
% Examples
% --------
% One maximum-growth condition with measured growth 0.4 1/h and simulated growth 0.5 1/h; the
% error is 0.1 * 41:
%
%     fluxData.conds = {'glc'}; fluxData.grRate = 0.4;
%     sims = struct('feasible', true, 'growth', 0.5);
%     [rmse, rmseList] = datasetRmse(fluxData, sims, false, containers.Map())   % 4.1
%
% See also
% --------
% simulateConditions, computeCarbonWeights, experimentDistance

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
