function [rmse, rmseList] = datasetRmse(fluxData, sims, constrain, carbonWeights, bioRxnId, varargin)
% datasetRmse  Carbon-weighted distance between a simulated and a measured dataset.
%
% For each condition the root mean square error is computed over the
% carbon-weighted growth rate, and if constrain is true also over the
% carbon-weighted measured exchange fluxes and the zero-flux reactions. A
% zero-flux reaction that is simulated at exactly 0 is left out so that
% correct predictions do not dilute the error. An infeasible condition gets
% a fixed penalty. The distance of the dataset is the mean of the errors of
% its conditions, each multiplied by the weight of the condition if the data
% has one. A condition with weight 0 is left out: its error is NaN and it is
% not part of the mean.
%
% Parameters
% ----------
% fluxData : struct
%     dataset that was simulated, as returned by loadFluxData.
% sims : struct
%     simulation result of the dataset, from simulateConditions. The elements
%     of the conditions with weight 0 are not read.
% constrain : logical
%     must be the value that was given to simulateConditions: true includes
%     the exchange flux terms (flux data), false scores the growth alone
%     (maximum growth data).
% carbonWeights : containers.Map
%     weights from computeCarbonWeights. A reaction that is not in the map
%     gets the weight 1.
% bioRxnId : char
%     id of the biomass reaction. It is not read at the moment; the growth
%     rate of the simulation is taken from sims.
%
% Name-Value Arguments
% --------------------
% penalty : double
%     error of an infeasible condition (default 99). It must dominate the
%     errors of feasible conditions.
%
% Returns
% -------
% rmse : double
%     weighted mean of the errors of the conditions with a weight other than
%     0, NaN if every condition has weight 0.
% rmseList : double
%     error of each condition, with the penalty but before weighting; NaN for
%     a condition with weight 0.
%
% Notes
% -----
% A sims with another number of elements than fluxData has conditions raises
% datasetRmse:wrongSimCount.
%
% Examples
% --------
%     [rmse, rmseList] = datasetRmse(fluxData, sims, true, carbonWeights, 'growth');
%
% See also
% --------
% computeCarbonWeights, simulateConditions, experimentDistance

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
