function [minFluxes, maxFluxes, exitFlags]=getAllowedBounds(model,varargin)
% getAllowedBounds  Return the minimal and maximal fluxes through reactions.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% rxns : cell or logical or double
%     either a cell array of reaction IDs, a logical vector with the same
%     number of elements as reactions in the model, or a vector of
%     reaction indexes (default model.rxns).
% runParallel : logical
%     run the calculations in parallel (default true). Parallel processing
%     requires the MATLAB Parallel Computing Toolbox; without it, the
%     calculations run serially regardless of this setting. For a few
%     reactions, the overhead of parallel processing exceeds the time it
%     saves.
%
% Returns
% -------
% minFluxes : double
%     minimal allowed fluxes.
% maxFluxes : double
%     maximal allowed fluxes.
% exitFlags : double
%     a matrix with one row per reaction. Column 1 holds the solver status
%     of the minimisation and column 2 the status of the maximisation
%     (1 if a flux could be calculated).
%
% Notes
% -----
% In cases where no solution can be calculated, NaN is returned.
%
% Examples
% --------
%     [minFluxes, maxFluxes, exitFlags] = getAllowedBounds(model, rxns, runParallel);

p=parseRAVENargs(varargin, {'rxns',[]; 'runParallel',true});
rxns=p.rxns;
runParallel=p.runParallel;
if isempty(rxns)
    rxns = 1:numel(model.rxns);
elseif ~islogical(rxns) && ~isnumeric(rxns)
    rxns = convertCharArray(rxns);
    rxns = getIndexes(model,rxns, 'rxns');
end

nW = parallelWorkersRAVEN(runParallel);

minFluxes = zeros(numel(rxns),1);
maxFluxes = zeros(numel(rxns),1);
exitFlags = zeros(numel(rxns),2);
c = zeros(numel(model.rxns),1);

PB = progressReport(numel(rxns),'Running getAllowedBounds');
parfor (i = 1:numel(rxns), nW)
    count(PB)
    tmpModel = model;
    tmpModel.c = c;

    % Get maximal flux
    tmpModel.c(rxns(i)) = 1;
    [solMax,hsSol]=solveLP(tmpModel);
    if ~isempty(solMax.f)
        maxFluxes(i) = solMax.x(rxns(i));
    else
        maxFluxes(i) = NaN;
    end

    % Get minimal flux
    tmpModel.c(rxns(i)) = -1;
    solMin = solveLP(tmpModel,[],[],hsSol);
    if ~isempty(solMin.f)
        minFluxes(i) = solMin.x(rxns(i));
    else
        minFluxes(i) = NaN;
    end
    exitFlags(i,:) = [solMin.stat solMax.stat];

end
end
