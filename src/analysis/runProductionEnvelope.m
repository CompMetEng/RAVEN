function [biomassValues, targetValues] = runProductionEnvelope(model, targetRxn, biomassRxn, varargin)
% runProductionEnvelope  Calculate the byproduct secretion envelope.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% targetRxn : char
%     identifier of target metabolite production reaction.
% biomassRxn : char
%     identifier of biomass reaction.
%
% Name-Value Arguments
% --------------------
% nPts : double
%     number of biomass values at which the target bounds are calculated
%     (default 20).
%
% Returns
% -------
% biomassValues : double
%     nPts biomass values, evenly spaced from the minimum to the maximum
%     biomass flux.
% targetValues : double
%     a matrix with two rows: the lower bounds (row 1) and the upper bounds
%     (row 2) of the target reaction flux at each biomass value. NaN where
%     the problem has no solution.
%
% Examples
% --------
%     [biomassValues, targetValues] = runProductionEnvelope(model, ...
%         targetRxn, biomassRxn, nPts);
%
% Notes
% -----
% Modified from COBRA Toolbox productionEnvelope.m.

p=parseRAVENargs(varargin, {'nPts',20});
nPts=p.nPts;

% Run FBA to get upper bound for biomass
model   = setParam(model,'obj',biomassRxn,1);
[solMax,hsSol]  = solveLP(model);
solMax  = solMax.x(logical(model.c));
model.c = -model.c;
solMin  = solveLP(model,0,[],hsSol);
solMin  = solMin.x(logical(model.c));

% Create biomass range vector
biomassValues = linspace(solMin,solMax,nPts);
targetUpperBound = nan(1,numel(biomassValues));
targetLowerBound = nan(1,numel(biomassValues));

PB = progressReport(length(biomassValues),'Creating production envelope');
% Max/min for target production
model = setParam(model,'obj',targetRxn,1);
for i = 1:length(biomassValues)
    model1 = setParam(model,'eq',biomassRxn,biomassValues(i));
    sol = solveLP(model1,0,[],hsSol);
    if (sol.stat > 0)
        targetUpperBound(i) = sol.x(logical(model.c));
    else
        targetUpperBound(i) = NaN;
    end
    model1.c = -model1.c;
    sol = solveLP(model1,0,[],hsSol);
    if (sol.stat > 0)
        targetLowerBound(i) = sol.x(logical(model1.c));
    else
        targetLowerBound(i) = NaN;
    end
    count(PB);
end

% Plot results
plot([biomassValues fliplr(biomassValues)],[targetUpperBound fliplr(targetLowerBound)],'blue','LineWidth',2);
axis tight;
ylabel([strrep(targetRxn,'_','-') ' (mmol/gDW h)']);
xlabel('Growth rate (1/h)');
targetValues = [targetLowerBound; targetUpperBound];
end
