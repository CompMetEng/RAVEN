function model=rescaleModelForINIT(model, varargin)
% rescaleModelForINIT  Rescale a model to avoid large stoichiometric differences.
%
% Rescales the stoichiometric matrix for the MILP problem in ftINIT to reduce
% differences in flux magnitudes between reactions. Such differences cause,
% among other effects, difficulties with the tolerances for integer variables.
% In each reaction, each coefficient with an absolute value above
% maxStoichDiff times the smallest absolute coefficient of the reaction is
% set to that limit (keeping its sign). The reaction is then multiplied by a
% factor that sets the mean absolute coefficient to 1.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% maxStoichDiff : double
%     largest allowed ratio between the largest and the smallest absolute
%     stoichiometric coefficient in a reaction (default 25).
%
% Returns
% -------
% model : struct
%     the rescaled model structure.

p=parseRAVENargs(varargin, {'maxStoichDiff',25});
maxStoichDiff=p.maxStoichDiff;

%Define maxMinRatio to be the ratio between the largest and smallest 
%stoichiometric coefficients in each reaction.
%Scale all rxns with maxMinRatio > maxStoichDiff - just set all coeffs that are
%larger than maxStoichDiff*minVal to that value, and center the mean coeff to 1 for all rxns.
SAbs = abs(model.S);
for i = 1:numel(model.rxns)
    tmp = SAbs(:,i);
    tmp = tmp(tmp ~= 0);
    mn = min(tmp);
    %modify matrix
    sign = ones(numel(model.mets),1);
    sign(model.S(:,i) < 0) = -1;
    sel = SAbs(:,i) > maxStoichDiff*mn;
    model.S(sel,i) = sign(sel) .* maxStoichDiff*mn;
end
%center around 1
boolMat = model.S ~= 0;
absMat = abs(model.S);
rxnScales = sum(boolMat,1)./sum(absMat,1);

model.S = model.S .* rxnScales;

%{
if (nargin < 2)
    maxStoichVal = 250;
end


%find all reactions containing high stoichiometry coefficients
maxCoeffs = max(model.S, [], 1);
rxnsToCheck = maxCoeffs > maxStoichVal;

%for debugging:
%constructEquations(model, model.rxns(rxnsToCheck))

rxnsToCheckInd = find(rxnsToCheck);

%We just scale all reactions so they don't have a coefficient > maxStoichVal for now

%do this with a loop, doesn't really matter
for rxnInd = rxnsToCheckInd
    %find the lowest index which is not a pointless metabolite, such as H+, H2O etc.
    rxnS = model.S(:,rxnInd);
    %metInd = find(rxnS ~= 0);
    
    largestCoeff = max(abs(rxnS));
    scaleFactor = maxStoichVal/largestCoeff;
    model.S(:,rxnInd) = model.S(:,rxnInd) .* scaleFactor;
end
%}
end
