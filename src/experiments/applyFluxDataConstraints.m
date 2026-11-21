function model = applyFluxDataConstraints(model, fluxData, bioRxnId, varargin)
% applyFluxDataConstraints  Constrain the growth rate and exchange fluxes of a model to one measured condition.
%
% Sets the bounds of the growth reaction to the measured growth rate and the bounds of each
% exchange reaction measured in the condition to the measured flux; unmeasured exchanges (NaN) are
% left unchanged. The result can be used for flux variability analysis, sampling or other
% analyses of that condition.
%
% looseStrictFlux sets how tightly an exchange is constrained. For a measured uptake of -10 and
% excretion of 5, the lower and upper bounds are:
%
%     looseStrictFlux     uptake -10       excretion 5
%     'loose'             -10 and 0        0 and 5
%     10 (percent)        -10.5 and -9.5   4.75 and 5.25
%
% A measured flux of -1000 or 1000 means unconstrained whatever looseStrictFlux is: bounds -1000
% and 0 (free uptake) or 0 and 1000 (free excretion).
%
% Parameters
% ----------
% model : struct
%     model with the exchange reactions of fluxData. Not changed; the constrained model is
%     returned.
% fluxData : struct
%     flux data, as returned by loadFluxData.
% bioRxnId : char
%     id of the growth reaction.
%
% Name-Value Arguments
% --------------------
% condition : double or char
%     index (from 1) or name of the condition in fluxData.conds (default 1).
% maxMinGrowth : char
%     'max' bounds the growth reaction to 0 and the growth rate; 'min' bounds it to the growth
%     rate and 1000, which suits minimising a resource such as the protein pool (default 'max').
% looseStrictFlux : char or double
%     'loose' sets one bound to 0 and the other to the measured flux, which allows zero flux. A
%     number p between 0 and 100 sets the bounds to the flux times (1 -/+ p/200), a band of p
%     percent in total, which forces flux through the reaction (default 'loose').
% cSource : char
%     id of a carbon-source uptake reaction to close (bounds 0 and 0) before the data is applied,
%     so that the data alone determines the carbon source. [] or '' closes none; a reaction that
%     is not in the model is ignored (default []).
%
% Returns
% -------
% model : struct
%     model with the constrained fluxes.
%
% Notes
% -----
% Errors: applyFluxDataConstraints:invalidMaxMinGrowth, :invalidLooseStrictFlux,
% :unknownCondition, :conditionOutOfRange and :missingReaction (an exchange reaction of fluxData
% or the growth reaction that is not in the model).
%
% Examples
% --------
%     model = applyFluxDataConstraints(model, fluxData, 'growth');
%     model = applyFluxDataConstraints(model, fluxData, 'growth', ...
%         'condition', 'glucose', 'looseStrictFlux', 10, 'cSource', 'EX_glc');
%
% See also
% --------
% loadFluxData, simulateConditions

p = parseRAVENargs(varargin, { ...
    'condition', 1; ...
    'maxMinGrowth', 'max'; ...
    'looseStrictFlux', 'loose'; ...
    'cSource', []});
condition = p.condition;
maxMinGrowth = p.maxMinGrowth;
looseStrictFlux = p.looseStrictFlux;
cSource = p.cSource;

if ~(ischar(maxMinGrowth) && any(strcmp(maxMinGrowth, {'max','min'})))
    error('applyFluxDataConstraints:invalidMaxMinGrowth', ...
        'maxMinGrowth must be ''max'' or ''min''.');
end
isLoose = ischar(looseStrictFlux) && strcmp(looseStrictFlux, 'loose');
if ~(isLoose || (isnumeric(looseStrictFlux) && isscalar(looseStrictFlux) ...
        && looseStrictFlux >= 0 && looseStrictFlux <= 100))
    error('applyFluxDataConstraints:invalidLooseStrictFlux', ...
        'looseStrictFlux must be ''loose'' or a number in [0, 100].');
end

if ischar(condition) || isstring(condition)
    condIdx = find(strcmp(fluxData.conds, char(condition)), 1);
    if isempty(condIdx)
        error('applyFluxDataConstraints:unknownCondition', ...
            'Condition %s is not in fluxData.conds.', char(condition));
    end
else
    if ~(isscalar(condition) && condition >= 1 && condition <= numel(fluxData.conds) ...
            && condition == round(condition))
        error('applyFluxDataConstraints:conditionOutOfRange', ...
            'Condition index %g is out of range (fluxData has %d condition(s)).', ...
            condition, numel(fluxData.conds));
    end
    condIdx = condition;
end

[found, exchIdx] = ismember(fluxData.exchRxnIDs, model.rxns);
if ~all(found)
    missing = fluxData.exchRxnIDs(~found);
    error('applyFluxDataConstraints:missingReaction', ...
        '%d exchange reaction(s) of fluxData are not in the model (examples: %s).', ...
        numel(missing), strjoin(missing(1:min(5, end)), ', '));
end
bioIdx = find(strcmp(model.rxns, bioRxnId));
if isempty(bioIdx)
    error('applyFluxDataConstraints:missingReaction', ...
        'The biomass reaction %s is not in the model.', bioRxnId);
end

if ~isempty(cSource)
    cIdx = find(strcmp(model.rxns, cSource));
    model.lb(cIdx) = 0;
    model.ub(cIdx) = 0;
end

growth = fluxData.grRate(condIdx);
if strcmp(maxMinGrowth, 'max')
    model.lb(bioIdx) = 0;
    model.ub(bioIdx) = growth;
else
    model.lb(bioIdx) = growth;
    model.ub(bioIdx) = 1000;
end

fluxes = fluxData.exchFluxes(condIdx, :);
for k = 1:numel(fluxes)
    flux = fluxes(k);
    if isnan(flux)
        continue
    end
    idx = exchIdx(k);
    if abs(flux) == 1000
        if flux == -1000
            model.lb(idx) = -1000;
            model.ub(idx) = 0;
        else
            model.lb(idx) = 0;
            model.ub(idx) = 1000;
        end
    elseif isLoose
        if flux < 0
            model.lb(idx) = flux;
            model.ub(idx) = 0;
        else
            model.lb(idx) = 0;
            model.ub(idx) = flux;
        end
    else
        lo = flux * (1 - looseStrictFlux / 200);
        hi = flux * (1 + looseStrictFlux / 200);
        model.lb(idx) = min(lo, hi);
        model.ub(idx) = max(lo, hi);
    end
end
end
