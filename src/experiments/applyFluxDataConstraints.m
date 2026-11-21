function model = applyFluxDataConstraints(model, fluxData, bioRxnId, varargin)
% applyFluxDataConstraints  Constrain the growth rate and exchange fluxes of a model to one measured condition.
%
% Sets the bounds of the biomass reaction to the measured growth rate and the bounds of each
% exchange reaction of fluxData that was measured in the condition to the measured flux. An
% exchange that was not measured (NaN) is left as it is. The constrained model can then be used
% for flux variability analysis, sampling or other analyses of that condition.
%
% How tightly an exchange is constrained is set by looseStrictFlux. For a measured uptake of -10
% and a measured excretion of 5, the lower and upper bounds are:
%
%     looseStrictFlux     uptake -10       excretion 5
%     'loose'             -10 and 0        0 and 5
%     10 (percent)        -10.5 and -9.5   4.75 and 5.25
%
% A measured flux of -1000 or 1000 means unconstrained, whatever looseStrictFlux is: -1000 gives
% the bounds -1000 and 0 (free uptake), and 1000 gives 0 and 1000 (free excretion).
%
% Parameters
% ----------
% model : struct
%     model with the exchange reactions of fluxData. Not changed; the constrained model is
%     returned.
% fluxData : struct
%     flux data, as returned by loadFluxData.
% bioRxnId : char
%     id of the biomass reaction.
%
% Name-Value Arguments
% --------------------
% condition : double or char
%     index of the condition in fluxData.conds (counted from 1), or its name (default 1).
% maxMinGrowth : char
%     'max' sets the bounds of the biomass reaction to 0 and the growth rate, so that the model
%     grows at most as fast as measured. 'min' sets them to the growth rate and 1000, so that the
%     model grows at least that fast, which suits the minimisation of a resource such as the
%     protein pool (default 'max').
% looseStrictFlux : char or double
%     how tightly the exchange fluxes are constrained (default 'loose'):
%
%     - 'loose' : one bound is 0 and the other the measured flux, which allows a flux of zero.
%     - 0-100 : the bounds are set to the measured flux times (1 -/+ p/200), a band of p percent
%       in total, so 10 gives 0.95 and 1.05 times the flux. This forces a flux through the
%       reaction.
% cSource : char
%     id of the preferred carbon source uptake reaction. Its bounds are set to 0 and 0 before the
%     flux data is applied, so that the data alone determines the carbon source. [] or '' blocks
%     no reaction, and a reaction that is not in the model is ignored (default []).
%
% Returns
% -------
% model : struct
%     model with the constrained fluxes.
%
% Notes
% -----
% The function raises applyFluxDataConstraints:invalidMaxMinGrowth,
% applyFluxDataConstraints:invalidLooseStrictFlux, applyFluxDataConstraints:unknownCondition (a
% name that is not in fluxData.conds), applyFluxDataConstraints:conditionOutOfRange and
% applyFluxDataConstraints:missingReaction (an exchange reaction of fluxData that is not in the
% model, or the biomass reaction).
%
% Examples
% --------
%     model = applyFluxDataConstraints(model, fluxData, 'growth');
%     model = applyFluxDataConstraints(model, fluxData, 'growth', ...
%         'condition', 'glucose', 'looseStrictFlux', 10, 'cSource', 'EX_glc');
%
% See also
% --------
% loadFluxData : reads the measurements.
% simulateConditions : simulates all conditions of a data set, without changing the model.

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
