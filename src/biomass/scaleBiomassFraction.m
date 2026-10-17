function model = scaleBiomassFraction(model, biomassConfig, componentName, newValue, varargin)
% scaleBiomassFraction  Rescale a biomass component to a target value.
%
% Rescale a biomass component to a target g/gDW value, optionally
% rescaling a second component so that the biomass total becomes
% targetTotal (1 g/gDW by default), whatever the total was before.
% The MATLAB counterpart of yeast-GEM's legacy scaleBioMass.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% biomassConfig : struct
%     Struct (see getBiomassFractions).
% componentName : char
%     Component to rescale.
% newValue : double
%     Target fraction in g/gDW.
%
% Name-Value Arguments
% --------------------
% balanceOut : char
%     Second component name to adjust so the biomass total becomes
%     targetTotal. It must be large enough to absorb the change: an error
%     is raised if its fraction would become zero or negative, or if it is
%     0 to begin with. Empty / omit to skip balancing.
% targetTotal : double
%     Biomass total in g/gDW that balanceOut brings the components to
%     (default 1).
%
% Returns
% -------
% model : struct
%     Modified model.
%
% Examples
% --------
%     model = scaleBiomassFraction(model, biomassConfig, 'protein', 0.5, 'carbohydrate');
%     model = scaleBiomassFraction(model, biomassConfig, 'protein', 0.5, ...
%         'balanceOut', 'carbohydrate', 'targetTotal', 0.95);
%
% See also
% --------
% getBiomassFractions

p=parseRAVENargs(varargin, {'balanceOut',''; 'targetTotal',1});
balanceOut=p.balanceOut;
targetTotal=p.targetTotal;

fractions = getBiomassFractions(model, biomassConfig);
if ~isfield(fractions, componentName)
    error('scaleBiomassFraction:unknownComponent', ...
        'biomassConfig has no component named %s', componentName);
end
current = fractions.(componentName);
if current == 0
    error('scaleBiomassFraction:zeroCurrent', ...
        ['Cannot scale %s to %g: current fraction is 0 ' ...
         '(pseudoreaction missing or empty).'], componentName, newValue);
end
factor = newValue / current;
model = scaleBiomassPseudoreaction(model, biomassConfig, componentName, factor);

if ~isempty(balanceOut)
    fractions = getBiomassFractions(model, biomassConfig);
    total = fractions.total;
    if ~isfield(fractions, balanceOut)
        error('scaleBiomassFraction:unknownComponent', ...
            'biomassConfig has no component named %s', balanceOut);
    end
    balanceCurrent = fractions.(balanceOut);
    if balanceCurrent == 0
        error('scaleBiomassFraction:zeroCurrent', ...
            ['Cannot balance with %s: current fraction is 0 ' ...
             '(pseudoreaction missing or empty).'], balanceOut);
    end
    balanceNew = balanceCurrent + (targetTotal - total);
    if balanceNew <= 0
        error('scaleBiomassFraction:negativeBalance', ...
            ['Cannot balance with %s: its fraction would become %g g/gDW. ' ...
             'It is %g g/gDW, but the other components total %g g/gDW.'], ...
            balanceOut, balanceNew, balanceCurrent, total - balanceCurrent);
    end
    balanceFactor = balanceNew / balanceCurrent;
    model = scaleBiomassPseudoreaction(model, biomassConfig, balanceOut, balanceFactor);
end
end
