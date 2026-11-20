function model = setGAM(model, value, biomassRxn, cofactorMetNames, varargin)
% setGAM  Set the growth-associated maintenance (GAM) coefficient.
%
% Set the growth-associated maintenance (GAM) coefficient in the biomass
% pseudoreaction, and optionally fix the non-growth maintenance (NGAM)
% reaction's bounds. The MATLAB counterpart of yeast-GEM's legacy
% changeGAM.
%
% For every metabolite in the biomass pseudoreaction whose model.metNames
% entry is in cofactorMetNames, the stoichiometric coefficient is set to
% value with the sign of the current coefficient. With keepOffsetsOf the
% coefficients keep their differences instead: the coefficient of the
% metabolite of that name becomes value (with its sign) and every other
% cofactor changes by the same amount in the direction of its sign. A proton
% coefficient that balances the charge of the other biomass terms then keeps
% its offset from the ATP coefficient. An error is raised before the model
% is changed if ngamRxn is given and is not in the model, if ngamValue is
% given without ngamRxn (RAVEN:badInput), and if keepOffsetsOf is not one
% of cofactorMetNames or does not name exactly one metabolite of the
% biomass pseudoreaction.
% For yeast-GEM, the cofactors are ATP, ADP, H2O, H+ and phosphate (ATP and
% H2O are substrates; ADP, H+ and phosphate are products).
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% value : double
%     New GAM value (mmol ATP / gDW per growth unit).
% biomassRxn : char
%     Reaction id of the biomass pseudoreaction.
% cofactorMetNames : cell
%     Cell array of metabolite names (matched against model.metNames, not
%     metabolite ids) to rescale, for example
%     {'ATP','ADP','H2O','H+','phosphate'}.
%
% Name-Value Arguments
% --------------------
% keepOffsetsOf : char
%     Name of the metabolite whose coefficient becomes value while the other
%     cofactors shift by the same amount. Must be one of cofactorMetNames
%     and in the biomass pseudoreaction. Empty (default) sets every
%     cofactor to value.
% ngamRxn : char
%     NGAM reaction id. Required when ngamValue is supplied.
% ngamValue : double
%     NGAM flux to fix. Both bounds of the NGAM reaction are set to
%     ngamValue.
%
% Returns
% -------
% model : struct
%     Modified model.
%
% Examples
% --------
%     model = setGAM(model, 80, 'r_4041', {'ATP','ADP','H2O','H+','phosphate'});
%     model = setGAM(model, 80, 'r_4041', {'ATP','ADP','H2O','H+','phosphate'}, ...
%         'keepOffsetsOf', 'ATP');
%
% See also
% --------
% getGAM : reads the value back.
% fitParameters : derives such a value from measured fluxes; setGAM then
%     applies it.

if nargin < 4
    error('setGAM:missingArgs', ...
        'biomassRxn and cofactorMetNames are required.');
end

p=parseRAVENargs(varargin, {'keepOffsetsOf',''; 'ngamRxn',[]; 'ngamValue',[]});
ngamRxn=p.ngamRxn;
ngamValue=p.ngamValue;
keepOffsetsOf=char(p.keepOffsetsOf);
%An empty id counts as not given
if ~isempty(ngamRxn)
    ngamRxn=cellstr(ngamRxn);
    ngamRxn=ngamRxn(~cellfun(@isempty,ngamRxn));
end
if ~isempty(ngamValue) && isempty(ngamRxn)
    error('RAVEN:badInput', ...
        'ngamValue is given without ngamRxn. Give the NGAM reaction id as ''ngamRxn''.');
end

bioPos = strcmp(model.rxns, biomassRxn);
if ~any(bioPos)
    error('setGAM:missingBiomassRxn', ...
        'Reaction %s not found in model.', biomassRxn);
end

setNgam = ~isempty(ngamRxn) && ~isempty(ngamValue);
if ~isempty(ngamRxn) && ~all(ismember(ngamRxn, model.rxns))
    error('setGAM:missingNgamRxn', ...
        'NGAM reaction %s not found in model.', strjoin(ngamRxn, ', '));
end

shift = 0;
if ~isempty(keepOffsetsOf)
    if ~any(strcmp(cofactorMetNames, keepOffsetsOf))
        error('setGAM:badKeepOffsetsOf', ...
            'keepOffsetsOf ''%s'' is not one of the cofactor names.', keepOffsetsOf);
    end
    shift = value - gamCoefficient(model, bioPos, biomassRxn, keepOffsetsOf);
end

for i = 1:length(model.mets)
    S_ix = model.S(i, bioPos);
    isCofactor = any(strcmp(cofactorMetNames, model.metNames{i}));
    if S_ix ~= 0 && isCofactor
        if isempty(keepOffsetsOf)
            model.S(i, bioPos) = sign(S_ix) * value;
        else
            model.S(i, bioPos) = S_ix + sign(S_ix) * shift;
        end
    end
end

if setNgam
    model = setParam(model, 'eq', ngamRxn, ngamValue);
end
end

function coef = gamCoefficient(model, bioPos, biomassRxn, metName)
% Size of the coefficient of the one metabolite named metName in the biomass
% pseudoreaction.
hit = find(model.S(:, bioPos) ~= 0 & strcmp(model.metNames, metName));
if numel(hit) ~= 1
    error('setGAM:badGamTerm', ...
        '%s has %d metabolites named ''%s''; the GAM term needs exactly one.', ...
        biomassRxn, numel(hit), metName);
end
coef = abs(full(model.S(hit, bioPos)));
end
