function model = setGAM(model, value, biomassRxn, cofactorMetNames, varargin)
% setGAM  Set the growth-associated maintenance (GAM), and optionally the non-growth maintenance (NGAM).
%
% A cell spends ATP on maintenance in two ways. GAM is the ATP hydrolysed per gram of biomass
% formed (mmol/gDCW). It scales with growth and sits in the coefficients of the cofactors of the
% ATP hydrolysis (ATP and water are consumed; ADP, phosphate and a proton are produced) in the
% biomass pseudoreaction. NGAM is the ATP hydrolysed per hour and gram of biomass (mmol/gDCW/h) at
% any growth rate, and is the fixed flux of a maintenance reaction.
%
% For every metabolite in the biomass pseudoreaction whose model.metNames entry is in
% cofactorMetNames, the stoichiometric coefficient is set to value with the sign of the current
% coefficient. With keepOffsetsOf the coefficients keep their differences instead: the coefficient
% of the metabolite of that name (ATP) becomes value (with its sign) and every other cofactor
% changes by the same amount in the direction of its sign. A proton coefficient that balances the
% charge of the other biomass terms then keeps its offset from the ATP coefficient: for ATP, water,
% ADP and phosphate at 60 and a proton at 57, setting GAM to 45 gives 45 and a proton at 42,
% whereas without keepOffsetsOf the proton becomes 45.
%
% If ngamRxn and ngamValue are given, both bounds of the NGAM reaction are set to ngamValue. An
% error is raised before the model is changed if ngamRxn is given and is not in the model, if
% ngamValue is given without ngamRxn (RAVEN:badInput), and if keepOffsetsOf is not one of
% cofactorMetNames or does not name exactly one metabolite of the biomass pseudoreaction.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% value : double
%     new GAM, mmol ATP per gDCW of biomass formed.
% biomassRxn : char
%     id of the biomass pseudoreaction whose coefficients hold GAM.
% cofactorMetNames : cell
%     names (matched against model.metNames, not metabolite ids) of the metabolites of the GAM
%     term, for example {'ATP','ADP','H2O','H+','phosphate'}.
%
% Name-Value Arguments
% --------------------
% keepOffsetsOf : char
%     name of the metabolite whose coefficient becomes value while the other cofactors shift by
%     the same amount. Must be one of cofactorMetNames and in the biomass pseudoreaction. Empty
%     (default) sets every cofactor to value.
% ngamRxn : char
%     id of the NGAM reaction. Required when ngamValue is supplied.
% ngamValue : double
%     NGAM, mmol ATP per gDCW and hour: the flux at which the NGAM reaction is fixed.
%
% Returns
% -------
% model : struct
%     modified model.
%
% Examples
% --------
%     cofactors = {'ATP','ADP','H2O','H+','phosphate'};
%     model = setGAM(model, 45, 'r_4041', cofactors, 'keepOffsetsOf', 'ATP');
%     model = setGAM(model, 45, 'r_4041', cofactors, 'ngamRxn', 'r_4046', 'ngamValue', 0.75);
%
% See also
% --------
% getGAM : reads the value back.
% fitParameters : derives such a value from measured fluxes; setGAM then applies it.

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
