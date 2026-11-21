function model = setGAM(model, value, biomassRxn, cofactorMetNames, varargin)
% setGAM  Set the growth-associated maintenance (GAM), and optionally the non-growth maintenance (NGAM).
%
% GAM is set in the coefficients of the ATP-hydrolysis cofactors of the biomass reaction. For
% every metabolite of the biomass reaction whose model.metNames entry is in cofactorMetNames, the
% coefficient becomes value, keeping its sign. With keepOffsetsOf, the coefficient of the
% metabolite with that name (for example ATP) becomes value and every other cofactor shifts by the
% same amount in the direction of its sign, so that offsets such as a proton coefficient that
% balances charge are kept. For ATP, water, ADP and phosphate at 60 and a proton at 57, setting
% GAM to 45 gives 45 and a proton at 42, whereas without keepOffsetsOf the proton becomes 45.
%
% If ngamRxn and ngamValue are given, both bounds of the NGAM reaction are set to ngamValue.
% Errors are raised before the model is changed if ngamRxn is not in the model, if ngamValue is
% given without ngamRxn (RAVEN:badInput), or if keepOffsetsOf is not in cofactorMetNames or does
% not name exactly one metabolite of the biomass reaction.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% value : double
%     new GAM (mmol ATP/gDCW).
% biomassRxn : char
%     id of the biomass reaction.
% cofactorMetNames : cell
%     names (matched against model.metNames, not ids) of the GAM cofactors, for example
%     {'ATP','ADP','H2O','H+','phosphate'}.
%
% Name-Value Arguments
% --------------------
% keepOffsetsOf : char
%     name of the metabolite that gets value while the other cofactors shift with it. Must be in
%     cofactorMetNames and in the biomass reaction. Empty sets every cofactor to value (default
%     []).
% ngamRxn : char
%     id of the NGAM reaction. Required with ngamValue.
% ngamValue : double
%     NGAM flux (mmol ATP/gDCW/h).
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
% getGAM, fitParameters

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
