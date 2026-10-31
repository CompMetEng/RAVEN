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
% value with the sign of the current coefficient. An error is raised
% before the model is changed if ngamRxn is given and is not in the model,
% and if ngamValue is given without ngamRxn (RAVEN:badInput).
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
%
% See also
% --------
% fitParameters : derives such a value from measured fluxes; setGAM then
%     applies it.

if nargin < 4
    error('setGAM:missingArgs', ...
        'biomassRxn and cofactorMetNames are required.');
end

p=parseRAVENargs(varargin, {'ngamRxn',[]; 'ngamValue',[]});
ngamRxn=p.ngamRxn;
ngamValue=p.ngamValue;
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

for i = 1:length(model.mets)
    S_ix = model.S(i, bioPos);
    isCofactor = any(strcmp(cofactorMetNames, model.metNames{i}));
    if S_ix ~= 0 && isCofactor
        model.S(i, bioPos) = sign(S_ix) * value;
    end
end

if setNgam
    model = setParam(model, 'eq', ngamRxn, ngamValue);
end
end
