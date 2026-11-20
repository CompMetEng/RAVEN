function gam = getGAM(model, biomassRxn, varargin)
% getGAM  Get the growth-associated maintenance (GAM) of a model.
%
% GAM is the size of the ATP coefficient of the biomass pseudoreaction: the
% mmol ATP hydrolysed per gDW of biomass formed. The value that setGAM sets.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% biomassRxn : char
%     Reaction id of the biomass pseudoreaction.
%
% Name-Value Arguments
% --------------------
% atpMetName : char
%     Name of the ATP metabolite, matched against model.metNames. Default
%     'ATP'.
%
% Returns
% -------
% gam : double
%     The size of the ATP coefficient (mmol/gDW). An error is raised if the
%     biomass pseudoreaction is not in the model (getGAM:missingBiomassRxn)
%     or does not have exactly one metabolite with that name
%     (getGAM:badGamTerm).
%
% Examples
% --------
%     gam = getGAM(model, 'r_4041');
%
% See also
% --------
% setGAM : sets the value.

p = parseRAVENargs(varargin, {'atpMetName','ATP'});
bioPos = strcmp(model.rxns, biomassRxn);
if ~any(bioPos)
    error('getGAM:missingBiomassRxn', ...
        'Reaction %s not found in model.', biomassRxn);
end
hit = find(model.S(:, bioPos) ~= 0 & strcmp(model.metNames, p.atpMetName));
if numel(hit) ~= 1
    error('getGAM:badGamTerm', ...
        '%s has %d metabolites named ''%s''; the GAM term needs exactly one.', ...
        biomassRxn, numel(hit), p.atpMetName);
end
gam = abs(full(model.S(hit, bioPos)));
end
