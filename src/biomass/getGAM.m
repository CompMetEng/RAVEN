function gam = getGAM(model, biomassRxn, varargin)
% getGAM  Get the growth-associated maintenance (GAM) of a model.
%
% Returns the coefficient of ATP in the biomass reaction, the value that setGAM sets.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% biomassRxn : char
%     id of the biomass reaction.
%
% Name-Value Arguments
% --------------------
% atpMetName : char
%     name of ATP in model.metNames (default 'ATP').
%
% Returns
% -------
% gam : double
%     size of the ATP coefficient (mmol/gDCW). Raises getGAM:missingBiomassRxn if the biomass
%     reaction is not in the model and getGAM:badGamTerm if it has no or several metabolites with
%     that name.
%
% Examples
% --------
%     gam = getGAM(model, 'r_4041');
%
% See also
% --------
% setGAM

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
