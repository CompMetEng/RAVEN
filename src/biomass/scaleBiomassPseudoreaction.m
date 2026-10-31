function model = scaleBiomassPseudoreaction(model, biomassConfig, componentName, factor)
% scaleBiomassPseudoreaction  Rescale a biomass component pseudoreaction.
%
% Multiply the substrate coefficients of one biomass component
% pseudoreaction by factor and rebalance H+ to preserve charge neutrality.
% The MATLAB counterpart of yeast-GEM's legacy rescalePseudoReaction.
%
% "Substrate" means every metabolite in the pseudoreaction whose
% metabolite name differs from the component name; the coefficient of the
% component's product is unchanged. After rescaling, the coefficient of
% biomassConfig.proton_met is recomputed so that the total charge of the
% pseudoreaction is zero. An error is raised if the model has no
% metCharges field or if a metabolite other than proton_met in the
% pseudoreaction has an undefined (NaN) charge.
%
% Parameters
% ----------
% model : struct
%     RAVEN model struct.
% biomassConfig : struct
%     Struct (see getBiomassFractions).
% componentName : char
%     Name of the component to rescale. It must match
%     biomassConfig.components{i}.name for some i and the model.metNames
%     entry of the produced metabolite in the matching pseudoreaction.
% factor : double
%     Factor by which the substrate coefficients are multiplied.
%
% Returns
% -------
% model : struct
%     Modified model.
%
% Examples
% --------
%     model = scaleBiomassPseudoreaction(model, biomassConfig, 'protein', 0.9);
%
% See also
% --------
% getBiomassFractions

comp = findComponent(biomassConfig, componentName);
rxnPos = find(strcmp(model.rxnNames, comp.pseudoreaction_name));
if isempty(rxnPos)
    error('scaleBiomassPseudoreaction:missingPseudoreaction', ...
        'No reaction named %s in model.', comp.pseudoreaction_name);
end

Hc = find(strcmp(model.mets, biomassConfig.proton_met));
if isempty(Hc)
    error('scaleBiomassPseudoreaction:missingProton', ...
        'Proton metabolite %s not found in model.', biomassConfig.proton_met);
end

for i = 1:length(model.mets)
    S_ir = model.S(i, rxnPos);
    isProd = strcmp(model.metNames{i}, componentName);
    if S_ir ~= 0 && ~isProd
        model.S(i, rxnPos) = factor * S_ir;
    end
end

% Rebalance H+ to keep charge neutrality.
model.S(Hc, rxnPos) = 0;
model.S(Hc, rxnPos) = -chargeResidual(model, rxnPos);
end

function residual = chargeResidual(model, rxnPos)
% Sum of charges over the metabolites that take part in a reaction. Only the
% participants may be summed: S is sparse and 0*NaN is NaN, so any unset
% charge in the model would otherwise poison the result. Summing with
% 'omitnan' instead would treat an unset charge as neutral and silently
% rebalance the reaction to the wrong coefficient.
metIdx = find(model.S(:, rxnPos));
if ~isfield(model, 'metCharges')
    error('scaleBiomassPseudoreaction:noCharges', ...
        ['Cannot rebalance charge: the model has no metCharges field. ' ...
         'Provide metabolite charges, or rescale without charge balancing.']);
end
unknown = metIdx(isnan(model.metCharges(metIdx)));
if ~isempty(unknown)
    error('scaleBiomassPseudoreaction:unknownCharge', ...
        ['Cannot rebalance charge of %s: the following metabolites have no ' ...
         'charge, so the charge balance is unknown rather than zero:\n\t%s'], ...
        model.rxns{rxnPos}, strjoin(model.mets(unknown), ', '));
end
residual = sum(full(model.S(metIdx, rxnPos)) .* model.metCharges(metIdx));
end

function comp = findComponent(cfg, name)
for i = 1:numel(cfg.components)
    if strcmp(cfg.components{i}.name, name)
        comp = cfg.components{i};
        return;
    end
end
error('scaleBiomassPseudoreaction:unknownComponent', ...
    'biomassConfig has no component named %s', name);
end
