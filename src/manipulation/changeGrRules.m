function model = changeGrRules(model,rxns,grRules,varargin)
% changeGrRules  Change the grRules of multiple reactions.
%
% Genes in grRules that are not in the model are added to the model.
%
% Parameters
% ----------
% model : struct
%     a model structure to change the gene association.
% rxns : char or cell
%     reaction IDs.
% grRules : char or cell
%     the replacement or additional gene association. ' and ' indicates subunits, ' or '
%     indicates isoenzymes, and brackets '()' separate different instances. A single
%     rule applies to all reactions in rxns.
%
% Name-Value Arguments
% --------------------
% replace : logical
%     true to replace the old gene association with the new association. False to
%     combine them as "(old) or (new)" (default true).
%
% Returns
% -------
% model : struct
%     an updated model structure.
%
% Examples
% --------
%     model = changeGrRules(model, rxns, grRules, replace);

p=parseRAVENargs(varargin, {'replace',true});
replace=p.replace;

rxns=convertCharArray(rxns);
grRules=convertCharArray(grRules);

if isscalar(grRules) && ~isscalar(rxns)
    grRules = repmat(grRules,1,numel(rxns));
end
if ~(numel(grRules)==numel(rxns))
    error('Number of rxns and grRules should be identical')
end

% Add genes to model
geneList = getGenesFromGrRules(grRules);
genesToAdd.genes=setdiff(geneList,model.genes); % Only keep the genes that are not yet part of the model.genes.
if ~isempty(genesToAdd.genes)
    model=addGenesRaven(model,genesToAdd); % Add genes
end
    
% Find reaction and gene indices
idx=getIndexes(model,rxns,'rxns');

% Change gene associations
if replace==true % Replace old gene associations
    model.grRules(idx)=grRules;
else % Add gene associations, add new gene rules after 'OR'.
    % A reaction with no existing grRule yet gets the new rule directly:
    % wrapping it as '() or (...)' regardless left an empty first clause,
    % a string neither grRuleToDNF nor isDnfGrRule can parse.
    existingRules=model.grRules(idx);
    hasExisting=~cellfun(@isempty,existingRules);
    combined=grRules;
    combined(hasExisting)=strcat('(',existingRules(hasExisting),') or (',grRules(hasExisting),')');
    model.grRules(idx)=combined;
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(model,true);
model.grRules = grRules;
model.rxnGeneMat = rxnGeneMat;
end
