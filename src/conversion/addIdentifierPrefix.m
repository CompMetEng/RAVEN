function [model, hasChanged]=addIdentifierPrefix(model,varargin)
% addIdentifierPrefix  Add identifier prefixes required by SBML.
%
% If any reaction, metabolite, compartment, gene or model identifier does
% not start with a letter or _, which conflicts with SBML specifications,
% a prefix is added to all identifiers in the respective model field.
% The prefixes are:
%
%     "R_" for model.rxns,
%     "M_" for model.mets and model.id,
%     "C_" for model.comps;
%     "G_" for model.genes (and also represented in model.grRules)
%
% Parameters
% ----------
% model : struct
%     model whose identifiers should be modified.
%
% Name-Value Arguments
% --------------------
% fields : cell
%     cell array with the names of the model fields to check for the need
%     of a prefix, possible values: 'rxns', 'mets', 'comps', 'genes', 'id'
%     (default all of these fields).
%
% Returns
% -------
% model : struct
%     modified model.
% hasChanged : cell
%     cell array of char with one entry per field that received a prefix,
%     for example 'model.rxns (R_ prefix)'.
%
% Examples
% --------
%     [model, hasChanged] = addIdentifierPrefix(model, fields);

p=parseRAVENargs(varargin, {'fields',[]});
fields=p.fields;
if isempty(fields)
    fields = {'rxns','mets','comps','genes','id'};
end

modelFields = {'rxns','R_';
               'mets','M_';
               'comps','C_';
               'genes','G_';
               'id','M_'};

toChangeIdx = find(ismember(modelFields(:,1),fields));
hasChanged  = false(numel(modelFields(:,1)),1);
for i=1:numel(toChangeIdx)
    currName    = modelFields{toChangeIdx(i),1};
    currPrefix  = modelFields{toChangeIdx(i),2};
    if isfield(model,currName)
        currField   = model.(currName);
    else
        continue;
    end
    if ~all(startsWith(currField,regexpPattern('^[a-zA-Z_]')))
        currField = strcat(currPrefix, currField);
        hasChanged(toChangeIdx(i)) = true;

        if strcmp(currName,'genes') && isfield(model,'grRules')
            model.grRules = cellfun(@(r) prefixGenes(r,currPrefix), model.grRules, 'UniformOutput', false);
        end
        model.(currName) = currField;
    end
end

hasChanged = modelFields(hasChanged,:);
hasChanged = append('model.', hasChanged(:,1), ' (', hasChanged(:,2), ' prefix)');
end


function rule = prefixGenes(rule, prefix)
% Prefix every gene in a grRule, leaving "and", "or" and brackets alone. A
% gene is everything between spaces and brackets, so ids that contain
% "-", "." or ":" get one prefix, not one per word.
[tokens, separators] = regexp(rule, '[^\s()]+', 'match', 'split');
out = separators{1};
for i = 1:numel(tokens)
    if ~any(strcmpi(tokens{i}, {'and','or'}))
        tokens{i} = [prefix tokens{i}];
    end
    out = [out tokens{i} separators{i+1}]; %#ok<AGROW>
end
rule = out;
end
