function reducedModel=removeMets(model,metsToRemove,varargin)
% removeMets  Delete a set of metabolites from a model.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% metsToRemove : cell or logical or double
%     either a cell array of metabolite IDs, a logical vector with the same
%     number of elements as metabolites in the model, or a vector of
%     indexes to remove.
%
% Name-Value Arguments
% --------------------
% isNames : logical
%     true if the supplied mets are metabolite names, false if they are IDs.
%     All metabolites with a listed name are deleted, which removes
%     metabolites in several compartments at once. Requires metsToRemove to
%     be a cell array (default false).
% removeUnusedRxns : logical
%     remove reactions that involve no metabolite (default false).
% removeUnusedGenes : logical
%     remove genes that no remaining reaction uses. Takes effect only when
%     removeUnusedRxns is true (default false).
% removeUnusedComps : logical
%     remove compartments that contain no metabolite (default false).
%
% Returns
% -------
% reducedModel : struct
%     an updated model structure.
%
% Examples
% --------
%     reducedModel = removeMets(model, metsToRemove, isNames, ...
%         removeUnusedRxns, removeUnusedGenes, removeUnusedComps);
if ~islogical(metsToRemove) && ~isnumeric(metsToRemove)
    metsToRemove=convertCharArray(metsToRemove);
end

p=parseRAVENargs(varargin, {'isNames',false; 'removeUnusedRxns',false; 'removeUnusedGenes',false; 'removeUnusedComps',false});
isNames=p.isNames;
removeUnusedRxns=p.removeUnusedRxns;
removeUnusedGenes=p.removeUnusedGenes;
removeUnusedComps=p.removeUnusedComps;

%Check that metsToRemove is a cell array
if isNames==true && ~iscell(metsToRemove)
    error('Must supply a cell array of strings if isNames=true');
end

reducedModel=model;

if isNames==false
    indexesToDelete=getIndexes(model,metsToRemove,'mets');
else
    indexesToDelete=[];
    for i=1:numel(metsToRemove)
        indexesToDelete=[indexesToDelete;find(strcmp(metsToRemove(i),model.metNames))];
    end
end

%Remove metabolites
if ~isempty(indexesToDelete)
    reducedModel.mets(indexesToDelete)=[];
    reducedModel.S(indexesToDelete,:)=[];
    if isfield(reducedModel,'b')
        reducedModel.b(indexesToDelete,:)=[];
    end
    if isfield(reducedModel,'metNames')
        reducedModel.metNames(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metComps')
        reducedModel.metComps(indexesToDelete)=[];
    end
    if isfield(reducedModel,'inchis')
        reducedModel.inchis(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metSmiles')
        reducedModel.metSmiles(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metFormulas')
        reducedModel.metFormulas(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metMiriams')
        reducedModel.metMiriams(indexesToDelete)=[];
    end
    if isfield(reducedModel,'unconstrained')
        reducedModel.unconstrained(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metFrom')
        reducedModel.metFrom(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metCharges')
        reducedModel.metCharges(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metDeltaG')
        reducedModel.metDeltaG(indexesToDelete)=[];
    end
    if isfield(reducedModel,'metNotes')
        reducedModel.metNotes(indexesToDelete)=[];
    end    
end

%Remove unused reactions
if removeUnusedRxns==true
    %Get unused reactions
    [~, a]=find(reducedModel.S);
    rxnsToRemove=1:numel(reducedModel.rxns);
    rxnsToRemove(a)=[];
    reducedModel=removeReactions(reducedModel,rxnsToRemove,false,removeUnusedGenes);
end

%Remove unused compartments
if removeUnusedComps==true
    oldComps=reducedModel.comps;
    I=ismember(1:numel(oldComps),reducedModel.metComps);
    if ~all(I)
        reducedModel.comps(~I)=[];
        reducedModel.compNames(~I)=[];
        if isfield(reducedModel,'compOutside')
            reducedModel.compOutside(~I)=[];
        end
        if isfield(reducedModel,'compMiriams')
            reducedModel.compMiriams(~I)=[];
        end
        [~, J]=ismember(oldComps(reducedModel.metComps),reducedModel.comps);
        reducedModel.metComps=J;
    end
end
end
