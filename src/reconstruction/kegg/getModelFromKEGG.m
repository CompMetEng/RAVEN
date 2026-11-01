function [model,KOModel]=getModelFromKEGG(varargin)
% getModelFromKEGG  Load the global KEGG model.
%
% Loads the global KEGG reaction and gene model from keggModel.mat. When
% keggModel.mat is not present, the underlying artefacts (the gene-free
% reference model plus the KO, reaction and organism-gene relational
% tables, published as a toolbox-assets release) are downloaded and
% assembled, and the result is saved to keggModel.mat for later calls. The
% first build takes longer because the organism-gene table covers every
% KEGG organism, and needs a few hundred MB of disk.
%
% Name-Value Arguments
% --------------------
% keepSpontaneous : logical
%     include reactions labelled as "spontaneous" (default true).
% keepUndefinedStoich : logical
%     include reactions in the form n A <=> n+1 A (default true).
% keepIncomplete : logical
%     include reactions labelled as "incomplete", "erroneous" or "unclear"
%     (default true).
% keepGeneral : logical
%     include reactions labelled as "general reaction" (default false).
%
% Returns
% -------
% model : struct
%     a model structure with the reactions that the keep* options retain (by
%     default all except those labelled as "general reaction") and the
%     metabolites used in them.
% KOModel : struct
%     a model structure representing KEGG Orthology ids and their associated
%     genes. KO ids are stored as reactions.
%
% Notes
% -----
% The model output can be used as a template for fillGaps after the genes
% and rxnGeneMat fields are removed:
%     model = rmfield(model, 'genes') and model = rmfield(model, 'rxnGeneMat').
%
% See also
% --------
% getKEGGModelForOrganism, getPhylDist, fillGaps

ravenPath=findRAVENroot();

p=parseRAVENargs(varargin, {'keepSpontaneous',true; 'keepUndefinedStoich',true; ...
    'keepIncomplete',true; 'keepGeneral',false});
keepSpontaneous=p.keepSpontaneous;
keepUndefinedStoich=p.keepUndefinedStoich;
keepIncomplete=p.keepIncomplete;
keepGeneral=p.keepGeneral;

modelFile=fullfile(ravenPath,'src','reconstruction','kegg','keggModel.mat');
if isfile(modelFile)
    fprintf(['Importing the global KEGG model from ' strrep(modelFile,'\','/') '... ']);
    load(modelFile,'model','KOModel','isSpontaneous','isUndefinedStoich','isIncomplete','isGeneral');
    fprintf('COMPLETE\n');
else
    [model,KOModel,isSpontaneous,isUndefinedStoich,isIncomplete,isGeneral]=buildGlobalKEGGModel(ravenPath);
    fprintf(['Saving the global KEGG model to ' strrep(modelFile,'\','/') ' for future use... ']);
    save(modelFile,'model','KOModel','isSpontaneous','isUndefinedStoich','isIncomplete','isGeneral','-v7.3');
    fprintf('COMPLETE\n');
end

if keepSpontaneous==false
    model=removeReactions(model,intersect(isSpontaneous,model.rxns),true,true);
end
if keepUndefinedStoich==false
    model=removeReactions(model,intersect(isUndefinedStoich,model.rxns),true,true);
end
if keepIncomplete==false
    model=removeReactions(model,intersect(isIncomplete,model.rxns),true,true);
end
if keepGeneral==false
    model=removeReactions(model,intersect(isGeneral,model.rxns),true,true);
end
end
