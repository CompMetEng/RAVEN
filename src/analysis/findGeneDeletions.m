function [genes, fluxes, originalGenes, details, grRatioMuts]=findGeneDeletions(model,varargin)
% findGeneDeletions  Delete genes and track the resulting fluxes.
%
% Deletes single genes or gene pairs, solves the model by FBA after each
% deletion, and returns the resulting fluxes. The output identifies gene
% deletion targets.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% testType : char
%     type of gene deletion (default "sgd"):
%
%     - "sgd" : single gene deletion
%     - "dgd" : double gene deletion
%
% Returns
% -------
% genes : double
%     a matrix with the genes that were deleted in each simulation (the
%     gene indexes in originalGenes). Each row corresponds to a column in
%     fluxes.
% fluxes : double
%     a matrix with the resulting fluxes. Double deletions that result in
%     an unsolvable problem have all zero flux. Single deletions that
%     result in an unsolvable problem are indicated in details instead.
% originalGenes : cell
%     the genes in the input model.
% details : double
%     a vector with one element per gene in originalGenes that states
%     whether the gene was deleted and why. Genes that only belong to
%     dead-end reactions, for example, are not deleted:
%
%     - 1 : deleted, and the problem was solved (a lethal deletion
%       with zero growth also gets 1)
%     - 2 : the single-deletion problem had no optimal solution
%       (infeasible or failed); fluxes and grRatioMuts are 0
%     - 3 : redundant (not assigned)
%     - 4 : involved in dead-end reaction
% grRatioMuts : double
%     growth rate ratio between mutated strain and wild type, matching the
%     originalGenes(genes) mutants. The elements do not map directly to
%     model.genes, unlike the output of COBRA getEssentialGenes. The
%     mapping to model.genes is obtained with:
%
%         grRatio=zeros(1,numel(model.genes));
%         grRatio(genes)=grRatioMuts;
%
% Examples
% --------
%     [genes, fluxes, originalGenes, details, grRatioMuts]=...
%         findGeneDeletions(model,testType);

originalModel=model;
p=parseRAVENargs(varargin, {'testType',[]});
testType=p.testType;
if isempty(testType)
    testType='sgd';
else
    testType=char(testType);
end

%Check that the test type is correct
if ~strcmpi(testType,'sgd') && ~strcmpi(testType,'dgd')
    EM='Incorrect test type';
    error('RAVEN:badInput', '%s', EM);
end

originalGenes=model.genes;
details=zeros(numel(model.genes),1);

%First simplify the model to reduce the size
model=simplifyModel(model,true,false,true,true);
model=removeReactions(model,{},true,true); %Removes unused genes
details(~ismember(originalGenes,model.genes))=4;

[~, geneMapping]=ismember(model.genes,originalGenes);
growthWT=solveLP(model);
growthWT=growthWT.f;

%Do single deletion. This is done here since the double deletion depends
%on which single deletions prove lethal (to reduce the size of the system)
if strcmpi(testType,'sgd') || strcmpi(testType,'dgd')
    fluxes=zeros(numel(model.rxns),numel(model.genes));
    grRatioMuts=zeros(1,numel(model.genes));
    for i=1:numel(model.genes)
        %Constrain all reactions involving the gene to 0
        tempModel=removeGenes(model,i,false,false,false);
        sol=solveLP(tempModel);

        %If the optimization terminated successfully
        if sol.stat==1
            fluxes(:,i)=sol.x;
            grRatioMuts(i)=sol.f/growthWT;
            details(geneMapping(i))=1;
        else
            fluxes(:,i)=0;
            grRatioMuts(i)=0;
            details(geneMapping(i))=2;
        end
    end
    genes=geneMapping;
end

%For double deletions
if strcmpi(testType,'dgd')
    %This is a little lazy but it is fine. Check which genes have already
    %been deleted in "sgd" analysis.
    [~, I]=ismember(originalGenes(details==1),model.genes);
    genesToModify=nchoosek(I,2);
    genes=geneMapping(genesToModify);
    grRatioMuts=zeros(1,numel(genes));
    fluxes=sparse(numel(model.rxns),size(genesToModify,1));
    for i=1:size(genesToModify,1)
        tempModel=removeGenes(model,genesToModify(i,:),false,false,false);
        sol=solveLP(tempModel);

        if sol.stat==1
            fluxes(:,i)=sol.x;
            grRatioMuts(i)=sol.f/growthWT;
        end
    end
end

%Map back to the old model
[~, I]=ismember(model.rxns,originalModel.rxns);
temp=fluxes;
fluxes=sparse(numel(originalModel.rxns),size(temp,2));
fluxes(I,:)=temp;
end
