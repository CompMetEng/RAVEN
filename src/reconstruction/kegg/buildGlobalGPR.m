function [genes,rxnGeneMat]=buildGlobalGPR(rxns,koReaction,organismGeneKO)
% buildGlobalGPR  Join the genes of all KEGG organisms to reference reactions by KO id.
%
% Joins every KEGG organism's genes onto the reference reactions through their
% shared KO (KEGG Orthology) ids. The function builds two sparse incidence
% matrices (gene to KO and KO to reaction) that share one KO axis, and
% multiplies them. The join is not expanded row by row, because organismGeneKO
% can contain millions of rows (every gene of every KEGG organism).
%
% Parameters
% ----------
% rxns : cell array
%     the reference model's reaction ids (model.rxns), in model order.
% koReaction : table
%     the ko_reaction table (columns 'ko', 'reaction') from readKEGGTable.
% organismGeneKO : table
%     the organism_gene_ko table (columns 'organism', 'gene', 'ko') from
%     readKEGGTable.
%
% Returns
% -------
% genes : cell array
%     sorted unique 'organism:gene' identifiers.
% rxnGeneMat : sparse double
%     numel(rxns) by numel(genes); rxnGeneMat(i,j)=1 when gene j shares at
%     least one KO with reaction i.

%Map each ko_reaction row onto the reference model's own reaction order;
%rows naming a reaction absent from rxns (should not happen for a
%consistent artefact set, but be defensive) are dropped.
[inModel,rxnIdx]=ismember(koReaction.reaction,rxns);
koCol=koReaction.ko(inModel);
rxnIdxCol=rxnIdx(inModel);

%Shared KO axis for both incidence matrices below.
[kos,~,koGroupForRxnRow]=unique(koCol);
numKOs=numel(kos);
numRxns=numel(rxns);
koRxn=sparse(koGroupForRxnRow,rxnIdxCol,1,numKOs,numRxns);

%'organism:gene' identifiers, one row per (organism,gene,ko) triple. A
%gene can carry more than one KO, so genes are de-duplicated separately.
geneId=strcat(organismGeneKO.organism,':',organismGeneKO.gene);
[genes,~,geneRowIdx]=unique(geneId);
numGenes=numel(genes);

%Rows whose KO is not linked to any reaction cannot contribute a gene-
%reaction edge; drop them (organism_gene_ko is published pre-restricted
%to linked KOs, so this is normally a no-op).
[inKOs,koIdxForGeneRow]=ismember(organismGeneKO.ko,kos);
geneKO=sparse(geneRowIdx(inKOs),koIdxForGeneRow(inKOs),1,numGenes,numKOs);

rxnGeneMat=double((geneKO*koRxn)'>0);
end
