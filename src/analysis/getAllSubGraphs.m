function subGraphs=getAllSubGraphs(model)
% getAllSubGraphs  Get all metabolic subgraphs in a model.
%
% Two metabolites are connected if they share a reaction.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Returns
% -------
% subGraphs : logical
%     a logical matrix with one row per metabolite and one column per
%     subgraph. An element is true if the metabolite belongs to the
%     subgraph. The columns are ordered by decreasing subgraph size.
%
% Examples
% --------
%     subGraphs = getAllSubGraphs(model);

%Generate the connectivity graph. Metabolites are connected through
%reactions. This is not a bipartite graph with the reactions.
G=model.S*model.S';
G(G~=0)=1;

%Keeps track of which mets have been assigned to a subgraph
isAssigned=false(numel(model.mets),1);

%Allocate space for subgraphs, initially one graph for each met
subGraphs=false(numel(model.mets),numel(model.mets));

%Main loop continues until all mets have been assigned to a subgraph
counter=1;
while ~all(isAssigned)
    currentSG=false(numel(model.mets),1);
    %Find the first non-assigned metabolite and assign it to the current SG
    I=find(isAssigned==false,1);
    currentSG(I)=true;
    isAssigned(I)=true;
    %Then iteratively add all mets that are connected to each other, until
    %no more such mets can be found
    while true
        J=sum(G(:,currentSG),2);
        
        %Get the new mets for this SG. Also assign them to the current SG
        newOnes=J~=0 & currentSG==false;
        isAssigned(newOnes)=true;
        currentSG(newOnes)=true;
        
        %If there are no new mets, then abort
        if ~any(newOnes)
            subGraphs(currentSG,counter)=true;
            counter=counter+1;
            break;
        end
    end
end
subGraphs=subGraphs(:,1:counter);

[~, I]=sort(sum(subGraphs),'descend');

subGraphs=subGraphs(:,I);

%Also remove empty subgraphs (can happen when metabolites are never used)
subGraphs(:,sum(subGraphs)==0)=[];
end
