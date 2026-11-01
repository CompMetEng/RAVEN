function [notProduced, notProducedNames, neededForProductionMat,minToConnect,model]=checkProduction(model,varargin)
% checkProduction  Check which metabolites can be produced from a model.
%
% Checks which metabolites can be produced from a model under the current
% constraints.
%
% The function identifies which metabolites must be connected in order to
% have a fully connected network. It first identifies which metabolites
% can have a net production in the network. It then calculates which other
% metabolites must be able to have net production in order to have
% production of all metabolites in the network. If a network contains the
% equations A[external]->B, C->D, and D->E, production of C connects the
% metabolites D and E.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% checkNeededForProduction : logical
%     true to calculate, for each metabolite that could not be produced,
%     which of the other such metabolites can be produced when an uptake
%     reaction for it is added (default false).
% excretionFromCompartments : cell
%     cell array with compartment ids from which metabolites can be excreted
%     (default model.comps).
% printDetails : logical
%     true to print the metabolites that could not be produced to the command
%     window (default true).
%
% Returns
% -------
% notProduced : double
%     indices of metabolites that could not be produced, counted among the
%     metabolites in the compartments of excretionFromCompartments, in model
%     order.
% notProducedNames : cell
%     cell array with names and compartments for metabolites that could not
%     be produced.
% neededForProductionMat : logical
%     matrix where n x m is true if metabolite n allows for production of
%     metabolite m. The rows and columns follow the order of notProduced.
%     It is empty if checkNeededForProduction is false.
% minToConnect : cell
%     metabolites that need to be connected in order to be able to produce
%     all other metabolites, each with the number of metabolites it
%     connects. The metabolites are selected iteratively, each time the one
%     that connects the most remaining metabolites, so the number of
%     metabolites need not be the global minimum. It is empty if
%     checkNeededForProduction is false.
% model : struct
%     model structure with an excretion reaction added for each metabolite
%     in the compartments of excretionFromCompartments. All added reactions
%     are excretion reactions, also when checkNeededForProduction is true.
%
% Examples
% --------
%     [notProduced, notProducedNames, neededForProductionMat, minToConnect, model] = ...
%         checkProduction(model, checkNeededForProduction, excretionFromCompartments, printDetails);

p=parseRAVENargs(varargin, {'checkNeededForProduction',false; ...
    'excretionFromCompartments',[]; ...
    'printDetails',true});
checkNeededForProduction=p.checkNeededForProduction;
excretionFromCompartments=p.excretionFromCompartments;
printDetails=p.printDetails;

if isempty(excretionFromCompartments)
    excretionFromCompartments=model.comps;
else
    excretionFromCompartments=convertCharArray(excretionFromCompartments);
end

%Add an exchange reaction for each metabolite in the allowed compartments
%and see if it can carry a flux
allowedMetIds=ismember(model.comps(model.metComps),excretionFromCompartments);
allowedMetIndexes=find(allowedMetIds);
[model, addedRxns]=addExchangeRxns(model,'out',allowedMetIndexes);

canProduce=haveFlux(model,10^-5,addedRxns);

notProduced=find(~canProduce);
minToConnect={};
if checkNeededForProduction==true
    %For each of the metabolites that couldn't be produced allow uptake and
    %check which of the other metabolites that couldn't be produced that
    %can be produced
    neededForProductionMat=false(numel(notProduced));
    for i=1:numel(notProduced)
        %Add uptake for this metabolite
        if i>1
            %Reset last iteration
            model.S(:,numel(model.rxns)-numel(addedRxns)+notProduced(i-1))=model.S(:,numel(model.rxns)-numel(addedRxns)+notProduced(i-1))*-1;
        end
        %Change the production reaction to an uptake reaction
        model.S(:,numel(model.rxns)-numel(addedRxns)+notProduced(i))=model.S(:,numel(model.rxns)-numel(addedRxns)+notProduced(i))*-1;
        
        %Test which of the metabolites that couldn't be produced that can
        %be produced now
        neededForProductionMat(i,:)=haveFlux(model,10^-5,addedRxns(notProduced));
    end
    if ~isempty(notProduced)
        %Reset the last iteration, so that the returned model has excretion
        %reactions only
        lastCol=numel(model.rxns)-numel(addedRxns)+notProduced(end);
        model.S(:,lastCol)=model.S(:,lastCol)*-1;
    end
    %Calculate the smallest number of metabolites that must be connected to
    %make everything connected and return their names
    
    %The algorithm is relatively straight forward. It finds the metabolite
    %that connects the most unconnected metabolites (iteratively), adds it
    %and removes the now connected metabolites until all are connected.
    %This is not guaranteed to find the global minimum
    neededForProdTemp=neededForProductionMat;
    while 1==1
        %Get all metabolites a metabolite makes connected
        totalConnected=false(size(neededForProdTemp));
        for i=1:numel(notProduced)
            totalConnected(i,:)=neededForProdTemp(i,:);
            
            lastIter=0;
            while 1==1
                [~, a]=find(neededForProdTemp(totalConnected(i,:),:));
                totalConnected(i,a)=true;
                if numel(a)==lastIter
                    break; %No more connections were possible
                else
                    lastIter=numel(a);
                end
            end
        end
        [connections, mostConnected]=max(sum(totalConnected,2));
        
        if connections>0
            %Add the most connected metabolite to the list and remove all
            %metabolites that it's connected to
            metID=allowedMetIndexes(notProduced(mostConnected));
            entry=[model.metNames{metID},'[',model.comps{model.metComps(metID)},'] (connects ' num2str(connections) ' metabolites)'];
            minToConnect=[minToConnect;entry];
            neededForProdTemp(totalConnected(mostConnected,:),:)=false;
            neededForProdTemp(:,totalConnected(mostConnected,:))=false;
        else
            break;
        end
    end
else
    neededForProductionMat=[];
end

notProducedNames=strcat(model.metNames(allowedMetIndexes(notProduced)),'[',model.comps(model.metComps(allowedMetIndexes(notProduced))),']');

if printDetails==true
    fprintf('The following metabolites could not be produced:\n');
    [notProducedNamesTemp,perm]=sort(notProducedNames);
    
    if checkNeededForProduction==true
        neededForProdTemp=neededForProductionMat(:,perm);
        neededForProdTemp=neededForProdTemp(perm,:);
        fprintf('\tIf the production of a metabolite is dependent on some other metabolites then those are printed under the name\n\n');
    end
    for i=1:numel(notProducedNamesTemp)
        fprintf([notProducedNamesTemp{i} '\n']);
        neededForProdTemp(i,i)=false; %Not neat to do this here. Prevent printing itself
        if checkNeededForProduction==true
            enablesProduction=find(neededForProdTemp(:,i));
            if any(enablesProduction)
                for j=1:numel(enablesProduction)
                    fprintf(['\t' notProducedNamesTemp{enablesProduction(j)} '\n']);
                end
            end
        end
    end
end
end
