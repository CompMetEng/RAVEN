function [noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, minToConnect,...
    neededForProductionMat, canProduceWithoutInput, canConsumeWithoutOutput, ...
    connectedFromTemplates, addedFromTemplates]=gapReport(model, varargin)
% gapReport  Perform a gap analysis and summarize the results.
%
% Runs a series of connectivity checks on the model - which reactions are
% blocked, how many isolated sub-networks exist, which metabolites cannot
% be produced (and the minimal set of metabolites that would need to be
% connected to fix that), and which metabolites can be produced/consumed
% without any input/output, with every boundary reaction closed. If
% templateModels is supplied, it also attempts to gap-fill using those
% models. A formatted report is printed to the command window; the same
% results are returned as outputs for further use.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% templateModels : cell
%     a cell array of template models to use for gap filling. If omitted,
%     the automated gap-filling section is skipped.
%
% Returns
% -------
% noFluxRxns : cell
%     reactions that cannot carry flux.
% noFluxRxnsRelaxed : cell
%     reactions that cannot carry flux even if the mass balance
%     constraint is relaxed so that net production of all metabolites is
%     allowed.
% subGraphs : struct
%     the metabolites in each of the isolated sub-networks.
% notProducedMets : cell
%     the metabolites that could not have net production.
% minToConnect : struct
%     the minimal number of metabolites that need to be connected in
%     order to be able to produce all other metabolites, and which
%     metabolites each of them connects.
% neededForProductionMat : double
%     matrix where n x m is true if metabolite n allows for production of
%     metabolite m.
% canProduceWithoutInput : cell
%     metabolites that could be produced with every boundary reaction
%     closed.
% canConsumeWithoutOutput : cell
%     metabolites that could be consumed with every boundary reaction
%     closed.
% connectedFromTemplates : cell
%     the reactions that could be connected using the template models.
% addedFromTemplates : struct
%     the reactions that were added from the template models and which
%     model they were added from.
%
% Notes
% -----
% The boundary reactions are those getExchangeRxns returns: exchanges,
% sinks and demands. For a model that keeps its boundary metabolites (an
% unconstrained field), those metabolites are deleted for the connectivity
% checks and kept, balanced, for the closed model.
%
% The model must be feasible under its own bounds and with its boundary
% closed. Otherwise gapReport stops with a RAVEN:infeasibleModel error that
% names the reactions with a forced flux; set the bounds of ATP maintenance
% and any other forced flux to zero first.
%
% Examples
% --------
%     [noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, ...
%         minToConnect, neededForProductionMat, connectedFromTemplates, ...
%         addedFromTemplates] = gapReport(model, templateModels);

p=parseRAVENargs(varargin, {'templateModels',[]});
templateModels=p.templateModels;

if isempty(templateModels)
    templateModels=[];
    connectedFromTemplates=[];
    addedFromTemplates=[];
end

%Escapes model text before it is spliced into an fprintf template, since
%a literal "%" in e.g. a model id or metabolite name would otherwise be
%misread as a format directive and swallow the remainder of the line.
esc=@(s) strrep(s,'%','%%');

fprintf(['Gap analysis for ' esc(model.id) ' - ' esc(model.name) '\n\n']);
if isfield(model,'unconstrained')
    closedModel=model;
    model=simplifyModel(model);
else
    closedModel=model;
    [~,boundary]=getExchangeRxns(closedModel);
    closedModel.lb(boundary)=0;
    closedModel.ub(boundary)=0;
end

requireFeasible(model,'under its own bounds');
requireFeasible(closedModel,'once its boundary reactions are closed');

model2=model;
model2.b=[model2.b inf(numel(model2.mets),1)];
I=haveFlux(model);
noFluxRxns=model.rxns(~I);
J=haveFlux(model2);
noFluxRxnsRelaxed=model2.rxns(~J);
bModel=removeReactions(model,~I,true,true);
cModel=removeReactions(model2,~J,true,true);
fprintf('***Overview\n');
fprintf([num2str(numel(model.rxns)-sum(I)) ' out of ' num2str(numel(model.rxns))...
    ' reactions cannot carry flux (' num2str(numel(model.rxns)-sum(J)) ' if net production of all metabolites is allowed)\n']);
fprintf([num2str(numel(model.mets)-numel(bModel.mets)) ' out of ' num2str(numel(model.mets))...
    ' metabolites are unreachable (' num2str(numel(model.mets)-numel(cModel.mets)) ' if net production of all metabolites is allowed)\n']);

fprintf('\n***Isolated subnetworks\n');
subGraphs=getAllSubGraphs(model);
fprintf(['A total of ' num2str(size(subGraphs,2)) ' isolated sub-networks are present in the model\n']);
for i=1:size(subGraphs,2)
    fprintf(['\t' num2str(i) '. ' num2str(sum(subGraphs(:,i))) ' metabolites\n']);
end

fprintf('\n***Metabolite connectivity\n');
[notProducedMets, ~, neededForProductionMat,minToConnect]=checkProduction(model,true,model.comps,false);
fprintf(['To enable net production of all metabolites, a total of ' num2str(numel(minToConnect)) ' metabolites must be connected\n']);
fprintf('Top 10 metabolites to connect:\n');
for i=1:min(10,numel(minToConnect))
    fprintf(['\t' num2str(i) '. ' esc(minToConnect{i}) '\n']);
end

fprintf('\n***Mass balancing\n');
produced=canExchange(closedModel,'produce');
canProduceWithoutInput=closedModel.mets(produced);
consumed=canExchange(closedModel,'consume');
canConsumeWithoutOutput=closedModel.mets(consumed);
fprintf([num2str(numel(canConsumeWithoutOutput)) ' metabolites could be consumed without any outputs\n' num2str(numel(canProduceWithoutInput)) ' metabolites could be produced without any inputs\n']);

if ~isempty(templateModels)
    fprintf('\n***Automated gap-filling\n');
    [connectedFromTemplates, ~, addedFromTemplates]=fillGaps(model,templateModels);
    t=esc(templateModels{1}.id);
    for i=2:numel(templateModels)
        t=[t ', ' esc(templateModels{i}.id)];
    end
    fprintf([num2str(numel(connectedFromTemplates)) ' unconnected reactions can be connected by including ' num2str(numel(addedFromTemplates)) ' reactions from\n' t '\n']);
end
end

function requireFeasible(model,condition)
%Errors, naming the forced-flux reactions, if the model is infeasible
model.c(:)=0;
sol=solveLP(model);
if sol.stat~=-1
    return
end
forced=model.rxns(model.lb>0 | model.ub<0);
hint='';
if ~isempty(forced)
    hint=[' Reactions with a forced flux: ' strjoin(forced,', ') '.'];
end
error('RAVEN:infeasibleModel',['The model is infeasible ' condition ...
    ', so gapReport cannot analyse it.' strrep(hint,'%','%%') ...
    ' Set the bounds of ATP maintenance and any other forced flux to zero first.']);
end
