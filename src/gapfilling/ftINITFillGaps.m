function [addedRxns, newModel, exitFlag]=ftINITFillGaps(tModel, origModel, tRefModel,allowNetProduction,supressWarnings,rxnScores,params,verbose,resolveTies)
% ftINITFillGaps
%   Variant of fillGaps adapted to speed up the generation of ftINIT models.
%
%   tModel              model that contains the task-specific rxns.
%   origModel           model without task-specific rxns.         
%   tRefModel           reference tModel: the full tModel, containing all rxns
%                       used for gap-filling of tModel + the task-specific rxns
%   allowNetProduction  true if net production of all metabolites is
%                       allowed. A reaction can be unable to carry flux because one of
%                       the reactants is unavailable or because one of the
%                       products cannot be processed further. If this
%                       parameter is true, only the first cause of blocked
%                       flux is considered (optional, default false)
%   supressWarnings     false if warnings should be displayed (optional, default
%                       false)
%   rxnScores           scores for each of the reactions in the
%                       reference tModel. The solver maximises the sum of
%                       the scores of the included reactions (optional,
%                       default -1 for all reactions)
%   params              solver parameters, passed to getMinNrFluxes (optional)
%   verbose             true to show the MILP progression.
%   resolveTies         true to fix the degenerate optimum of the gap-fill MILP
%                       to a canonical answer, which otherwise depends on the
%                       solver Seed. See resolveTies of getMinNrFluxes
%                       (optional, default false)
%
%   addedRxns           the rxns added
%   newModel            the tModel with reactions added to fill gaps
%   exitFlag            1: optimal solution found
%                      -1: no solution found, either because the problem is
%                          infeasible or because the solver reached its time
%                          limit before finding one
%                      -2: a solution was found but is not proven optimal,
%                          because the solver stopped before reaching
%                          optimality
%
%   This method works on the reference model, which already contains the
%   tModel reactions. Reactions constrained to zero flux are removed
%   (simplifyModel(tRefModel,false,false,true), deleteZeroInterval). It then
%   solves the MILP problem of including as few additional reactions as
%   possible from the reference model so that the constraints already set
%   on tModel (for example a required biomass flux, or exchange bounds) can
%   be satisfied.
%
%   fillGaps merges the reference model with the model on every call, so it
%   does not need a pre-merged reference model.
%
%   See also: fillGaps
%
% Usage: [addedRxns, newModel, exitFlag]=...
%           ftINITFillGaps(tModel,origModel,tRefModel,allowNetProduction,...
%           supressWarnings,rxnScores,params,verbose,resolveTies)

if nargin<9 || isempty(resolveTies)
    resolveTies=false;
end
if isempty(rxnScores)
    rxnScores=ones(numel(tRefModel.rxns),1)*-1;
end

%Simplify the template models to remove constrained rxns. At the same time,
%check that the id of the template models isn't the same as the tModel. That
%would cause an error further down
tRefModel.rxnScores=rxnScores;
tRefModel=simplifyModel(tRefModel,false,false,true);

%This is a rather ugly solution to the issue that it's a bit tricky to keep
%track of which scores belong to which reactions. This requires that
%removeReactions and mergeModels are modified to check for the new field.
tModel.rxnScores=zeros(numel(tModel.rxns),1);

%First merge all models into one big one
fullModel = tRefModel;%mergeModels([{tModel};models],'metNames',true);

%Add that net production is ok
if allowNetProduction==true
    %A second column in tModel.b means that the b field is lower and upper
    %bound on the RHS.
    tModel.b=[tModel.b(:,1) inf(numel(tModel.mets),1)];
    fullModel.b=[fullModel.b(:,1) inf(numel(fullModel.mets),1)];
end


fullModel.c(:)=0;

%Then minimize for the number of fluxes used. The fixed rxns doesn't need
%to participate
templateRxns = find(~ismember(fullModel.rxns, tModel.rxns)); %Check if this is slow, in that case keep track of this in fitTasksOpt and send it in

if isempty(templateRxns)
    %The reference model holds no reaction that the model does not already
    %have, so no set of additions can satisfy the constraints. This has to
    %be caught first: an empty toMinimize means "all reactions" to
    %getMinNrFluxes, and its solution would then be indexed against an
    %empty templateRxns.
    addedRxns = {};
    newModel = origModel;
    exitFlag = -1;
    return;
end

%The reversible formulation uses one binary per reaction rather than one
%per irreversible reaction, which is what makes this fast enough to run
%once per task on a genome-scale reference model
[~, J, exitFlag]=getMinNrFluxes(fullModel,templateRxns,params,fullModel.rxnScores(templateRxns),'formulation','reversible','verbose',verbose,'resolveTies',resolveTies);%only the scores from the template rxns are used, so the others doesn't matter

%Remove everything except for the added ones
addedRxns = fullModel.rxns(templateRxns(J));
if isempty(addedRxns)
    %No reactions were selected, either because the MILP found no solution
    %or because none were needed. Handing an empty set to addRxns below
    %errors out, which would reach the caller as a thrown exception rather
    %than as the exitFlag it keys on.
    addedRxns = {};
    newModel = origModel;
    return;
end
rxnsToAdd.rxns = addedRxns;
rxnsToAdd.equations = constructEquations(fullModel, addedRxns);
rxnsToAdd.ub = fullModel.ub(templateRxns(J));
rxnsToAdd.lb = fullModel.lb(templateRxns(J));
newModel = addRxns(origModel,rxnsToAdd, 3, [], true); %we ignore gene rules etc here, they are not needed - we regenerate the model in the end by subtracting rxns from the full model.

end
