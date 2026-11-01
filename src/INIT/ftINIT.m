function [model, metProduction, addedRxnsForTasks, deletedRxnsInINIT, fullMipRes] = ftINIT(prepData, tissue, celltype, hpaData, varargin)
% ftINIT  Generate a model using the ftINIT algorithm.
%
% Generates a model from proteomics, transcriptomics, metabolomics and
% metabolic tasks, in any combination except metabolomics alone. The function
% prepINITModel must be run first on the template model (such as the generic
% Human-GEM). It runs once per template and precalculates the parts that are
% independent of the omics data. Its output, prepData, is the input to this
% function.
%
% Parameters
% ----------
% prepData : struct
%     the prepData for the model, as returned by prepINITModel.
% tissue : char
%     tissue to score for. Must exist in hpaData.tissues or
%     transcrData.tissues.
% celltype : char
%     cell type to score for. Must exist in hpaData.celltypes or
%     transcrData.celltypes for this tissue. If empty, the maximum value
%     among all the cell types of the tissue is used.
% hpaData : struct
%     HPA data structure from parseHPA (required if transcrData is not
%     supplied, default []).
%
% Name-Value Arguments
% --------------------
% transcrData : struct
%     gene expression data structure (required if hpaData is not supplied,
%     default []). With fields:
%
%     - genes : cell array with the unique gene names.
%     - tissues : cell array with the tissue names. The list can contain
%       duplicates, as there can be multiple cell types per tissue.
%     - celltypes : cell array with the cell type names for each tissue.
%     - levels : GENESxTISSUES array with the expression level for each gene in
%       each tissue and cell type. NaN marks a gene that was not measured.
%     - threshold : a single value or a vector of gene expression thresholds,
%       above which genes are considered to be "expressed". If omitted, the
%       mean expression level of each gene across all tissues in transcrData
%       is used as the threshold.
% metabolomicsData : cell
%     cell array with metabolite names that the model should produce (default
%     []).
% INITSteps : cell
%     specifies the steps in the algorithm. For more info, see INITStepDesc and
%     getINITSteps (default getINITSteps([],'1+1'), the standard ftINIT).
% removeGenes : logical
%     if true, low-score genes are removed from grRules, unless they
%     are the only gene associated with a reaction, or a subunit of an enzyme
%     complex (see "removeLowScoreGenes" function for details). If false,
%     grRules are not modified; genes that were associated only with removed
%     reactions are not present in the final model (default true).
% useScoresForTasks : logical
%     true if the calculated reaction scores are used as weights when
%     fitting to tasks. Scores above -0.1 are set to -0.1 for this step
%     (default true).
% paramsFT : struct
%     solver parameters for the gap-filling step that fits the tasks, passed
%     to fitTasks as params (default []).
% verbose : logical
%     if true, the MILP progression is shown (default false).
% resolveTies : logical
%     if true, the degenerate MILP of each step is solved to a canonical
%     answer in addition to the solver Seed: the score objective is held at
%     its optimum, then the count of removable ("on") reactions is
%     minimised, then, among the sparsest solutions, their summed
%     reaction-id rank is minimised. This applies to every main-extraction
%     step and to the final gap-filling step (resolveTies of fitTasks); see
%     raven-gecko-parity#104 and resolve_ties in raven-toolbox. In each
%     main-extraction step except the last (ftINITInternalAlg), the fluxes
%     set the direction and magnitude at which a reaction is forced in the
%     next step. They are taken from the least-total-flux solution for that
%     on/off state and, among those, from the solution with the least
%     id-rank-weighted backward flux (a merged reaction is ranked by the
%     smallest original id in its group). The last step skips these flux
%     phases, as no later step reads its fluxes. A step whose primary solve
%     stops at the time limit keeps its primary solution without any of
%     these phases. The count and id-rank phases get half of the step's
%     TimeLimit each. A warning names each phase that does not prove optimal
%     or finds no solution. Measured with raven-toolbox's ftinit, whose phase
%     handling differs from ftINIT, on one cell line (Human-GEM v2.0.0-38,
%     seeds 1234 and 7, time limit 1800 s): two seeds with the default
%     escalation differ by 110 reactions and 17 predicted essential genes,
%     and with resolveTies alone by 41 reactions and 1 gene, at about 3.6
%     times the run time (about 14 min against 4 min) (default false).
% proveAbsGap : double
%     if set, each step is proven to this fixed absolute MILP gap in a single
%     solve, replacing the relative-gap escalation in stp.MILPParams. The
%     final near-zero-objective run of the escalation accepts an arbitrary
%     incumbent within the gap. A warning is raised when a step reaches the
%     time limit before the proof. The figures below were measured with
%     raven-toolbox's ftinit, whose phase handling differs from ftINIT, on one
%     cell line (Human-GEM v2.0.0-38, time limit 1800 s) with resolveTies on
%     and without gap-filling, in a sweep over four seeds (1234, 7, 1 and 2).
%     The step 1 objective of the escalation ended 4.0 to 12.0 units below the
%     proven optimum; with seeds 1234 and 7 alone it ended 7.0 and 12.0
%     below. A value of 5.0 stopped 2.0 below the optimum for one seed (seed
%     2), and values from 2.0 to 0.05 reached it for all four seeds and gave
%     the same models. Every primary solve of the steps (48 solves) ended
%     optimal; the flux phases phase2c-minflux (3 of 4 seeds) and
%     phase2d-direction (4 of 4 seeds) reached the time limit and were adopted
%     as unproven incumbents. With 1.0 and resolveTies the models of two seeds
%     differ by 1 reaction and 0 predicted essential genes, at about 26 times
%     the run time of the default escalation (about 100 min against 4 min);
%     the mean of a whole build was 5766 s at 0.05 and 5641 s at 1.0. The
%     study recommends 1.0 together with resolveTies when the variation of
%     resolveTies alone is too large and the run time is acceptable. Its
%     full-pipeline runs used a time limit of 600 s per step and its
%     determinism builds 1800 s (default [], which uses the escalation).
%
% Returns
% -------
% model : struct
%     the resulting model structure.
% metProduction : logical
%     array from the last MILP step that indicates which of the metabolites in
%     metabolomicsData could be produced (true) or not (false). Metabolites
%     with a producing reaction that is essential are not part of the MILP
%     and are not in the array. The array is determined before the
%     gap-filling for the defined tasks. To see which metabolites can be
%     produced in the final model, use canExchange.
% addedRxnsForTasks : cell
%     cell array of the reactions which were added in order to perform the
%     tasks.
% deletedRxnsInINIT : cell
%     cell array of the reactions of the template model that are not in the
%     final model, excluding exchange reactions (these are all kept). It
%     contains the reactions removed by the MILP steps that gap-filling did
%     not add back.
% fullMipRes : struct
%     the solver results from the last MILP step run.
%
% Notes
% -----
% The template model keeps its exchange metabolites,
% since checkTasks and fitTasks put specific constraints on them for each
% task. The INIT algorithm removes exchange metabolites if any are
% present. Use importModel(file,false) to import a model with exchange
% metabolites remaining.
%
% ftINIT passes a prodWeight of 5 (the reward for each produced metabolite of
% metabolomicsData) to ftINITInternalAlg, whose own default is 0.5.
%

p=parseRAVENargs(varargin, {'transcrData',[]; 'metabolomicsData',[]; 'INITSteps',[]; 'removeGenes',[]; 'useScoresForTasks',[]; 'paramsFT',[]; 'verbose',false; 'resolveTies',false; 'proveAbsGap',[]});
transcrData=p.transcrData;
metabolomicsData=p.metabolomicsData;
INITSteps=p.INITSteps;
removeGenes=p.removeGenes;
useScoresForTasks=p.useScoresForTasks;
paramsFT=p.paramsFT;
verbose=p.verbose;
resolveTies=p.resolveTies;
proveAbsGap=p.proveAbsGap;
if isempty(INITSteps)
    INITSteps = getINITSteps([],'1+1');
end
if isempty(removeGenes)
    removeGenes = true;
end
if isempty(useScoresForTasks)
    useScoresForTasks = true;
end
%Handle detected mets:
%Giving a bonus for secreting these metabolites does not work, since the
%metabolite secretion and uptake can be lost when we merge linearly dependent
%reactions. Instead, we need to figure out which reactions either produce or
%take up the mets, and give a bonus if any of them carry flux. To simplify
%things, we focus on reactions that produce the metabolite (since there must
%be one such reaction). We further reason that the direction does not matter -
%we can force one of these reactions in any direction - if it becomes a consumer, it will
%automatically force another producer on as well (otherwise we will have a net consumption).

if (~isempty(metabolomicsData))
    if length(unique(upper(metabolomicsData))) ~= length(metabolomicsData)
        error('RAVEN:badInput', '%s', 'Metabolomics contains the same metabolite multiple times');
    end
    metData = false(numel(metabolomicsData), length(prepData.minModel.rxns)); %one row per metabolite that is a boolean vector
    for i=1:numel(metabolomicsData)
        %Get the matching mets
        metSel = ismember(upper(prepData.refModel.metNames),upper(metabolomicsData{i}));
        prodRxnsSel = any(prepData.refModel.S(metSel,:) > 0,1) | ... %direct producers
                     (any(prepData.refModel.S(metSel,:) < 0,1) & prepData.refModel.rev.'); %reversible reactions that are consumers
        %convert the production rxns from refModel to minModel
        prepData.groupIds;
        [~,ia,ib] = intersect(prepData.minModel.rxns,prepData.refModel.rxns);
        grpIdsMerged = nan(length(prepData.minModel.rxns),1);
        grpIdsMerged(ia) = prepData.groupIds(ib);
        
        groupIdsPos = unique(prepData.groupIds(prodRxnsSel));%gets all group ids which includes a production rxn
        groupIdsPos = groupIdsPos(groupIdsPos ~= 0);%remove the 0 id, it means there is no group
        %the other option is that there is a direct match between the rxn id in minModel and refModel:
        posRxns = prepData.refModel.rxns(prodRxnsSel);
        directMatch = ismember(prepData.minModel.rxns, posRxns).';
        
        metData(i,:) = ismember(grpIdsMerged, groupIdsPos).' | directMatch;
    end
    metData = sparse(metData);
else
    metData = [];
end

% Get rxn scores and adapt them to the minimized model
origRxnScores = scoreModel(prepData.refModel,hpaData,transcrData,tissue,celltype);
origRxnScores(origRxnScores > -0.1 & origRxnScores <= 0) = -0.1;%we don't want reaction scores that are exactly 0 (or close), this causes problems in the milp
origRxnScores(origRxnScores < 0.1 & origRxnScores > 0) = 0.1;

rxnsTurnedOn = false(length(prepData.minModel.rxns),1);
fluxes = zeros(length(prepData.minModel.rxns),1);

rxnsToIgnoreLastStep = [1;1;1;1;1;1;1;1];

%We assume that all essential rxns are irrev - this is taken care of in
%prepINITModel. We then use an initial flux "from last run" of 0.1 for all 
%reactions. This is used for knowing what flux should be forced through an
%essential rxn.
fluxes = ones(length(prepData.minModel.rxns), 1).*0.1;

%resolveTies ranks a merged reaction by its group's smallest original id,
%which does not depend on which member names the merged reaction.
rankKeys = prepData.minModel.rxns;
[inRef, refInd] = ismember(rankKeys, prepData.refModel.rxns);
for i = find(inRef)'
    g = prepData.groupIds(refInd(i));
    if g ~= 0
        members = sort(prepData.refModel.rxns(prepData.groupIds == g));
        rankKeys{i} = members{1};
    end
end

for initStep = 1:length(INITSteps)
    disp(['ftINIT: Running step ' num2str(initStep)])
    stp = INITSteps{initStep};
    
    if any ((rxnsToIgnoreLastStep - stp.RxnsToIgnoreMask) < 0)
        error('RAVEN:badInput', '%s', 'RxnsToIgnoreMask may not cover rxns not covered in previous steps, but the other way around is fine.');
    end
    rxnsToIgnoreLastStep = stp.RxnsToIgnoreMask;
    %The fluxes of a step set the forced directions of the next one, so the
    %last step skips the flux phases of resolveTies.
    pinFluxes = initStep < numel(INITSteps);
    
    mm = prepData.minModel;
    
    if (~isempty(stp.MetsToIgnore))
        if (~isempty(stp.MetsToIgnore.simpleMets))
            %Here, we remove simple metabolites that will not really affect the milp but 
            %are very common in the S matrix. For example H2O, H+, etc.
            %It is also possible to leave compartments untouched, for example the i compartment in the mitochondria (for H+).
            metsToRem = ismember(mm.metNames,stp.MetsToIgnore.simpleMets.mets);
            compsToKeep = find(ismember(mm.comps, stp.MetsToIgnore.simpleMets.compsToKeep));
            metsToRem = metsToRem & ~ismember(mm.metComps, compsToKeep);
            mm.S(metsToRem,:) = 0;
        end    
    end

    %Set up the reaction scores and essential rxns
    rxnsToIgnore = getRxnsFromPattern(stp.RxnsToIgnoreMask, prepData);
    rxnScores = groupRxnScores(prepData.minModel, origRxnScores, prepData.refModel.rxns, prepData.groupIds, rxnsToIgnore);

    essentialRxns = prepData.essentialRxns;
    toRev = false(numel(mm.rxns),1);
    %Handle the results from previous steps ("ignore", "exclude", "essential")
    if strcmp(stp.HowToUsePrevResults, 'exclude')
        rxnScores(rxnsTurnedOn) = 0; %This is not used anymore in any step setup.
    elseif strcmp(stp.HowToUsePrevResults, 'essential')
        %Make all reversible reactions turned on in previous steps reversible
        %in the direction that they were previously carrying flux
        
        %first reverse the reactions that need to be reversed
        rev = mm.rev == 1;
        toRev = rxnsTurnedOn & rev & fluxes < 0;
        mm = reverseRxns(mm, mm.rxns(toRev));
        
        %Then make them irreversible
        mm.rev(rxnsTurnedOn) = 0;
        mm.lb(rxnsTurnedOn) = 0;

        essentialRxns = unique([prepData.essentialRxns;mm.rxns(rxnsTurnedOn)]);
    end

    
    fullMipRes = [];
    lastError = [];
    if ~isempty(proveAbsGap)
        %One solve proven to this fixed absolute gap, replacing the relative-gap
        %escalation below. The final near-zero-objective run of the escalation
        %accepts an arbitrary incumbent within the gap; on genome-scale Human-GEM
        %the step 1 objective of the escalation ended 4.0 to 12.0 units below the
        %proven optimum (raven-toolbox's prove_abs_gap, raven-gecko-parity#104).
        if ~isempty(stp.MILPParams) && isfield(stp.MILPParams{1}, 'TimeLimit')
            params = struct('TimeLimit', stp.MILPParams{1}.TimeLimit);
        else
            params = struct('TimeLimit', 5000);
        end
        params.MIPGap = 0;
        params.MIPGapAbs = proveAbsGap;
        try
            [deletedRxnsInINIT1, metProduction,fullMipRes,rxnsTurnedOn1,fluxes1] = ftINITInternalAlg(mm,rxnScores,metData,essentialRxns,5,stp.AllowMetSecr,stp.PosRevOff,params, [], fluxes, verbose, 0.1, resolveTies, rankKeys, pinFluxes);
            fluxes1(toRev) = -fluxes1(toRev);
            success = true;
            if fullMipRes.stat ~= 1
                %stat==2 (checkSolution.m): feasible but not proven within the time
                %limit -- an unproven incumbent, same caveat as the escalation path
                %below, but disclosed since proveAbsGap explicitly asks for a proof.
                EM=['ftINIT step ' num2str(initStep) ' hit the time limit before proving proveAbsGap; its kept reactions are an unproven incumbent. Raise TimeLimit, or loosen proveAbsGap.'];
                warning('RAVEN:warning', '%s', EM);
            end
        catch e
            success = false;
            lastError = e;
        end
    else
        mipGap = 1;
        first = true;
        success = false;
        for rn = 1:length(stp.MILPParams)
            params = stp.MILPParams{rn};
            if ~isfield(params, 'MIPGap')
                params.MIPGap = 0.0004;
            end

            if ~isfield(params, 'TimeLimit')
                params.TimeLimit = 5000;
            end

            if ~first
                %There is sometimes a problem with that the objective function becomes close to zero,
                %which leads to that a small percentage of that (which is the MIPGap sent in) is very small
                %and the MILP hence takes a lot of time to finish. We also therefore use an absolute MIP gap,
                %converted to a percentage using the last value of the objective function.
                params.MIPGap = min(max(params.MIPGap, stp.AbsMIPGaps{rn}/abs(lastObjVal)),1);
                params.seed = 1234;%use another seed, may work better

                if mipGap <= params.MIPGap
                    success = true;
                    break; %we're done - this will not happen the first time
                else
                    disp(['MipGap too high, trying with a different run. MipGap = ' num2str(mipGap) ' New MipGap Limit = ' num2str(params.MIPGap)])
                end
            end

            first = false;

            %now run the MILP
            try
                %The prodweight for metabolomics is set to 5. This value has not
                %been evaluated, but is assumed in the test cases - if changed, update the test case
                startVals = [];
                if ~isempty(fullMipRes)
                    startVals = fullMipRes.full;
                end
                [deletedRxnsInINIT1, metProduction,fullMipRes,rxnsTurnedOn1,fluxes1] = ftINITInternalAlg(mm,rxnScores,metData,essentialRxns,5,stp.AllowMetSecr,stp.PosRevOff,params, startVals, fluxes, verbose, 0.1, resolveTies, rankKeys, pinFluxes);
                %This is a bit tricky - since we reversed some reactions, those fluxes also need to be reversed
                fluxes1(toRev) = -fluxes1(toRev);

                mipGap = fullMipRes.mipgap;
                lastObjVal = fullMipRes.obj;
            catch e
                mipGap = Inf;
                lastObjVal = Inf; %we need to set something here, Inf leads to that this doesn't come into play
                lastError = e; %kept so a throw is not reported as a time-limit miss
            end

            success = mipGap <= params.MIPGap;
        end
    end

    if ~success
        if ~isempty(lastError)
            %The MILP threw rather than returning a poor gap; reporting this
            %as a time-limit miss sends users off to raise params.TimeLimit
            %for what is often a solver or input problem.
            error('RAVEN:badInput', '%s', ['Failed to find good enough solution within the time frame. The last MILP attempt failed with: ' lastError.message]);
        else
            error('RAVEN:badInput', '%s', ['Failed to find good enough solution within the time frame. MIPGap: ' num2str(mipGap)]);
        end
    end
    
    %save the reactions turned on and their fluxes for the next step
    rxnsTurnedOn = rxnsTurnedOn | rxnsTurnedOn1.';
    %The fluxes are a bit tricky - what if they change direction between the steps?
    %The fluxes are used to determine the direction in which reactions are forced on 
    %(to simplify the problem it is good if they are unidirectional).
    %We use the following strategy:
    %1. Use the fluxes from the most recent step.
    %2. If any flux is very low there (i.e. basically zero), use the flux from the previous steps
    %This could in theory cause problems, but seems to work well practically
    fluxesOld = fluxes;
    fluxes = fluxes1;
end


%get the essential rxns
essential = ismember(prepData.minModel.rxns,prepData.essentialRxns);
%So, we only add reactions where the linearly merged scores are zero for all linearly dependent reactions 
% (this cannot happen by chance, taken care of in the function groupRxnScores)
rxnsToIgn = rxnScores == 0; 
deletedRxnsInINITSel = ~(rxnsTurnedOn | rxnsToIgn | essential);
deletedRxnsInINIT = prepData.minModel.rxns(deletedRxnsInINITSel);

%Here we need to figure out which original reactions (before the linear merge) 
%that were removed. These are all reactions with the same group ids as the removed reactions
groupIdsRemoved = prepData.groupIds(ismember(prepData.refModel.rxns, deletedRxnsInINIT)); %can improve this slightly, use sel above
groupIdsRemoved = groupIdsRemoved(groupIdsRemoved ~= 0);%zero means that the reaction was not grouped, all with zeros are not a group!
rxnsToRem = union(prepData.refModel.rxns(ismember(prepData.groupIds,groupIdsRemoved)), deletedRxnsInINIT);%make a union here to include the ungrouped (unmerged) as well

initModel = removeReactions(prepData.refModel,rxnsToRem,false,true);

% remove metabolites separately to avoid removing those needed for tasks
unusedMets = initModel.mets(all(initModel.S == 0,2));
initModel = removeMets(initModel, setdiff(unusedMets, prepData.essentialMetsForTasks));

%if printReport == true
%    printScores(initModel,"INIT model statistics",hpaData,transcrData,tissue,celltype);
%    printScores(removeReactions(cModel,setdiff(cModel.rxns,rxnsToRem),true,true),"Reactions deleted by INIT",hpaData,transcrData,tissue,celltype);
%end

%The full model has exchange reactions in it. fitTasks with gapFillMode
%'preMerged' calls ftINITFillGaps, which automatically removes exchange
%metabolites (because it
%assumes that the reactions are constrained when appropriate). In this case the
%uptakes/outputs are retrieved from the task sheet instead. To prevent
%exchange reactions being used to fill gaps, they are deleted from the
%reference model here.
initModel.id = 'INITModel';

%If gaps in the model should be filled using a task list
if ~isempty(prepData.taskStruct)
    %Remove exchange reactions and reactions already included in the INIT
    %model. Exchange rxns that are essential are kept.

    %At this stage the model is fully connected and most of the genes with
    %good scores should have been included. The final gap-filling should
    %take the scores of the genes into account, so that "rather bad"
    %reactions are preferred to "very bad" reactions. However, reactions
    %with positive scores will be included even if they are not connected
    %in the current formulation. Therefore, such reactions will have to be
    %assigned a small negative score instead.
    exchRxns = getExchangeRxns(prepData.refModel);
    refModelNoExc = removeReactions(prepData.refModelWithBM,exchRxns,false,true);
    exchRxns = getExchangeRxns(initModel);
    initModelNoExc = removeReactions(closeModel(initModel),exchRxns,false,true);
    
    if useScoresForTasks == true
        %map the rxn scores to the model without exchange rxns
        [~,ia,ib] = intersect(refModelNoExc.rxns,prepData.refModel.rxns);
        rxnScores2nd = NaN(length(refModelNoExc.rxns),1);
        rxnScores2nd(ia) = origRxnScores(ib);
        %all(rxnScores2nd == refRxnScores);%should be the same, ok!
        [outModel,addedRxnMat] = fitTasks(initModelNoExc,refModelNoExc,[], ...
            'printOutput',true,'rxnScores',min(rxnScores2nd,-0.1), ...
            'taskStructure',prepData.taskStruct,'gapFillMode','preMerged', ...
            'params',paramsFT,'verbose',verbose,'resolveTies',resolveTies);
    else
        [outModel,addedRxnMat] = fitTasks(initModelNoExc,refModelNoExc,[], ...
            'printOutput',true,'rxnScores',[], ...
            'taskStructure',prepData.taskStruct,'gapFillMode','preMerged', ...
            'params',paramsFT,'verbose',verbose,'resolveTies',resolveTies);
    end
    %if printReport == true
    %    printScores(outModel,"Functional model statistics",hpaData,transcrData,tissue,celltype);
    %    printScores(removeReactions(outModel,intersect(outModel.rxns,initModel.rxns),true,true),"Reactions added to perform the tasks",hpaData,transcrData,tissue,celltype);
    %end
    
    addedRxnsForTasks = refModelNoExc.rxns(any(addedRxnMat,2));
else
    outModel = initModel;
    addedRxnMat = [];
    addedRxnsForTasks = {};
end

% The model can now perform all the tasks defined in the task list.
model = outModel;


% At this stage the model will contain some exchange reactions but probably
% not all (and maybe zero). This can be inconvenient, so all exchange
% reactions from the reference model are added, except for those which
% involve metabolites that are not in the model.

%Start from the original model, and just remove the reactions that are no longer there (and keep exchange rxns). The model we got out
%from the problem is not complete, it does not have GRPs etc.
%The logic below is a bit complicated. We identify the reactions that should be removed from the full model as 
%reactions that have been removed in the init model except the ones that were added back. In addition, we make 
%sure that no exchange rxns are removed - they can be removed in the init model if they were linearly merged with other
%reactions that were decided to be removed from the model. We want to keep all exchange rxns to make sure the tasks can
%be performed also without manipulating the b vector in the model (which is what is done in the gap-filling).
exchRxns = getExchangeRxns(prepData.refModel);
deletedRxnsInINIT = setdiff(prepData.refModel.rxns,union(union(initModel.rxns, addedRxnsForTasks), exchRxns));
outModel = removeReactions(prepData.refModel, deletedRxnsInINIT, true, true);

% If requested, attempt to remove negative-score genes from the model, 
% depending on their role (isozyme or complex subunit) in each grRule.
% See the "removeLowScoreGenes" function more more details, and to adjust
% any default parameters therein.
if ( removeGenes )
    [~, geneScores] = scoreModel(outModel,hpaData,transcrData,tissue,celltype);
    outModel = removeLowScoreGenes(outModel,geneScores);
end


model = outModel;

end

%This is for printing a summary of a model
function [rxnS, geneS] = printScores(model,name,hpaData,transcrData,tissue,celltype)
    [a, b] = scoreModel(model,hpaData,transcrData,tissue,celltype);
    rxnS = mean(a);
    geneS = mean(b,'omitnan');
    fprintf([name ':\n']);
    fprintf(['\t' num2str(numel(model.rxns)) ' reactions, ' num2str(numel(model.genes)) ' genes\n']);
    fprintf(['\tMean reaction score: ' num2str(rxnS) '\n']);
    fprintf(['\tMean gene score: ' num2str(geneS) '\n']);
    fprintf(['\tReactions with positive scores: ' num2str(100*sum(a>0)/numel(a)) '%%\n\n']);
end

function rxnsToIgnore = getRxnsFromPattern(rxnsToIgnorePattern, prepData)
    rxnsToIgnore = false(length(prepData.toIgnoreExch),1);
    if rxnsToIgnorePattern(1) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreExch; end;
    if rxnsToIgnorePattern(2) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreImportRxns; end;
    if rxnsToIgnorePattern(3) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreSimpleTransp; end;
    if rxnsToIgnorePattern(4) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreAdvTransp; end;
    if rxnsToIgnorePattern(5) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreSpont; end;
    if rxnsToIgnorePattern(6) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreS; end;
    if rxnsToIgnorePattern(7) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreCustomRxns; end;
    if rxnsToIgnorePattern(8) rxnsToIgnore = rxnsToIgnore | prepData.toIgnoreAllWithoutGPRs; end;
end

