function [newModel, removedRxns]=removeBadRxns(model,varargin)
% removeBadRxns  Remove reactions that enable production or consumption from nothing.
%
% Iteratively removes reactions that enable production or consumption of a metabolite
% without any uptake or excretion.
%
% Parameters
% ----------
% model : struct
%     a model structure. The model must not allow any uptake or excretion; the function
%     warns about open exchange reactions when the model has no unconstrained field.
%     importModel("filename", false) imports a model without uptake and excretion.
%
% Name-Value Arguments
% --------------------
% rxnRules : double
%     which reactions may be removed (default 1):
%
%     - 1 : remove only unbalanced reactions
%     - 2 : also remove reactions that could not be checked for mass
%       balance
%     - 3 : any reaction can be removed
% ignoreMets : cell or logical or double
%     a cell array of metabolite ids, a logical vector with one element per metabolite
%     of the model, or a vector of metabolite indexes to exclude from the analysis
%     (default []).
% isNames : logical
%     true if ignoreMets holds metabolite names (false: ids), which excludes the
%     metabolites of that name in all compartments. Applies only if ignoreMets is a
%     cell array (default false).
% balanceElements : cell
%     a cell array with the elements for which to balance the reactions, chosen from
%     the elements defined in parseFormulas (default {"C";"P";"S";"N";"O"}).
% refModel : struct
%     a reference model, intended as a copy of model with uptake and excretion allowed
%     and some objectives (such as production of biomass) constrained to a non-zero
%     flux. A supplied reference model triggers a warning that the feature is not
%     supported, followed by an error.
% ignoreIntBounds : logical
%     true to ignore internal bounds (including reversibility). Exchange reactions are
%     not affected. Ignoring the bounds finds unbalanced solutions that the default
%     constraints do not allow (default false).
% printReport : logical
%     true to print a report (default false).
% seed : double
%     seed for the random number generator. A seed makes the choice among several
%     candidate reactions reproducible (default: the current rng state is not
%     changed).
%
% Returns
% -------
% newModel : struct
%     the model structure without the problematic reactions.
% removedRxns : cell
%     a cell array with the ids of the removed reactions.
%
% Notes
% -----
% The function removes reactions that enable production or consumption of metabolites
% when no exchange reaction is used. Models inferred automatically from databases can
% contain unbalanced reactions that allow net production or consumption of metabolites
% without any consumption or excretion. One cause is a general compound that has a
% different meaning in different reactions (as DNA has in these two reactions):
%
%     dATP + dGTP + dCTP + dTTP <=> DNA + 4 PPi
%     0.25 dATP + 0.25 dGTP + 0.25 dCTP + 0.25 dTTP <=> DNA + PPi
%
% Such reactions are elementally unbalanced, but the model can contain further
% unbalanced reactions that are not problematic. The function removes as few reactions
% as possible so that the model cannot produce or consume anything from nothing. It
% calls findLeakMetabolite repeatedly, selects at random one of the elementally
% unbalanced reactions that carry flux in the solution, removes it, and iterates until
% no metabolite can be produced or consumed. The 'produce' direction runs before the
% 'consume' direction.
%
% Examples
% --------
%     [newModel, removedRxns] = removeBadRxns(model, rxnRules, ignoreMets, ...
%         isNames, balanceElements, refModel, ignoreIntBounds, printReport);

p=parseRAVENargs(varargin, {'rxnRules',1; 'ignoreMets',[]; 'isNames',false; 'balanceElements',{'C';'P';'S';'N';'O'}; 'refModel',[]; 'ignoreIntBounds',false; 'printReport',false; 'seed',[]});
if ~isempty(p.seed)
    rng(p.seed);
end
rxnRules=p.rxnRules;
ignoreMets=p.ignoreMets;
if ~isempty(ignoreMets) && ~islogical(ignoreMets) && ~isnumeric(ignoreMets)
    ignoreMets=convertCharArray(ignoreMets);
end
isNames=p.isNames;
balanceElements=convertCharArray(p.balanceElements);
refModel=p.refModel;
if ~isempty(refModel)
    EM='The feature to supply a reference model is currently not supported';
    warning('RAVEN:warning', '%s', EM);
end
ignoreIntBounds=p.ignoreIntBounds;
printReport=p.printReport;

%Check if the model has open exchange reactions and print a warning in that
%case
if ~isfield(model,'unconstrained')
    [~, I]=getExchangeRxns(model);
    if any(I)
        if any(model.lb(I)~=0) || any(model.ub(I)~=0)
            EM='The model contains open exchange reactions. This is not the intended use of this function. Consider importing your model using importModel(filename,false)';
            warning('RAVEN:warning', '%s', EM);
        end
    end
end

%Check that the model is feasible
sol=solveLP(model);
if isempty(sol.f)
    EM='The model is not feasible. Consider removing lower bounds (such as ATP maintenance)';
    error('RAVEN:badInput', '%s', EM);
end

%Check that the reference model is feasible
if any(refModel)
    sol=solveLP(refModel);
    if isempty(sol.f)
        EM='The reference model is not feasible';
        error('RAVEN:badInput', '%s', EM);
    end
end

%Initialize some stuff
removedRxns={};
balanceStructure=getElementalBalance(model);

%Check which elements to balance for
if ~isempty(setdiff(balanceElements,balanceStructure.elements.abbrevs))
    EM='Could not recognize all elements to balance for';
    error('RAVEN:badInput', '%s', EM);
end
bal=ismember(balanceStructure.elements.abbrevs,balanceElements);

%Main loop. First run for production, second for consumption
warned=false(2,1); %This is to prevent the same warning being printed multiple times if rxnRules==3
for i=1:2
    while 1
        %Make some metabolite using as few reactions as possible
        if i==1
            [solution, metabolite]=findLeakMetabolite(model,'produce', ...
                'ignoreMets',ignoreMets,'isNames',isNames,'minNrFluxes',false, ...
                'allowExcretion',true,'ignoreIntBounds',ignoreIntBounds);
            if ~isempty(solution)
                if printReport
                    fprintf(['Can make: ' model.metNames{metabolite(1)} '\n']);
                end
            else
                %If no solution could be found, then finish
                break;
            end
        else
            [solution, metabolite]=findLeakMetabolite(model,'consume', ...
                'ignoreMets',ignoreMets,'isNames',isNames,'minNrFluxes',false, ...
                'ignoreIntBounds',ignoreIntBounds);
            if ~isempty(solution)
                if printReport
                    fprintf(['Can consume: ' model.metNames{metabolite(1)} '\n']);
                end
            else
                %If no solution could be found, then finish
                break;
            end
        end
        
        %Find all reactions that are unbalanced and still carry flux
        I=find(abs(solution)>10^-8 & balanceStructure.balanceStatus>=0 & ~all(balanceStructure.leftComp(:,bal)==balanceStructure.rightComp(:,bal),2));
        
        %If there are unbalanced rxns then delete one of them and iterate
        if any(I)
            rxnToRemove=I(randperm(numel(I),1));
        else
            %If there are no unbalanced rxns in the solution
            if rxnRules==1
                if i==1
                    EM=['No unbalanced reactions were found in the solution, but the model can still make "' model.metNames{metabolite} '". Aborting search. Consider setting rxnRules to 2 or 3 for a more exhaustive search'];
                    warning('RAVEN:warning', '%s', EM);
                else
                    EM=['No unbalanced reactions were found in the solution, but the model can still consume "' model.metNames{metabolite} '". Aborting search. Consider setting rxnRules to 2 or 3 for a more exhaustive search'];
                    warning('RAVEN:warning', '%s', EM);
                end
                break;
            else
                %Find reactions which are not checked for mass balancing,
                %but that still carry flux
                I=find(abs(solution)>10^-8 & balanceStructure.balanceStatus<0);
                
                %If there are any such reactions, remove one of them and
                %iterate
                if any(I)
                    rxnToRemove=I(randperm(numel(I),1));
                else
                    if rxnRules==2
                        %This happens when all reactions used are balanced
                        %according to the metabolite formulas. This cannot
                        %be, and indicates that one or more of the formulas
                        %are wrong. Print a warning and delete any reaction
                        %with flux
                        if i==1
                            EM=['No unbalanced or unparsable reactions were found in the solution, but the model can still make "' model.metNames{metabolite} '". Aborting search. Consider setting rxnRules to 3 for a more exhaustive search'];
                            warning('RAVEN:warning', '%s', EM);
                        else
                            EM=['No unbalanced or unparsable reactions were found in the solution, but the model can still consume "' model.metNames{metabolite} '". Aborting search. Consider setting rxnRules to 3 for a more exhaustive search'];
                            warning('RAVEN:warning', '%s', EM);
                        end
                        break;
                    else
                        if i==1
                            if warned(1)==false
                                EM=['No unbalanced or unparsable reactions were found in the solution, but the model can still make "' model.metNames{metabolite} '". This indicates some error in the metabolite formulas. Removing random reactions in the solution'];
                                warning('RAVEN:warning', '%s', EM);
                                warned(1)=true;
                            end
                        else
                            if warned(2)==false
                                EM=['No unbalanced or unparsable reactions were found in the solution, but the model can still consume "' model.metNames{metabolite} '". This indicates some error in the metabolite formulas. Removing random reactions in the solution'];
                                warning('RAVEN:warning', '%s', EM);
                                warned(2)=true;
                            end
                        end
                        I=find(abs(solution)>10^-8);
                        rxnToRemove=I(randperm(numel(I),1));
                    end
                end
            end
        end
        removedRxns=[removedRxns;model.rxns(rxnToRemove)];
        if printReport
            fprintf(['\tRemoved: '  model.rxns{rxnToRemove} '\n']);
        end
        model=removeReactions(model,rxnToRemove);
        balanceStructure.balanceStatus(rxnToRemove)=[];
        balanceStructure.leftComp(rxnToRemove,:)=[];
        balanceStructure.rightComp(rxnToRemove,:)=[];
    end
end

newModel=model;
end
