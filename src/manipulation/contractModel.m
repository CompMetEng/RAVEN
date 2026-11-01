function [reducedModel, removedRxns, indexedDuplicateRxns]=contractModel(model,varargin)
% contractModel  Contract a model by merging identical reactions.
%
% simplifyModel runs contractModel for its deleteDuplicates option. If the duplicated
% reactions have "_EXP_*" suffixes (where * is a digit), the model has passed through
% expandModel, and the function removes these suffixes.
%
% The gene associations and EC numbers of merged reactions are combined. The merged
% reaction keeps the widest bounds and the largest objective coefficient, and a warning
% is given when the duplicates differ in these. findDuplicateRxns lists duplicate
% reactions without changing the model.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% distReverse : logical
%     true to treat reactions with the same metabolites but different reversibility as
%     different reactions (default true).
% mets : char or cell
%     metabolite ids as a string or a cell array of strings; only the reactions that
%     involve these metabolites are checked for duplication (default: all reactions).
%     replaceMets uses this option.
%
% Returns
% -------
% reducedModel : struct
%     a model structure without duplicate reactions.
% removedRxns : cell
%     cell array with the ids of the removed duplicate reactions.
% indexedDuplicateRxns : cell
%     cell array with one entry per reaction of reducedModel, holding the ids of the
%     duplicates removed into that reaction (multiple ids separated by semicolon).
%
% Examples
% --------
%     [reducedModel, removedRxns, indexedDuplicateRxns] = contractModel(model, distReverse, mets);
%
% Notes
% -----
% The function may not work for grRules strings with nested expressions of "and" and "or".

p=parseRAVENargs(varargin, {'distReverse',true; 'mets',[]});
distReverse=p.distReverse;
mets=p.mets;

%First sort the model so that reversible reactions are in the same
%direction
modelS=sortModel(model,true);

%Get a list of duplicate reactions
if distReverse
    x=[modelS.S; model.rev']';
else
    x=modelS.S';
end
if ~isempty(mets)
    toCheck = getIndexes(model,mets,'mets');
    [~, toCheck] = find(modelS.S(toCheck,:));
    x = [x, (1:size(x,1))']; % Make all other rxns unique by additional column 
    x(toCheck,end) = 0;
end

[~,I,J] = unique(x,'rows','first');

%Initialize cell array of indexedDuplicateRxns
indexedDuplicateRxns=cell(numel(model.rxns),1);
indexedDuplicateRxns(:)={''};

duplicateRxns=setdiff(1:numel(model.rxns),I);
mergeTo=I(J(duplicateRxns));

mergedRxns=unique(mergeTo);

%Now add all the info from this one. Print a warning if they have different
%bounds or objective function coefficients. Uses the widest bounds and
%largest magnitude of objective coefficient
for i=1:numel(mergedRxns)
    duplRxn=transpose([mergedRxns(i),duplicateRxns(mergeTo==mergedRxns(i))]);
    if numel(unique(model.lb(duplRxn)))>1
        EM=['Duplicates of reaction ' model.rxns{mergedRxns(i)} ' have different lower bound. Uses the most negative/smallest lower bound'];
        warning('RAVEN:warning', '%s', EM);
        model.lb(mergedRxns(i))=min(model.lb(duplRxn));
    end
    if numel(unique(model.ub(duplRxn)))>1
        EM=['Duplicates of reaction ' model.rxns{mergedRxns(i)} ' have different upper bound. Uses the most positive/largest upper bound'];
        warning('RAVEN:warning', '%s', EM);
        model.ub(mergedRxns(i))=max(model.ub(duplRxn));
    end
    if numel(unique(model.c(duplRxn)))>1
        EM=['Duplicates of reaction ' model.rxns{mergedRxns(i)} ' has a different objective function coefficient. Uses the largest coefficient'];
        warning('RAVEN:warning', '%s', EM);
        model.c(mergedRxns(i))=max(model.c(duplRxn));
    end
    if isfield(model,'grRules') && any(~isempty(model.grRules(duplRxn)))
        rules=model.grRules(duplRxn);
        allRules={};
        for j=1:numel(rules)
            rules{j}=ignoreORinComplex(rules{j});
            allRules=[allRules regexp(rules{j},' or ','split')];
        end
        allRules=unique(allRules);
        allRules=strrep(allRules,'__OR__',' or ');
        andRules=contains(allRules,' and ');
        allRules(andRules)=strcat('(',allRules(andRules),')');
        if numel(allRules)==1
            model.grRules{mergedRxns(i)}=allRules{1};
        else
            model.grRules{mergedRxns(i)}=strjoin(allRules,' or ');
        end
    end    
    if isfield(model,'eccodes') && any(~isempty(model.eccodes(duplRxn)))
        codes=model.eccodes(duplRxn);
        allCodes={};
        for j=1:numel(codes)
            allCodes=[allCodes regexp(codes{j},';','split')];
        end
        allCodes=unique(allCodes);
        if numel(allCodes)==1
            model.eccodes{mergedRxns(i)}=allCodes{1};
        else
            model.eccodes{mergedRxns(i)}=strjoin(allCodes,';');
        end
    end
    %Generate indexedDuplicateRxns cell array
    if numel(duplRxn)==2
        indexedDuplicateRxns{duplRxn(1)}=model.rxns{duplRxn(2)};
    else
        indexedDuplicateRxns{duplRxn(1)}=strjoin(model.rxns(duplRxn(2:end)),';');
    end
    %If all reactions have _EXP_* (where * is digit) as suffix, then the
    %model had passed through expandModel(), and these suffixes are removed.
    if all(~cellfun(@isempty,regexp(model.rxns(duplRxn),'_EXP_\d+$')))
        model.rxns(duplRxn(1))=regexprep(model.rxns(duplRxn(1)),'_EXP_\d+$','');
    end
end

%Delete the duplicate reactions
reducedModel=removeReactions(model,duplicateRxns);
removedRxns=model.rxns(duplicateRxns);
[~, index]=ismember(reducedModel.rxns,model.rxns);
indexedDuplicateRxns=indexedDuplicateRxns(index);

if isfield(reducedModel,'rxnGeneMat')
    %Fix grRules and reconstruct rxnGeneMat
    [grRules,rxnGeneMat] = standardizeGrRules(reducedModel,true);
    reducedModel.grRules = grRules;
    reducedModel.rxnGeneMat = rxnGeneMat;
end
end

function grRule = ignoreORinComplex(grRule)
%In a grRule, if OR relationship is nested in an AND relationship, then
%obfuscate the OR before splitting isoenzymes
grRule=['(' grRule ')'];
brOpen=strfind(grRule,'(');
brClose=strfind(grRule,')');
andPos=strfind(grRule,' and ');
%Find opening bracket closest before AND
stillCapturing = 0;
for i=1:numel(andPos)
    searchPos = andPos(i);
    while stillCapturing == 0
        closestOpen = brOpen(max(find(brOpen<searchPos)));
        inbetweenClose = brClose(brClose<searchPos & brClose>closestOpen);
        if ~isempty(inbetweenClose)
            searchPos=max(inbetweenClose);
        else
            stillCapturing = 1;
            beginPos = closestOpen;
        end
    end
    stillCapturing = 0;
    searchPos = andPos(i);
    while stillCapturing == 0
        closestClose = brClose(min(find(brClose>searchPos)));
        inbetweenOpen = brOpen(brOpen>searchPos & brOpen<closestOpen);
        if ~isempty(inbetweenOpen)
            searchPos=min(closestClose);
        else
            stillCapturing = 1;
            endPos = closestClose;
        end
    end
    replacePart=regexprep(grRule(beginPos:endPos),' or ','__OR__');
    grRule=strrep(grRule,grRule(beginPos:endPos),replacePart);
end
grRule=grRule(2:end-1);
end
