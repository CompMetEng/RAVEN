function newRxnScores = groupRxnScores(model, origRxnScores, origRxnIds, groupIds, origRxnsToZero)
% groupRxnScores
% Sums the reaction scores of all reactions that were merged into one by the
% linear merge.
%
% model          The model with linearly merged rxns.
% origRxnScores  The rxnScores from the model before the linear merge.
% origRxnIds     The rxn ids of the model before the linear merge.
% groupIds       The groupIds vector output from linearMerge.
%                There is one integer for each rxn in origRxnIds. 0 means
%                the reaction was not merged with any other reaction. A
%                non-zero integer means the reaction was merged with all
%                other rxns sharing the same integer.
% origRxnsToZero A logical vector saying which of the original rxns should not
%                be part of the problem; the score for these reactions is set to 0.
%                Since 0 means "excluded", an original score of 0 and a group
%                score that sums to 0 while at least one score of the group
%                is not 0 are set to 0.01.


newRxnScores = zeros(length(model.rxns),1);
[~,ia,ib] = intersect(model.rxns,origRxnIds);
grpIdsMerged = nan(length(model.rxns),1);
grpIdsMerged(ia) = groupIds(ib);
%check if any of the original scores are 0, in that case change them to 0.01 (unlikely)
origRxnScores(origRxnScores == 0) = 0.01;
%Then set the rxn scores for rxns to zero to 0
origRxnScores(origRxnsToZero) = 0;

%fill in original scores
newRxnScores(ia) = origRxnScores(ib);

for i = 1:length(model.rxns)
    %for reactions that are not merged with anything, just keep score as it is
    if grpIdsMerged(i) ~= 0
        %find all original rxns in the group
        sel = groupIds == grpIdsMerged(i);
        newRxnScores(i) = sum(origRxnScores(sel));
        if (newRxnScores(i) == 0 && any(origRxnScores(sel) ~= 0))
           %special unfortunate case, where the reactions happen to sum to 0 while some of them are nonzero
           %set to 0.01 in this case (proabably pretty unusual)
           newRxnScores(i) = 0.01;
        end
    end
end

end
