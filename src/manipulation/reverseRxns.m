function model=reverseRxns(model, rxns)
% reverseRxns  Reverse the direction of reactions in a model.
%
% Negates the stoichiometric coefficients and the objective coefficient of
% the reactions, and sets the lower bound to the negated upper bound and the
% upper bound to the negated lower bound. A flux in the new positive
% direction equals a flux in the original negative direction. Reaction IDs
% in rxns that are not in the model are ignored.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% rxns : cell
%     IDs of the reactions to reverse.
%
% Returns
% -------
% model : struct
%     an updated model structure.
%
% Examples
% --------
%     model = reverseRxns(model, rxns);

rxnInd = find(ismember(model.rxns,rxns));

%reverse the reactions
model.S(:,rxnInd) = model.S(:,rxnInd)*-1;
%swap the bounds
ub = model.ub(rxnInd);
model.ub(rxnInd) = -model.lb(rxnInd);
model.lb(rxnInd) = -ub;
%flip the objective
model.c(rxnInd) = -model.c(rxnInd);

end
