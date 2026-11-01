function issues=findPotentialErrors(model)
% findPotentialErrors  Identify reactions with non-DNF GPR rules.
%
% Detects grRules that are not in disjunctive normal form (DNF), that is, rules in
% which an OR operator is nested inside an AND context. The RAVEN format requires DNF:
% groups of genes that must all be present together (enzyme complexes) are written as
% conjunctions, and the reaction is catalysed by at least one complex (a disjunction of
% conjunctions).
%
% Parameters
% ----------
% model : struct
%     a model structure with grRules and rxns fields.
%
% Returns
% -------
% issues : struct
%     struct array with one element per problematic reaction and the fields:
%
%     - index  : double, 1-based index into model.rxns.
%     - rxn    : char, reaction id.
%     - grRule : char, the problematic grRule string.
%     - reason : char, short explanation.
%
% Notes
% -----
% Detection walks the parse tree from parseGrRule, so gene ids that contain substrings
% such as "and" or "or" are not misread as operators.
%
% A rule is non-DNF when an AND operator has an OR anywhere beneath it.
% Bracketing alone does not make a rule non-DNF: "((G1 and G2) or G3)" and
% "(G1 or G2) or (G3 and G4)" are both in DNF. standardizeGrRules attempts an
% automatic repair.
%
% Rules that cannot be parsed at all are reported as issues too, with the
% parser's message as the reason.
%
% Examples
% --------
%     issues = findPotentialErrors(model);
%     if ~isempty(issues)
%         fprintf("%d reactions with non-DNF grRules\n", numel(issues));
%     end
%     indexes2check = vertcat(issues.index);

issues=struct('index',{},'rxn',{},'grRule',{},'reason',{});

if ~isfield(model,'grRules') || ~isfield(model,'rxns')
    return
end

for i=1:numel(model.grRules)
    rule=strtrim(model.grRules{i});
    if isempty(rule)
        continue
    end

    try
        isDnf=isDnfGrRule(rule);
    catch ME
        if ~strcmp(ME.identifier,'RAVEN:badGrRule')
            rethrow(ME)
        end
        % A rule that cannot be parsed is itself worth reporting: it is
        % certainly not usable by the isoenzyme/complex reasoning, and
        % silently skipping it would hide the very thing this function
        % exists to surface.
        issues(end+1,1)=struct( ...
            'index',  i, ...
            'rxn',    model.rxns{i}, ...
            'grRule', rule, ...
            'reason', ['Cannot be parsed: ' ME.message]); %#ok<AGROW>
        continue
    end

    if ~isDnf
        issues(end+1,1)=struct( ...
            'index',  i, ...
            'rxn',    model.rxns{i}, ...
            'grRule', rule, ...
            'reason', 'Non-DNF: OR operator nested inside AND context'); %#ok<AGROW>
    end
end
end
