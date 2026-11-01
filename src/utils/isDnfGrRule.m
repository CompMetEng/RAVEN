function tf=isDnfGrRule(rule)
% isDnfGrRule  Test whether a grRule is in disjunctive normal form.
%
% A rule in disjunctive normal form (DNF, an OR of AND-complexes) has no OR
% operator beneath an AND operator. It is a single gene, an AND-complex, or an
% OR of genes and AND-complexes. The standard RAVEN format requires DNF: the
% subunits that must all be present are joined by AND into a complex, and the
% reaction is catalysed by at least one such complex.
%
% Parameters
% ----------
% rule : char or string or struct
%     a grRule, or a tree from parseGrRule. An empty rule is trivially DNF.
%
% Returns
% -------
% tf : logical
%     true if the rule is in disjunctive normal form.
%
% Examples
% --------
%     isDnfGrRule('(G1 and G2) or G3')   % true
%     isDnfGrRule('(G1 or G2) and G3')   % false
%
% See also: parseGrRule, grRuleToDNF, findPotentialErrors

if ~isstruct(rule)
    rule=parseGrRule(rule);
end
if isempty(rule)
    tf=true;
    return
end
tf=isDnfNode(rule);
end


function tf=isDnfNode(node)
switch node.type
    case 'gene'
        tf=true;
    case 'and'
        % An AND-complex may not contain an OR anywhere beneath it.
        tf=~any(cellfun(@containsOr,node.children));
    case 'or'
        % Every disjunct must itself be DNF.
        tf=all(cellfun(@isDnfNode,node.children));
    otherwise
        % Unknown node type: don't flag it as a problem.
        tf=true;
end
end


function tf=containsOr(node)
switch node.type
    case 'or'
        tf=true;
    case 'and'
        tf=any(cellfun(@containsOr,node.children));
    otherwise
        tf=false;
end
end
