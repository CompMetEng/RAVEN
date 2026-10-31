function rule=grRuleToString(tree)
% grRuleToString  Render a grRule syntax tree back to a string.
%
% Writes the standard RAVEN format, with " and " and " or " as operators. An
% AND node inside an OR node and an OR node inside an AND node are placed in
% parentheses. The result parses with parseGrRule to a rule with the same
% logic.
%
% Parameters
% ----------
% tree : struct
%     a tree from parseGrRule. An empty tree yields ''.
%
% Returns
% -------
% rule : char
%     the rendered grRule.
%
% Examples
% --------
%     grRuleToString(parseGrRule('(G1 and G2) or G3'))   % '(G1 and G2) or G3'
%
% See also: parseGrRule, standardizeGrRules

if isempty(tree)
    rule='';
    return
end
rule=nodeToString(tree);
end


function str=nodeToString(node)
switch node.type
    case 'gene'
        str=node.id;
    case 'and'
        % An OR child inside an AND must be bracketed to keep precedence.
        parts=cellfun(@(c) bracketIf(c,'or'),node.children,'UniformOutput',false);
        str=strjoin(parts,' and ');
    case 'or'
        % An AND child inside an OR needs no brackets for correctness, but
        % RAVEN writes complexes bracketed, matching standardizeGrRules and
        % the joinOR helper in removeLowScoreGenes.
        parts=cellfun(@(c) bracketIf(c,'and'),node.children,'UniformOutput',false);
        str=strjoin(parts,' or ');
    otherwise
        error('RAVEN:badGrRule',['Unexpected grRule node type: ' node.type]);
end
end


function str=bracketIf(node,type)
str=nodeToString(node);
if strcmp(node.type,type)
    str=['(' str ')'];
end
end
