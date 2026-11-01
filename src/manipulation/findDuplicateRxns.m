function pairs = findDuplicateRxns(model, varargin)
% findDuplicateRxns  Find reactions that share identical stoichiometry.
%
% Only the stoichiometry is compared; bounds, GPRs and annotations are ignored. With
% the default ignoreDirection, A→B and B→A are duplicates. The pairs identify
% reactions that could be merged (contractModel merges duplicate reactions).
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% ignoreDirection : logical
%     true to treat A→B and B→A as duplicates (default true).
%
% Returns
% -------
% pairs : double
%     Nx2 numeric array of reaction-index pairs (i, j) where reactions i
%     and j share the same (possibly negated) stoichiometry, with i < j.
%     Empty if the model has no duplicates.
%
% Examples
% --------
%     pairs = findDuplicateRxns(model);
%     pairs = findDuplicateRxns(model, false);

p=parseRAVENargs(varargin, {'ignoreDirection',true});
ignoreDirection=p.ignoreDirection;

S=model.S';  % nRxns x nMets, kept sparse: unique(...,'rows') accepts sparse input directly

if ignoreDirection
    % Canonical form: flip rows whose first nonzero entry is negative so
    % that A->B and B->A hash to the same fingerprint.
    for i=1:size(S,1)
        k=find(S(i,:),1);
        if ~isempty(k) && S(i,k)<0
            S(i,:)=-S(i,:);
        end
    end
end

[~,~,ic]=unique(S,'rows');  % ic(i) = group label for reaction i

pairs=zeros(0,2);
for g=1:max(ic)
    idx=sort(find(ic==g));
    for ii=1:numel(idx)-1
        for jj=ii+1:numel(idx)
            pairs(end+1,:)=[idx(ii),idx(jj)]; %#ok<AGROW>
        end
    end
end
if ~isempty(pairs)
    pairs=sortrows(pairs);
end
end
