function newIds=generateNewIds(model,type,prefix,varargin)
% generateNewIds  Generate a list of new metabolite or reaction ids.
%
% The ids are sequentially numbered with a defined prefix. The model is
% queried for the highest existing number of that type of id. A suffix
% counts as a number when str2double reads it as a finite integer value;
% text, Inf, NaN, fractions (r_2.5) and numbers that overflow (r_1e400) do
% not count.
%
% Parameters
% ----------
% model : struct
%     model structure.
% type : char
%     type of id, 'rxns' or 'mets'.
% prefix : char
%     prefix to be used in all ids, e.g. 's_' or 'r_'.
%
% Name-Value Arguments
% --------------------
% quantity : double
%     number of new ids that should be generated (default 1).
% numLength : double
%     length of the numerical part of the id. E.g. 4 gives ids like r_0001
%     and 6 gives ids like r_000001. If the prefix is already used in the
%     model, then the model-defined length will be used instead
%     (default 4).
%
% Returns
% -------
% newIds : cell
%     cell array with the generated ids.
%
% Examples
% --------
%     newIds = generateNewIds(model, type, prefix, quantity, numLength);
type=char(type);
prefix=char(prefix);

if type=='rxns'
    existingIds=model.rxns;
elseif type=='mets'
    existingIds=model.mets;
else
    error('type should be either ''rxns'' or ''mets''.')
end
p=parseRAVENargs(varargin, {'quantity',1; 'numLength',4});
quantity=p.quantity;
numLength=p.numLength;

% Subset only existingIds that have the prefix. Escaped: prefix is
% arbitrary caller-supplied text, not a regex pattern.
escapedPrefix=regexptranslate('escape',prefix);
existingIds=existingIds(~cellfun(@isempty,regexp(existingIds,['^' escapedPrefix])));

if ~isempty(existingIds)
    existingIds=regexprep(existingIds,['^' escapedPrefix],'');
    numbers=str2double(existingIds);
    notInteger=~isfinite(numbers) | imag(numbers)~=0 | real(numbers)~=floor(real(numbers));
    numbers=real(numbers);
    numbers(notInteger)=NaN;
    if all(isnan(numbers))
        [~,I]=sort(existingIds);
        lastId=existingIds{I(end)};
        numLength=length(lastId);
        lastId=NaN;
    else
        [lastId,I]=max(numbers);
        numLength=length(existingIds{I});
    end
else
    lastId=0;
    fprintf(['No ' type ' ids with prefix "' prefix ...
        '" currently exist in the model. The first new id will be "' ...
        [prefix,num2str(1,['%0' num2str(numLength) 'd'])] '"\n'],'%s')
end

if isnan(lastId)
    lastId = 0;
end

newIds=cell(quantity,1);

for k=1:quantity
    newIds{k}=[prefix,num2str(k+lastId,['%0' num2str(numLength) 'd'])];
end
end