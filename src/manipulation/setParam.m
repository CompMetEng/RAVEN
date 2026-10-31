function model=setParam(model, paramType, rxnList, params, var)
% setParam  Set parameters for reactions.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% paramType : char
%     the type of parameter to set. The value is matched case-insensitively,
%     except for the reset of the other objective coefficients under "obj"
%     (below):
%
%     - "lb" : lower bound.
%     - "ub" : upper bound.
%     - "eq" : both upper and lower bound (equality constraint).
%     - "obj" : objective coefficient. The objective coefficients of all
%       other reactions are set to 0 only when every element of paramType
%       contains "obj" in lower case; with another capitalisation, or when
%       paramType also holds other types, the other coefficients keep their
%       values.
%     - "rev" : reversibility (only changes the model.rev field, does not
%       affect model.lb and model.ub). A non-zero value sets the reaction
%       reversible.
%     - "var" : variance around measured bound.
%     - "unc" : unconstrained, sets the lower and upper bound to
%       model.annotation.defaultLB and model.annotation.defaultUB if these
%       are defined, otherwise to -1000 and 1000.
% rxnList : cell or double
%     a cell array of reaction IDs or a vector with their corresponding
%     indexes. IDs that are not in the model are ignored with a warning.
% params : double
%     a vector of the corresponding values, or a single value that is
%     applied to all reactions in rxnList.
% var : double, optional
%     percentage of variance around measured value, if "var" is set as
%     paramType. Defining "var" as 5 results in lb and ub at 97.5% and
%     102.5% of the provided params value (if params value is negative,
%     then lb and ub are 102.5% and 97.5%).
%
% Returns
% -------
% model : struct
%     an updated model structure.
%
% Examples
% --------
%     model = setParam(model, paramType, rxnList, params, var);

paramType=convertCharArray(paramType);
if ~any(strcmpi(paramType,{'lb','ub','eq','obj','rev','var','unc'}))
    EM=['Incorrect parameter type: "' paramType{1} '"'];
    error('RAVEN:badInput', '%s', EM);
end
if isnumeric(rxnList) || islogical(rxnList)
    rxnList=model.rxns(rxnList);
else
    rxnList=convertCharArray(rxnList);
end
if isempty(rxnList)
    return;
end

%Allow to set several parameters to the same value
if numel(rxnList)~=numel(params) && numel(params)~=1
    EM='The number of parameter values and the number of reactions must be the same';
    error('RAVEN:badInput', '%s', EM);
end

if length(rxnList)>1
    if length(paramType)==1
        paramType(1:length(rxnList))=paramType;
    end
    if length(params)==1
        params(1:length(rxnList))=params;
    end
end

%Find the indexes for the reactions in rxnList. Do not use getIndexes
%as we do not want to throw errors if matches fail
indexes=zeros(numel(rxnList),1);

[Lia,Locb] = ismember(rxnList,model.rxns);
indexes(Lia) = Locb;
if any(~Lia)
    params(~Lia)=[];
    indexes(~Lia)=[];
    paramType(~Lia)=[];
    warning('RAVEN:warning', '%s', ravenList('Reactions not present in model, will be ignored:', rxnList(~Lia)));
end

%Change the parameters
if ~isempty(indexes)
    if contains(paramType,'obj')
        model.c=zeros(numel(model.c),1); % parameter is changed, not added
    end
    for j=1:length(paramType)
        if strcmpi(paramType{j},'eq')
            model.lb(indexes(j))=params(j);
            model.ub(indexes(j))=params(j);
        end
        if strcmpi(paramType{j},'lb')
            model.lb(indexes(j))=params(j);
        end
        if strcmpi(paramType{j},'ub')
            model.ub(indexes(j))=params(j);
        end
        if strcmpi(paramType{j},'obj')
            model.c(indexes(j))=params(j);
        end
        if strcmpi(paramType{j},'rev')
            %Non-zero values are interpreted as reversible
            model.rev(indexes(j))=params(j)~=0;
        end
        if strcmpi(paramType{j},'var')
            if params(j) < 0
                model.lb(indexes(j)) = params(j) * (1+var/200);
                model.ub(indexes(j)) = params(j) * (1-var/200);
            else
                model.lb(indexes(j)) = params(j) * (1-var/200);
                model.ub(indexes(j)) = params(j) * (1+var/200);
            end
        end
        if strcmpi(paramType{j},'unc')
            if isfield(model,'annotation') && isfield(model.annotation,'defaultLB')
                lb = model.annotation.defaultLB;
            else
                lb = -1000;
            end
            if isfield(model,'annotation') && isfield(model.annotation,'defaultUB')
                ub = model.annotation.defaultUB;
            else
                ub = 1000;
            end
            model.lb(indexes(j)) = lb;
            model.ub(indexes(j)) = ub;
        end
    end
end
if any(ismember(paramType,{'lb','ub','unc'}))
    invalidBound = model.lb(indexes) > model.ub(indexes);
    if any(invalidBound)
        error(['Invalid set of bounds for reaction(s): ', strjoin(model.rxns(indexes(invalidBound)),', '), '.'])
    end
end
