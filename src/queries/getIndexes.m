function indexes=getIndexes(model, objects, type, varargin)
% getIndexes  Retrieve the indexes for a list of reactions or metabolites.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% objects : cell or logical or double
%     either a cell array of IDs, a logical vector with the same number of
%     elements as objects of the given type in the model, or a vector of
%     indexes. An error is raised if an ID is not found in the model, except
%     for type 'metnames', where a name that is not found gives an empty
%     index vector.
% type : char
%     'rxns', 'mets' or 'genes', depending on what to retrieve. 'metnames'
%     queries metabolite names and returns, for each name, the indexes of all
%     metabolites with that name. 'metcomps' queries specific
%     metabolites and their compartments in the format metaboliteName[comp].
%     An error is raised if such a name matches no metabolite, or more than
%     one.
%     If a model.ec structure exists (GECKO 3), 'ecenzymes', 'ecrxns' and
%     'ecgenes' are also allowed.
%
% Name-Value Arguments
% --------------------
% returnLogical : logical
%     sets whether to return a logical array or an array with the indexes
%     (default false). It cannot be true for type 'metnames', which raises
%     the error RAVEN:badInput.
%
% Returns
% -------
% indexes : logical or double
%     a logical array if returnLogical is true, otherwise a double array. For
%     'metnames' with several names, a cell array of double vectors.
%
% Notes
% -----
% If 'ecenzymes', 'ecrxns' or 'ecgenes' are used with a GECKO 3 model, then
% the indexes are from the model.ec.enzymes, model.ec.rxns or model.ec.genes
% fields, respectively.
%
% Examples
% --------
%     indexes = getIndexes(model, objects, type, returnLogical);

p=parseRAVENargs(varargin, {'returnLogical', false});
returnLogical=p.returnLogical;

if ~islogical(objects) && ~isnumeric(objects)
    objects=convertCharArray(objects);
end
type=char(type);

indexes=[];

switch type
    case 'rxns'
        searchIn=model.rxns;
    case 'mets'
        searchIn=model.mets;
    case 'genes'
        searchIn=model.genes;
    case 'metnames'
        if returnLogical==true
            error('RAVEN:badInput', ['The combination of type ''metnames'' and ' ...
                'returnLogical true is not supported, because each name can ' ...
                'match several metabolites. Use type ''mets'' or ''metcomps'' ' ...
                'for a logical vector.']);
        end
        searchIn=model.metNames;
    case 'metcomps'
        % If provided as metaboliteName[comp], then index search is quite
        % different from general approach, and therefore coded separately.
        mets=regexprep(objects,'(^.+)\[(.+)\]$','$1');
        comps=regexprep(objects,'(^.+)\[(.+)\]$','$2');
        indexes=zeros(1,numel(objects));
        for i=1:numel(objects)
            index=find(strcmp(mets(i),model.metNames));
            index=index(strcmp(comps(i),model.comps(model.metComps(index))));
            if isscalar(index)
                indexes(i)=index;
            elseif isempty(index)
                error(['Could not find object ''' objects{i} ''' in the model']);
            else
                error(['''' objects{i} ''' is ambiguous: it matches the ' ...
                    'metabolites ' strjoin(model.mets(index),', ') ...
                    '. Use metabolite ids (type ''mets'') instead.']);
            end
        end
        indexes=indexes(:);
        if returnLogical==true
            tempIndexes=false(numel(model.mets),1);
            tempIndexes(indexes)=true;
            indexes=tempIndexes;
        end
        return %None of the remaining function needs to run if metcomps
    case 'ecrxns'
        if ~isfield(model,'ec')
            error('Type %s cannot be used if no model.ec structure is present.',type)
        end
        searchIn=model.ec.rxns;
    case 'ecenzymes'
        if ~isfield(model,'ec')
            error('Type %s cannot be used if no model.ec structure is present.',type)
        end
        searchIn=model.ec.enzymes;
    case 'ecgenes'
        if ~isfield(model,'ec')
            error('Type %s cannot be used if no model.ec structure is present.',type)
        end
        searchIn=model.ec.genes;
    otherwise
        if contains(objects,{'rxns','mets','metnames','metcomps','genes','ecgenes','ecenzymes','ecrxns'})
            error('The second and third parameter provided to run getIndexes are likely switched. Note that the second parameter is the object to find the index for, while the third parameter is the object type (''rxns'', ''mets'', ''metnames'', ''metcomps'', ''genes'', ''ecgenes'', ''ecenzymes'' or ''ecrxns'').')
        else
            error('Incorrect value of the ''type'' parameter. Allowed values are ''rxns'', ''mets'', ''metnames'', ''metcomps'', ''genes'', ''ecgenes'', ''ecenzymes'' or ''ecrxns''.');
       end
end

if iscell(objects)
    if strcmpi(type,'metnames')
        % metnames allows multiple matches per name (same metabolite in different
        % compartments), so return a cell array of index vectors.
        for i=1:numel(objects)
            indexes{i}=find(strcmp(objects(i),searchIn));
        end
    else
        % Build a hash map once (O(m)) for O(1) per-object lookup instead of
        % O(m) strcmp per object.
        searchMap = containers.Map(searchIn, 1:numel(searchIn));
        for i=1:numel(objects)
            if isKey(searchMap, objects{i})
                indexes(i) = searchMap(objects{i});
            else
                error(['Could not find object ''' objects{i} ''' in the model']);
            end
        end
    end
else
    %Now it's either a logical mask or an array with indexes. We want indexes.
    if islogical(objects)
        indexes=find(objects);
    else
        indexes=objects;
    end
end

if returnLogical==true
    tempIndexes=false(numel(searchIn),1);
    tempIndexes(indexes)=true;
    indexes=tempIndexes;
end

indexes=indexes(:);
if iscell(indexes) && length(indexes)==1
    indexes=indexes{1};
end
end
