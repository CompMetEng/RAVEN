function orderedModel=standardizeModelFieldOrder(model)
% standardizeModelFieldOrder  Order RAVEN model structure fields.
%
% Orders the fields of a RAVEN or COBRA model structure. The RAVEN order is
% the one listed in orderRavenFields.csv, described at
% https://raven-docs.readthedocs.io/en/latest/guide/model-structure.html
% A structure with the field 'rules' is ordered as a COBRA structure, using
% COBRA_structure_fields.csv. Fields that are not listed in the file are
% placed last, in alphabetical order.
%
% The content of the fields is unchanged. To convert model fields between
% RAVEN and COBRA format, use ravenCobraWrapper().
%
% Parameters
% ----------
% model : struct
%     model structure, either RAVEN or COBRA format.
%
% Returns
% -------
% orderedModel : struct
%     model structure with ordered fields.
%
% Examples
% --------
%     orderedModel = standardizeModelFieldOrder(model);
%
% See also
% --------
% ravenCobraWrapper

ravenPath=findRAVENroot();

if ~isfield(model,'rules') % Check if model is RAVEN
    fid = fopen(fullfile(ravenPath,'src','conversion','orderRavenFields.csv'));
    fields = textscan(fid,'%s','Delimiter',',','HeaderLines',0);
    fields = fields{1};
    fclose(fid);
else % If model is COBRA
    fid = fopen(fullfile(ravenPath,'src','conversion','COBRA_structure_fields.csv')); % Taken from https://github.com/opencobra/cobratoolbox/blob/develop/src/base/io/definitions/COBRA_structure_fields.csv
    fields = textscan(fid,repmat('%s',1,15),'Delimiter','\t','HeaderLines',1);
    fields = fields{1};
    fclose(fid);
end

modelfields = fieldnames(model);
order = fields(ismember(fields(:,1),modelfields));
remainingOrder = sort(setdiff(modelfields,order));
overallOrder = [columnVector(order);columnVector(remainingOrder)];
orderedModel = orderfields(model,overallOrder);
end

function vecT = columnVector(vec)
% Converts a vector to a column vector
%
% USAGE:
%
%   vecT = columnVector(vec)
%
% INPUT:
%   vec:     a vector
%
% OUTPUT:
%   vecT:    a column vector
%
% .. Authors:
%     - Markus Herrgard (COBRA Toolbox)

[n, m] = size(vec);

if (m ~= 1 && n < m) || n == 1
    vecT = vec';
else
    vecT = vec;
end
end
