function model = loadDeltaGCSV(model, varargin)
% loadDeltaGCSV  Load metDeltaG and rxnDeltaG from CSV files.
%
% Each CSV is a two-column table: identifier, deltaG. On load, rows whose
% identifier is not in the model are skipped without a message, and every
% matched value is stored as it appears in the CSV, including yeast-GEM's
% "no measurement" placeholder (10000000.0). The function does not
% interpret CSV values, so a value that is to be treated as absent has to
% be filtered after loading. An empty string ('') for metCsv or rxnCsv
% skips that side.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% metCsv : char
%     path to the metabolite deltaG CSV, or '' to skip (default '').
% rxnCsv : char
%     path to the reaction deltaG CSV, or '' to skip (default '').
%
% Returns
% -------
% model : struct
%     model structure. The metDeltaG field is set only if metCsv is given and
%     the rxnDeltaG field only if rxnCsv is given.
%
% Examples
% --------
%     model = loadDeltaGCSV(model, 'data/met_dG.csv', 'data/rxn_dG.csv');
%     model = loadDeltaGCSV(model, 'metCsv', 'data/met_dG.csv');

p=parseRAVENargs(varargin, {'metCsv',''; 'rxnCsv',''});
metCsv=p.metCsv;
rxnCsv=p.rxnCsv;

model=loadField(model, 'metDeltaG', model.mets, metCsv);
model=loadField(model, 'rxnDeltaG', model.rxns, rxnCsv);
end

function model=loadField(model, fieldname, ids, csvPath)
if isempty(csvPath); return; end
if isfield(model, fieldname)
    disp(['Existing ' fieldname ' field will be overwritten.'])
else
    model.(fieldname)=nan(numel(ids),1);
end
G=readtable(csvPath);
values=G.(G.Properties.VariableNames{2});
[a, b]=ismember(ids, G.(G.Properties.VariableNames{1}));
model.(fieldname)(a)=values(b(a));
if any(~a)
    fprintf(['Not all identifiers are matched to %s; the latter\n' ...
             'file might have to be supplemented with deltaG values for new entries.\n'], ...
             csvPath);
end
end
