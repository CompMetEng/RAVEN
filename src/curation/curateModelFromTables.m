function newModel=curateModelFromTables(model,metsInfo,varargin)
% curateModelFromTables  Curate or add mets, rxns and genes from tables.
%
% Curate existing and/or add new metabolites, reactions and genes from
% tabular data files.
%
% If the *.tsv files contain metabolites, reactions and/or genes that are
% already present in the model, then information in the model will be
% overwritten. Note that this includes empty annotations in the *.tsv
% files! Metabolites are matched by metaboliteName[comp]; reactions by the
% stoichiometry of its reactants and products; genes by their gene name.
% This function can therefore be used to add new entities in the model, or
% curate those already existing in the model.
%
% Columns are identified by their header (case-insensitive), so their order
% does not matter and optional columns can be left out. Every field that a
% table supplies is written to new and existing entities alike; a field the
% model does not have yet is created, with its default value for the other
% entities, unless all supplied values for it are empty. An empty grRules
% cell leaves the rule of an existing reaction unchanged.
%
% Parameters
% ----------
% model : struct
%     RAVEN model structure to be curated.
% metsInfo : char
%     Path to a *.tsv file with metabolite information, or 'none' to skip
%     metabolite curation. Columns: metNames, comps (both required),
%     formula, charge, inchi, metNotes, and any number of MIRIAM-namespace
%     columns. Singular forms and the RAVEN field names (metFormulas,
%     metCharges, inchis) are accepted as headers too.
%
% Name-Value Arguments
% --------------------
% genesInfo : char
%     Path to a *.tsv file with gene information, or 'none'. Columns:
%     genes (required), geneShortNames (or geneShortName), and MIRIAM.
% rxnsCoeffs : char
%     Path to a *.tsv file with reaction stoichiometric coefficients, or
%     'none'. Columns (all required): rxnIdx (or index), rxnNames,
%     metNames, comps, coefficient. One row per (reaction, metabolite)
%     pair.
% rxnsInfo : char
%     Path to a *.tsv file with reaction information, or 'none'. Columns:
%     rxnIdx (or index), rxnNames (both required), grRules, lb, ub, rev,
%     subSystems, eccodes, rxnNotes, rxnReferences, rxnConfidenceScores,
%     and MIRIAM.
% metPrefix : char
%     Prefix used to mint fresh metabolite ids (e.g. 's_' for yeast-GEM,
%     'M_' for the BiGG default) (default 'M_').
% rxnPrefix : char
%     Prefix used to mint fresh reaction ids (default 'R_').
%
% Returns
% -------
% newModel : struct
%     Curated RAVEN model structure.
%
% Examples
% --------
%     newModel = curateModelFromTables(model, metsInfo, genesInfo, ...
%                     rxnsCoeffs, rxnsInfo, metPrefix, rxnPrefix);
%
% Notes
% -----
% In all three info tables, any column whose header is not one of the
% listed core fields is treated as a MIRIAM annotation namespace and
% stored on the matching entity.

p=parseRAVENargs(varargin, {'genesInfo','none'; 'rxnsCoeffs','none'; 'rxnsInfo','none'; 'metPrefix',[]; 'rxnPrefix',[]});
genesInfo=p.genesInfo;
rxnsCoeffs=p.rxnsCoeffs;
rxnsInfo=p.rxnsInfo;
metPrefix=p.metPrefix;
if isempty(metPrefix)
    metPrefix = 'M_';
end
rxnPrefix=p.rxnPrefix;
if isempty(rxnPrefix)
    rxnPrefix = 'R_';
end
if xor(strcmp(rxnsCoeffs,'none'), strcmp(rxnsInfo,'none'))
    error('Provide both a ''rxnsInfo'' and a ''rxnsCoeffs'' file')
end
newModel=model;

%% Metabolites
if ~strcmp(metsInfo,'none')
    [header,data] = readCurationTable(metsInfo);
    [col,miriamCols] = matchColumns(header,metsInfo,{'metNames','comps'},{ ...
        'metNames',    {'metNames','metName'}; ...
        'comps',       {'comps','comp'}; ...
        'metFormulas', {'formula','formulas','metFormula','metFormulas'}; ...
        'metCharges',  {'charge','charges','metCharge','metCharges'}; ...
        'inchis',      {'inchi','inchis'}; ...
        'metNotes',    {'metNotes','metNote'}});
    metsToAdd.metNames      = data(:,col.metNames);
    metsToAdd.compartments  = data(:,col.comps);
    textFields = {'metFormulas','inchis','metNotes'};
    for f = textFields(isfield(col,textFields))
        metsToAdd.(f{1}) = data(:,col.(f{1}));
    end
    if isfield(col,'metCharges')
        metsToAdd.metCharges = cellfun(@str2num,data(:,col.metCharges),'UniformOutput',false);
        emptyEntry = cellfun(@isempty,metsToAdd.metCharges);
        if all(emptyEntry)
            metsToAdd       = rmfield(metsToAdd,'metCharges');
        else % If some charges are given, assume 0 for those without charges specified
            metsToAdd.metCharges(emptyEntry) = {0};
            metsToAdd.metCharges = cell2mat(metsToAdd.metCharges);
        end
    end
    optFields = intersect({'metFormulas','metCharges','inchis','metNotes'},fieldnames(metsToAdd)','stable');

    % Check if metabolite already exists (check by metName[comp])
    existingMets            = [];
    existingMetsIdx         = [];
    newMets                 = [];
    for i=1:numel(metsToAdd.metNames)
        metIdx              = find(strcmp(newModel.metNames,metsToAdd.metNames{i}));
        existMet            = strcmp(metsToAdd.compartments{i},newModel.comps(newModel.metComps(metIdx)));
        if any(existMet)
            existingMets	= [existingMets, i];
            existingMetsIdx = [existingMetsIdx, metIdx(existMet)];
        else
            newMets	= [newMets, i];
        end
    end

    % Overwrite annotation of existing entries
    if any(existingMets)
        for f = optFields
            newModel = writeExistingField(newModel,f{1},existingMetsIdx,metsToAdd.(f{1})(existingMets));
        end
        newModel = extractAndAddMiriam(newModel,header(miriamCols),data(:,miriamCols),existingMets,existingMetsIdx,'met');
        warning(['The following metabolites are already present in the model, '...
            'their annotation will be overwritten to match the metsInfo file. '...
            'If you do not particularly want to curate their annotations, it '...
            'would be better to removes these metabolites from metsInfo:\n\t\t%s'],...
            strjoin(strcat(metsToAdd.metNames(existingMets),'[',metsToAdd.compartments(existingMets),']'),'\n\t\t'));
    end

    % Continue with new metabolites
    if numel(newMets)>0
        metsToAdd.metNames      = metsToAdd.metNames(newMets);
        metsToAdd.compartments  = metsToAdd.compartments(newMets);
        for f = optFields
            metsToAdd.(f{1}) = metsToAdd.(f{1})(newMets);
            if iscell(metsToAdd.(f{1})) && all(cellfun(@isempty,metsToAdd.(f{1})))
                metsToAdd = rmfield(metsToAdd,f{1});
            end
        end
        % Add metabolites
        newModel        = addMets(newModel,metsToAdd,true,metPrefix);
        addedIdx        = numel(newModel.mets)-numel(newMets)+1:numel(newModel.mets);
        newModel        = extractAndAddMiriam(newModel,header(miriamCols),data(:,miriamCols),newMets,addedIdx,'met');
    end
end
%% Genes
if ~strcmp(genesInfo,'none')
    [header,data] = readCurationTable(genesInfo);
    [col,miriamCols] = matchColumns(header,genesInfo,{'genes'},{ ...
        'genes',          {'genes','gene'}; ...
        'geneShortNames', {'geneShortNames','geneShortName'}});
    genesToAdd.genes            = data(:,col.genes);
    if isfield(col,'geneShortNames')
        genesToAdd.geneShortNames   = data(:,col.geneShortNames);
    end

    [~,notNewGene,existingGene] = intersect(genesToAdd.genes,newModel.genes);
    if ~isempty(notNewGene)
        warning(['The following genes are already present in the model, their annotation '...
            'will be overwritten to match the genesInfo file: \n\t\t%s'],...
            strjoin(genesToAdd.genes(notNewGene),'\n\t\t'));
        if isfield(genesToAdd,'geneShortNames')
            newModel = writeExistingField(newModel,'geneShortNames',existingGene,genesToAdd.geneShortNames(notNewGene));
        end
        newModel = extractAndAddMiriam(newModel,header(miriamCols),data(:,miriamCols),notNewGene,existingGene,'gene');
    end

    % Continue with new genes
    toAdd               = 1:numel(genesToAdd.genes);
    toAdd(notNewGene)   = [];
    if ~isempty(toAdd)
        genesToAdd.genes            = genesToAdd.genes(toAdd);
        if isfield(genesToAdd,'geneShortNames')
            genesToAdd.geneShortNames   = genesToAdd.geneShortNames(toAdd);
            if all(cellfun(@isempty,genesToAdd.geneShortNames))
                genesToAdd = rmfield(genesToAdd,'geneShortNames');
            end
        end
        % Add genes
        newModel = addGenesRaven(newModel,genesToAdd);
        addedIdx = numel(newModel.genes)-numel(toAdd)+1:numel(newModel.genes);
        newModel = extractAndAddMiriam(newModel,header(miriamCols),data(:,miriamCols),toAdd,addedIdx,'gene');
    end
end
%% Reactions
if ~any(strcmp({rxnsCoeffs,rxnsInfo},'none'))
    % Gather all data, first about reaction stoichiometries...
    [header,data] = readCurationTable(rxnsCoeffs);
    col = matchColumns(header,rxnsCoeffs,{'rxnIdx','rxnNames','metNames','comps','coefficient'},{ ...
        'rxnIdx',      {'rxnIdx','index'}; ...
        'rxnNames',    {'rxnNames','rxnName'}; ...
        'metNames',    {'metNames','metName'}; ...
        'comps',       {'comps','comp'}; ...
        'coefficient', {'coefficient','coefficients','coeffs','stoichCoeffs'}});
    rxnCheck.coeffsIdx = str2double(data(:,col.rxnIdx));
    rxns.rxnNames   = data(:,col.rxnNames);
    rxns.metNames   = data(:,col.metNames);
    rxns.comps      = data(:,col.comps);
    rxns.coeffs     = str2double(data(:,col.coefficient));
    rxnCheck.coeffs = strcat(data(:,col.rxnIdx),'***',rxns.rxnNames);

    % ... and then additional reaction-specific data.
    [header,data] = readCurationTable(rxnsInfo);
    [col,miriamCols] = matchColumns(header,rxnsInfo,{'rxnIdx','rxnNames'},{ ...
        'rxnIdx',              {'rxnIdx','index'}; ...
        'rxnNames',            {'rxnNames','rxnName'}; ...
        'grRules',             {'grRules'}; ...
        'lb',                  {'lb'}; ...
        'ub',                  {'ub'}; ...
        'rev',                 {'rev'}; ...
        'subSystems',          {'subSystems'}; ...
        'eccodes',             {'eccodes'}; ...
        'rxnNotes',            {'rxnNotes'}; ...
        'rxnReferences',       {'rxnReferences'}; ...
        'rxnConfidenceScores', {'rxnConfidenceScores'}});
    rxnCheck.rxnsIdx = str2double(data(:,col.rxnIdx));
    rxnCheck.rxns = strcat(data(:,col.rxnIdx),'***',data(:,col.rxnNames));

    notMatching = setxor(rxnCheck.coeffs,rxnCheck.rxns);
    if numel(notMatching)>1
        error(['The following reactions ánd/or their indices are not matched '...
            'between the rxnsInfo and rxnsCoeffs files:\n\t\t%s'],...
            strjoin(regexprep(notMatching,'^\d+***',''),'\n\t\t'))
    end

    rxnsToAdd.rxnNames              = data(:,col.rxnNames);
    for f = intersect({'grRules','subSystems','eccodes','rxnNotes','rxnReferences'},fieldnames(col)','stable')
        rxnsToAdd.(f{1}) = data(:,col.(f{1}));
    end
    for f = intersect({'lb','ub','rev'},fieldnames(col)','stable')
        rxnsToAdd.(f{1}) = cellfun(@str2num,data(:,col.(f{1})));
    end
    if isfield(col,'rxnConfidenceScores')
        rxnsToAdd.rxnConfidenceScores   = cellfun(@str2num,data(:,col.rxnConfidenceScores),'UniformOutput',false);
        emptyEntry = cellfun(@isempty,rxnsToAdd.rxnConfidenceScores);
        if all(emptyEntry)
            rxnsToAdd           = rmfield(rxnsToAdd,'rxnConfidenceScores');
        else % If some rxnConfidenceScores are given, assume 0 for those without rxnConfidenceScores specified
            rxnsToAdd.rxnConfidenceScores(emptyEntry) = {0};
            rxnsToAdd.rxnConfidenceScores = cell2mat(rxnsToAdd.rxnConfidenceScores);
        end
    end
    optFields = intersect({'subSystems','eccodes','rxnNotes','rxnReferences','rxnConfidenceScores'},fieldnames(rxnsToAdd)','stable');

    existingRxn=[];
    notNewRxn=[];
    for i=1:numel(rxnCheck.rxnsIdx)
        rxnRows                     = find(rxnCheck.rxnsIdx(i)==rxnCheck.coeffsIdx);
        rxnsToAdd.mets{i,1}         = cell(1,1);
        rxnsToAdd.stoichCoeffs{i,1} = cell(1,1);
        for j=1:numel(rxnRows)
            newMetComp = [rxns.metNames{rxnRows(j)}, '[', rxns.comps{rxnRows(j)}, ']'];
            try
                newMetComp = getIndexes(newModel,newMetComp,'metcomps');
            catch
                error(['Not all metabolites in reaction "%s" are present in the '...
                    'model or in the provided table with new metabolites.'],rxnsToAdd.rxnNames{i})
            end
            rxnsToAdd.mets{i}(j) = newModel.mets(newMetComp);
            rxnsToAdd.stoichCoeffs{i}{j} = rxns.coeffs(rxnRows(j));
        end
        rxnsToAdd.stoichCoeffs{i}=cell2mat(rxnsToAdd.stoichCoeffs{i});
        %Check if the reaction not already exists, by stoichiometry of its
        %products and reactants
        metIdx              = getIndexes(newModel,rxnsToAdd.mets{i},'mets');
        modelCoeffs         = transpose(newModel.S(metIdx,:));
        modelCoeffs2        = find(~all(modelCoeffs==0,2));
        % Only reactions built from exactly these metabolites can be
        % duplicates. Any model reaction that also involves other
        % metabolites projects onto the same row of coefficients, but is a
        % different reaction and must not be overwritten.
        sameSupport         = full(sum(newModel.S(:,modelCoeffs2)~=0,1))' == numel(metIdx);
        modelCoeffs2        = modelCoeffs2(sameSupport);
        modelCoeffs         = modelCoeffs(modelCoeffs2,:);
        [~, duplicateRxn]   = intersect(modelCoeffs,rxnsToAdd.stoichCoeffs{i},'rows');
        if ~isempty(duplicateRxn)
            existingRxn     = [existingRxn,modelCoeffs2(duplicateRxn)];
            notNewRxn       = [notNewRxn,i];
        end
    end

    % If some reactions already existed, then replace its information
    if ~isempty(existingRxn)
        warning(['The following reactions are the same as existing model '...
            'reactions, their annotation will be overwritten to match the '...
            'rxnCoeffs file: \n\t\t(Existing reaction ID): New reaction name\n\t\t%s'],...
            strjoin(strcat('(',newModel.rxns(existingRxn),{'): '},rxnsToAdd.rxnNames(notNewRxn)),'\n\t\t'));
        newModel.rxnNames(existingRxn) = rxnsToAdd.rxnNames(notNewRxn);
        for f = intersect({'lb','ub','rev'},fieldnames(rxnsToAdd)','stable')
            newModel.(f{1})(existingRxn) = rxnsToAdd.(f{1})(notNewRxn);
        end
        for f = optFields
            values = rxnsToAdd.(f{1})(notNewRxn);
            if strcmp(f{1},'subSystems') && isfield(newModel,'subSystems') ...
                    && any(cellfun(@iscell,newModel.subSystems))
                values = cellfun(@(x) {x}, values, 'UniformOutput', false);
            end
            newModel = writeExistingField(newModel,f{1},existingRxn,values);
        end
        if isfield(rxnsToAdd,'grRules')
            emptyEntries = cellfun(@isempty,rxnsToAdd.grRules(notNewRxn));
            if ~all(emptyEntries)
                newModel    = changeGrRules(newModel,newModel.rxns(existingRxn(~emptyEntries)),rxnsToAdd.grRules(notNewRxn(~emptyEntries)),true);
            end
        end
        newModel    = extractAndAddMiriam(newModel,header(miriamCols),data(:,miriamCols),notNewRxn,existingRxn,'rxn');
    end

    % Continue with new reactions
    toAdd               = 1:numel(rxnsToAdd.rxnNames);
    toAdd(notNewRxn)    = [];
    if ~isempty(toAdd)
        rxnsToAdd.rxnNames      = rxnsToAdd.rxnNames(toAdd);
        for f = intersect({'lb','ub','rev','grRules'},fieldnames(rxnsToAdd)','stable')
            rxnsToAdd.(f{1}) = rxnsToAdd.(f{1})(toAdd);
        end
        for f = optFields
            rxnsToAdd.(f{1}) = rxnsToAdd.(f{1})(toAdd);
            if iscell(rxnsToAdd.(f{1})) && all(cellfun(@isempty,rxnsToAdd.(f{1})))
                rxnsToAdd = rmfield(rxnsToAdd,f{1});
            end
        end
        rxnsToAdd.mets          = rxnsToAdd.mets(toAdd);
        rxnsToAdd.stoichCoeffs  = rxnsToAdd.stoichCoeffs(toAdd);
        rxnsToAdd.rxns          = generateNewIds(model,'rxns',rxnPrefix,numel(rxnsToAdd.rxnNames));

        newModel        = addRxns(newModel,rxnsToAdd,1,[],false,false);
        rxnsModelIdx    = numel(newModel.rxns)-numel(toAdd)+1:numel(newModel.rxns);
        newModel        = extractAndAddMiriam(newModel,header(miriamCols),data(:,miriamCols),toAdd,rxnsModelIdx,'rxn');
    end
end
end

function [header,data] = readCurationTable(fileName)
% Read a tab-delimited curation table into its header row and a cell
% array of text with one column per header field.
fid = fopen(fileName);
if fid==-1
    error('Could not open the file "%s"',fileName);
end
headerLine = fgetl(fid);
if ~ischar(headerLine)
    fclose(fid);
    error('The file "%s" is empty',fileName);
end
headerLine = strrep(headerLine,char(65279),''); % byte order mark
header = strtrim(strsplit(headerLine,'\t'));
header = regexprep(header,'^"(.*)"$','$1');
nCols  = numel(header);
raw    = textscan(fid,repmat('%q',1,nCols),'Delimiter','\t');
fclose(fid);
nRows  = max(cellfun(@numel,raw));
data   = repmat({''},nRows,nCols);
for i=1:nCols
    data(1:numel(raw{i}),i) = raw{i};
end
end

function [col,miriamCols] = matchColumns(header,fileName,required,spec)
% Locate the core columns of a curation table by header name (case
% insensitive). spec has one row per core field: {field, accepted headers}.
% Every other column is a MIRIAM namespace.
col        = struct();
isCore     = false(1,numel(header));
for i=1:size(spec,1)
    match = find(ismember(lower(header),lower(spec{i,2})));
    if ~isempty(match)
        col.(spec{i,1}) = match(1);
        isCore(match)   = true;
    end
end
missing = setdiff(required,fieldnames(col)','stable');
if ~isempty(missing)
    error('The file "%s" has no column with header: %s',fileName,strjoin(missing,', '));
end
miriamCols = find(~isCore & ~cellfun(@isempty,header));
end

function model = writeExistingField(model,field,modelIdx,values)
% Write the values of an optional field onto existing entities. A field
% the model does not have is created from its ravenModelFields default,
% unless all supplied values are empty.
if ~isfield(model,field)
    if iscell(values) && all(cellfun(@isempty,values))
        return
    end
    reg   = ravenModelFields();
    entry = reg(strcmp({reg.name},field));
    n     = numel(model.([entry.type 's']));
    if iscell(values)
        model.(field) = repmat({entry.default},n,1);
    else
        model.(field) = repmat(entry.default,n,1);
    end
end
model.(field)(modelIdx) = values;
end

function newModel = extractAndAddMiriam(model,miriamName,miriamData,inputIndex,modelIndex,type)
newModel=model;
miriamValues    = miriamData(inputIndex,:);
emptyMiriam     = all(cellfun(@isempty,miriamValues),1);
miriamName(emptyMiriam)     = [];
miriamValues(:,emptyMiriam) = [];
if ~isfield(newModel,[type 'Miriams'])
    newModel.([type 'Miriams'])=cell(numel(newModel.([type 's'])),1);
end
for i=1:numel(miriamName)
    newModel = editMiriam(newModel,type,modelIndex,miriamName{i},miriamValues(:,i),'replace');
end
end
