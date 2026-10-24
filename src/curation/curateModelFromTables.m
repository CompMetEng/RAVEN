function newModel=curateModelFromTables(model,metsInfo,varargin)
% curateModelFromTables  Curate or add mets, rxns and genes from tables.
%
% Curate existing and/or add new metabolites, reactions and genes from
% tabular data files.
%
% If the *.tsv files contain metabolites, reactions and/or genes that are
% already present in the model, then information in the model will be
% overwritten. Note that this includes empty annotations in the *.tsv
% files, unless 'overwriteWithEmpty' is false! Metabolites are matched by
% metaboliteName[comp]; reactions by the stoichiometry of its reactants and
% products; genes by their gene name.
% This function can therefore be used to add new entities in the model, or
% curate those already existing in the model.
%
% Columns are identified by their header (case-insensitive), so their order
% does not matter and optional columns can be left out. Every field that a
% table supplies is written to new and existing entities alike; a field the
% model does not have yet is created, with its default value for the other
% entities, unless all supplied values for it are empty. An empty grRules
% cell leaves the rule of an existing reaction unchanged, and so does an
% empty lb, ub or rev cell: an empty bound means that no bound is given.
% A new reaction without a lb gets a lower bound of -Inf, and one without an
% ub gets model.annotation.defaultUB (Inf if the model has none).
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
% overwriteWithEmpty : logical
%     if true, an empty cell in a row of an existing metabolite, gene or
%     reaction empties the corresponding field of the model. If false, an
%     empty cell leaves the existing value unchanged and only the cells with a
%     value are written (default true; [] also means the default). Must be a
%     scalar logical or 0/1. A cell with only blanks counts as empty in both
%     modes. Empty grRules cells never overwrite. An empty or NaN lb, ub or
%     rev cell always keeps the existing value. For the charge and the
%     confidence score, an empty or NaN cell is treated alike: when
%     overwriteWithEmpty is true and at least one other row of the column has
%     a value, it becomes 0; a column in which no row has a value is not
%     written, and the existing values stay. The same holds for a MIRIAM
%     column: a column in which no row has a value does not clear the
%     annotations of existing entities.
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
%
% Cell values: the charge, confidence score, lb, ub, rev, rxnIdx and
% coefficient cells are read as plain decimal numbers (digits with an
% optional sign, decimal point and exponent), Inf or NaN in any case. Only
% the lb, ub and rev cells also accept true and false in any case (read as
% 1 and 0). Any other text, including text with a comma, hexadecimal numbers
% and expressions, raises RAVEN:badInput; no expression is evaluated. A
% charge must be a finite whole number (or empty or NaN). A rev must be 0, 1,
% true or false (or empty or NaN). The rev column sets the reversibility of
% existing reactions only, independently of lb and ub, which are written as
% given; for new reactions the reversibility is derived from lb and ub.
%
% Gene rules: a grRules cell must be a well-formed gene rule (genes joined by
% and/or in any case, or by &, &&, | and ||, with balanced parentheses); a
% rule with an operator that lacks an operand, empty parentheses or two genes
% without an operator between them raises RAVEN:badInput before the model is
% changed. Whether the genes exist is checked when the rule is applied.
%
% Genes: a gene id that the genesInfo table lists more than once raises
% RAVEN:badInput, whether the model has the gene or not.
%
% Bounds: a lower bound above the upper bound is written as given. solveLP
% rejects such a model when it is run ("Invalid set of bounds for
% reaction(s): ...").
%
% Reaction tables: a row of the rxnsCoeffs table without a rxnIdx, rxnNames,
% metNames and comps value is ignored; a row with any of them needs a
% rxnIdx and a coefficient that is a finite number other than 0. A row of
% the rxnsInfo table without any value is ignored; a row with a value needs
% a rxnIdx. A rxnIdx is a whole number of at least 1; the indices (and
% names) of the two reaction tables must match exactly, compared by their
% numeric value (1, 1.0 and 1e0 are the same index). Otherwise RAVEN:badInput
% is raised.

p=parseRAVENargs(varargin, {'genesInfo','none'; 'rxnsCoeffs','none'; 'rxnsInfo','none'; 'metPrefix',[]; 'rxnPrefix',[]; 'overwriteWithEmpty',true});
overwriteWithEmpty=p.overwriteWithEmpty;
if isempty(overwriteWithEmpty)
    overwriteWithEmpty=true;
end
if ~((islogical(overwriteWithEmpty) || isnumeric(overwriteWithEmpty)) && isscalar(overwriteWithEmpty) ...
        && isreal(overwriteWithEmpty) && ~isnan(overwriteWithEmpty) && ismember(overwriteWithEmpty,[0 1]))
    error('RAVEN:badInput','overwriteWithEmpty must be a scalar logical or 0/1');
end
overwriteWithEmpty=logical(overwriteWithEmpty);
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
    error('RAVEN:badInput','Provide both a ''rxnsInfo'' and a ''rxnsCoeffs'' file')
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
        metsToAdd.metCharges = cellfun(@(x) parseChargeCell(x,'metsInfo'),data(:,col.metCharges),'UniformOutput',false);
        isNanCharge = cellfun(@(x) isscalar(x) && isnan(x),metsToAdd.metCharges);
        metsToAdd.metCharges(isNanCharge) = {[]}; % a NaN charge counts as not given
        emptyEntry = cellfun(@isempty,metsToAdd.metCharges);
        emptyCharge = emptyEntry;
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
            skip = [];
            if ~overwriteWithEmpty
                if strcmp(f{1},'metCharges')
                    skip = emptyCharge(existingMets);
                else
                    skip = emptyCells(metsToAdd.(f{1})(existingMets));
                end
            end
            newModel = writeExistingField(newModel,f{1},existingMetsIdx,metsToAdd.(f{1})(existingMets),skip);
        end
        newModel = extractAndAddMiriam(newModel,header(miriamCols),data(:,miriamCols),existingMets,existingMetsIdx,'met',~overwriteWithEmpty);
        warning(['The following metabolites are already present in the model, '...
            'their annotation will be overwritten to match the metsInfo file. '...
            'If you do not particularly want to curate their annotations, it '...
            'would be better to remove these metabolites from metsInfo:\n\t\t%s'],...
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
    checkGenesListedOnce(genesToAdd.genes);
    if isfield(col,'geneShortNames')
        genesToAdd.geneShortNames   = data(:,col.geneShortNames);
    end

    [~,notNewGene,existingGene] = intersect(genesToAdd.genes,newModel.genes);
    if ~isempty(notNewGene)
        warning(['The following genes are already present in the model, their annotation '...
            'will be overwritten to match the genesInfo file: \n\t\t%s'],...
            strjoin(genesToAdd.genes(notNewGene),'\n\t\t'));
        if isfield(genesToAdd,'geneShortNames')
            skip = [];
            if ~overwriteWithEmpty
                skip = emptyCells(genesToAdd.geneShortNames(notNewGene));
            end
            newModel = writeExistingField(newModel,'geneShortNames',existingGene,genesToAdd.geneShortNames(notNewGene),skip);
        end
        newModel = extractAndAddMiriam(newModel,header(miriamCols),data(:,miriamCols),notNewGene,existingGene,'gene',~overwriteWithEmpty);
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
    % A row without a rxnIdx, rxnNames, metNames and comps value is ignored
    rowNo = (1:size(data,1))';
    keyCells = data(:,[col.rxnIdx,col.rxnNames,col.metNames,col.comps]);
    ignoredRow = all(reshape(emptyCells(keyCells(:)),size(keyCells)),2);
    data(ignoredRow,:) = [];
    rowNo(ignoredRow) = [];
    rxnCheck.coeffsIdx = parseIndexCells(data(:,col.rxnIdx),'rxnIdx','rxnsCoeffs');
    rxns.rxnNames   = data(:,col.rxnNames);
    rxns.metNames   = data(:,col.metNames);
    rxns.comps      = data(:,col.comps);
    rxns.coeffs     = parseNumericCells(data(:,col.coefficient),'coefficient','rxnsCoeffs');
    checkCoefficientsGiven(rxnCheck.coeffsIdx,rxns.coeffs,rowNo);
    rxnCheck.coeffs = indexKeys(rxnCheck.coeffsIdx,rxns.rxnNames);

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
    % A row without any value is ignored; a row with a value needs a rxnIdx
    rowNo = (1:size(data,1))';
    ignoredRow = all(reshape(emptyCells(data(:)),size(data)),2);
    data(ignoredRow,:) = [];
    rowNo(ignoredRow) = [];
    rxnCheck.rxnsIdx = parseIndexCells(data(:,col.rxnIdx),'rxnIdx','rxnsInfo');
    noIdx = find(isnan(rxnCheck.rxnsIdx),1);
    if ~isempty(noIdx)
        error('RAVEN:badInput','Data row %d of the rxnsInfo file has values but its rxnIdx is empty or NaN',rowNo(noIdx));
    end
    rxnCheck.rxns = indexKeys(rxnCheck.rxnsIdx,data(:,col.rxnNames));

    notMatching = setxor(rxnCheck.coeffs,rxnCheck.rxns);
    if ~isempty(notMatching)
        error('RAVEN:badInput',['The following reactions and/or their indices are not matched '...
            'between the rxnsInfo and rxnsCoeffs files:\n\t\t%s'],...
            strjoin(regexprep(notMatching,'^(\d+)\*\*\*(.*)$','index $1 ($2)'),'\n\t\t'))
    end

    rxnsToAdd.rxnNames              = data(:,col.rxnNames);
    for f = intersect({'grRules','subSystems','eccodes','rxnNotes','rxnReferences'},fieldnames(col)','stable')
        rxnsToAdd.(f{1}) = data(:,col.(f{1}));
    end
    if isfield(col,'grRules')
        checkGrRules(rxnsToAdd.grRules,rxnsToAdd.rxnNames,rowNo);
    end
    for f = intersect({'lb','ub','rev'},fieldnames(col)','stable')
        rxnsToAdd.(f{1}) = parseNumericCells(data(:,col.(f{1})),f{1},'rxnsInfo');
    end
    if isfield(rxnsToAdd,'rev')
        checkRevValues(rxnsToAdd.rev,data(:,col.rev));
    end
    if isfield(col,'rxnConfidenceScores')
        rxnsToAdd.rxnConfidenceScores   = cellfun(@(x) parseNumberCell(x,'rxnConfidenceScores','rxnsInfo'),data(:,col.rxnConfidenceScores),'UniformOutput',false);
        isNanScore = cellfun(@(x) isscalar(x) && isnan(x),rxnsToAdd.rxnConfidenceScores);
        rxnsToAdd.rxnConfidenceScores(isNanScore) = {[]}; % a NaN score counts as not given
        emptyEntry = cellfun(@isempty,rxnsToAdd.rxnConfidenceScores);
        emptyScore = emptyEntry;
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
                error('RAVEN:badInput',['Not all metabolites in reaction "%s" are present in the '...
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
        nameValues = rxnsToAdd.rxnNames(notNewRxn);
        nameIdx    = existingRxn;
        if ~overwriteWithEmpty
            nameIdx    = nameIdx(~emptyCells(nameValues));
            nameValues = nameValues(~emptyCells(nameValues));
        end
        newModel.rxnNames(nameIdx) = nameValues;
        for f = intersect({'lb','ub','rev'},fieldnames(rxnsToAdd)','stable')
            values = rxnsToAdd.(f{1})(notNewRxn);
            given  = ~isnan(values);
            newModel.(f{1})(existingRxn(given)) = values(given);
        end
        for f = optFields
            values = rxnsToAdd.(f{1})(notNewRxn);
            skip = [];
            if ~overwriteWithEmpty
                if strcmp(f{1},'rxnConfidenceScores')
                    skip = emptyScore(notNewRxn);
                else
                    skip = emptyCells(values);
                end
            end
            if strcmp(f{1},'subSystems') && isfield(newModel,'subSystems') ...
                    && any(cellfun(@iscell,newModel.subSystems))
                values = cellfun(@(x) {x}, values, 'UniformOutput', false);
            end
            newModel = writeExistingField(newModel,f{1},existingRxn,values,skip);
        end
        if isfield(rxnsToAdd,'grRules')
            emptyEntries = emptyCells(rxnsToAdd.grRules(notNewRxn));
            if ~all(emptyEntries)
                newModel    = changeGrRules(newModel,newModel.rxns(existingRxn(~emptyEntries)),rxnsToAdd.grRules(notNewRxn(~emptyEntries)),true);
            end
        end
        newModel    = extractAndAddMiriam(newModel,header(miriamCols),data(:,miriamCols),notNewRxn,existingRxn,'rxn',~overwriteWithEmpty);
    end

    % Continue with new reactions
    toAdd               = 1:numel(rxnsToAdd.rxnNames);
    toAdd(notNewRxn)    = [];
    if ~isempty(toAdd)
        rxnsToAdd.rxnNames      = rxnsToAdd.rxnNames(toAdd);
        for f = intersect({'lb','ub','rev','grRules'},fieldnames(rxnsToAdd)','stable')
            rxnsToAdd.(f{1}) = rxnsToAdd.(f{1})(toAdd);
        end
        if isfield(rxnsToAdd,'lb')
            rxnsToAdd.lb(isnan(rxnsToAdd.lb)) = -inf;
        end
        if isfield(rxnsToAdd,'ub')
            defaultUB = inf;
            if isfield(newModel,'annotation') && isfield(newModel.annotation,'defaultUB')
                defaultUB = newModel.annotation.defaultUB;
            end
            rxnsToAdd.ub(isnan(rxnsToAdd.ub)) = defaultUB;
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
    error('RAVEN:badInput','Could not open the file "%s"',fileName);
end
headerLine = fgetl(fid);
if ~ischar(headerLine)
    fclose(fid);
    error('RAVEN:badInput','The file "%s" is empty',fileName);
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
blank = cellfun(@(x) ischar(x) && isempty(strtrim(x)),data);
data(blank) = {''};
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
    error('RAVEN:badInput','The file "%s" has no column with header: %s',fileName,strjoin(missing,', '));
end
miriamCols = find(~isCore & ~cellfun(@isempty,header));
end

function model = writeExistingField(model,field,modelIdx,values,skip)
% Write the values of an optional field onto existing entities. A field
% the model does not have is created from its ravenModelFields default,
% unless all supplied values are empty. Entries flagged in the logical
% vector skip are not written.
if nargin>4 && ~isempty(skip)
    modelIdx = modelIdx(~skip(:)');
    values   = values(~skip(:));
    if isempty(modelIdx)
        return
    end
end
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

function newModel = extractAndAddMiriam(model,miriamName,miriamData,inputIndex,modelIndex,type,skipEmpty)
% Replace the MIRIAM annotations of the entities in modelIndex with the
% values of the table rows in inputIndex. With skipEmpty, an empty cell
% leaves the annotation of that entity unchanged.
newModel=model;
miriamValues    = miriamData(inputIndex,:);
if nargin>6 && skipEmpty
    for i=1:numel(miriamName)
        given = ~emptyCells(miriamValues(:,i));
        if any(given)
            if ~isfield(newModel,[type 'Miriams'])
                newModel.([type 'Miriams'])=cell(numel(newModel.([type 's'])),1);
            end
            idx = modelIndex(:);
            newModel = editMiriam(newModel,type,idx(given),miriamName{i},miriamValues(given,i),'replace');
        end
    end
    return
end
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

function empty = emptyCells(values)
% Logical column vector that is true for the empty entries of a cell array
% of text (a cell with only blanks counts as empty) or of a numeric cell.
empty = cellfun(@(x) isempty(x) || (ischar(x) && isempty(strtrim(x))),values(:));
end

function values = parseIndexCells(cells,field,fileLabel)
% Reaction indices from table cells: NaN where the cell is empty or says NaN
% (the row is ignored), otherwise a finite whole number of at least 1.
values = parseNumericCells(cells,field,fileLabel);
bad = find(~isnan(values) & (~isfinite(values) | values~=round(values) | values<1),1);
if ~isempty(bad)
    error('RAVEN:badInput','The %s value "%s" in the %s file is not a whole number of at least 1',field,strtrim(cells{bad}),fileLabel);
end
end

function values = parseNumericCells(cells,field,fileLabel)
% Numbers from table cells, NaN where the cell is empty or says NaN. A cell
% that is not a number (or true/false) raises an error.
values = nan(numel(cells),1);
for i=1:numel(cells)
    v = parseNumberCell(cells{i},field,fileLabel);
    if ~isempty(v)
        values(i) = v;
    end
end
end

function v = parseChargeCell(text,fileLabel)
% Charge from a table cell: [] for an empty cell, NaN for NaN, otherwise a
% finite whole number. A fractional or infinite charge raises an error.
v = parseNumberCell(text,'charge',fileLabel);
if ~isempty(v) && ~isnan(v) && (~isfinite(v) || v~=round(v))
    error('RAVEN:badInput','The charge value "%s" in the %s file is not a finite whole number',strtrim(text),fileLabel);
end
end

function checkCoefficientsGiven(idx,coeffs,rowNo)
% Every row of the rxnsCoeffs table that is kept has a rxnIdx and a
% coefficient that is a finite number other than 0.
noIdx = find(isnan(idx),1);
if ~isempty(noIdx)
    error('RAVEN:badInput',['Data row %d of the rxnsCoeffs file has values but its rxnIdx is empty or NaN'],rowNo(noIdx));
end
missing = find(isnan(coeffs),1);
if ~isempty(missing)
    error('RAVEN:badInput',['The coefficient in data row %d of the rxnsCoeffs file is empty or NaN; ' ...
        'every row with a reaction, metabolite or compartment needs a coefficient'],rowNo(missing));
end
bad = find(~isfinite(coeffs) | coeffs==0,1);
if ~isempty(bad)
    error('RAVEN:badInput','The coefficient in data row %d of the rxnsCoeffs file must be a finite non-zero number',rowNo(bad));
end
end

function keys = indexKeys(idx,names)
% One key per table row from the numeric reaction index and the name, so that
% the indices of two tables are compared by value.
keys = cell(numel(idx),1);
for i=1:numel(idx)
    keys{i} = sprintf('%.0f***%s',idx(i),names{i});
end
end

function checkGenesListedOnce(genes)
% Raise RAVEN:badInput for a gene id that the genesInfo table lists more than
% once, whether the model has the gene or not: two rows for one gene would
% say two things about it, and the second row would add a second gene with
% the same id.
[ids,~,j] = unique(genes);
repeated = ids(accumarray(j,1) > 1);
if ~isempty(repeated)
    error('RAVEN:badInput',['The following genes are listed more than once '...
        'in the genesInfo file:\n\t\t%s'],strjoin(repeated,'\n\t\t'));
end
end


function checkGrRules(rules,rxnNames,rowNo)
% Raise RAVEN:badInput for a grRules cell that is not a well-formed gene rule
% (unbalanced parentheses, an operator without two operands, empty
% parentheses, two genes without an operator between them). Empty cells are
% not checked. The grammar is that of parseGrRule.
for i = find(~emptyCells(rules))'
    try
        parseGrRule(rules{i});
    catch err
        if strcmp(err.identifier,'RAVEN:badGrRule')
            error('RAVEN:badInput','The grRules value "%s" of reaction "%s" (data row %d of the rxnsInfo file) is not a valid gene rule: %s',...
                rules{i},rxnNames{i},rowNo(i),err.message);
        end
        rethrow(err);
    end
end
end


function checkRevValues(values,cells)
% A rev cell is 0, 1, true or false; an empty or NaN cell is not given.
bad = find(~isnan(values) & values~=0 & values~=1,1);
if ~isempty(bad)
    error('RAVEN:badInput','The rev value "%s" in the rxnsInfo file is not 0, 1, true or false',strtrim(cells{bad}));
end
end

function v = parseNumberCell(text,field,fileLabel)
% Real number from a table cell: [] for an empty cell, 1 or 0 for true or
% false in any case (lb, ub and rev cells only), NaN or +/-Inf for those words in any case, otherwise a
% plain decimal number (digits with an optional sign, decimal point and
% exponent). Any other text raises an error; nothing is evaluated.
text = strtrim(text);
if isempty(text)
    v = [];
    return
end
word = lower(text);
if any(strcmp(word,{'true','false'})) && ~ismember(field,{'lb','ub','rev'})
    error('RAVEN:badInput','The %s value "%s" in the %s file is not a number',field,text,fileLabel);
end
if strcmp(word,'true')
    v = 1;
elseif strcmp(word,'false')
    v = 0;
elseif ~isempty(regexp(word,'^[+-]?nan$','once'))
    v = NaN;
elseif ~isempty(regexp(word,'^[+-]?inf$','once'))
    v = Inf;
    if word(1)=='-'
        v = -Inf;
    end
elseif ~isempty(regexp(text,'^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$','once'))
    v = str2double(text);
else
    error('RAVEN:badInput','The %s value "%s" in the %s file is not a number',field,text,fileLabel);
end
end
