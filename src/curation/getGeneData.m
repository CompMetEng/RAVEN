function geneTable = getGeneData(accession, varargin)
% getGeneData  Build a gene-mapping table from NCBI genome annotation files.
%
% Parses the GFF3 annotation file of an NCBI genome assembly and produces
% a table mapping locus tags to gene symbols, protein IDs and other
% identifiers. The table can be passed to renameModelGenes to change gene
% identifiers in a RAVEN model. If the input is an assembly accession,
% downloadGenomeData downloads the GFF3 and protein FASTA (.faa) files to
% downloadDir, unless both are already present there.
%
% Parameters
% ----------
% accession : char
%     NCBI genome assembly accession, for example 'GCF_000002595.2', or the
%     path to a local GFF3 annotation file. An error is raised if the
%     input is neither an existing file nor an accession of the form
%     'GCF_' or 'GCA_' followed by digits and an optional version.
%
% Name-Value Arguments
% --------------------
% outputFile : char
%     Path of the tab-delimited .tsv file to which the table is written.
%     No file is written if omitted or empty (default '').
% downloadDir : char
%     Directory where genome files are downloaded if not already present.
%     Only used when the input is an accession (default: current working
%     directory).
%
% Returns
% -------
% geneTable : table
%     MATLAB table with one row per distinct combination of gene and
%     protein accession (a CDS feature with a protein accession, repeated
%     CDS rows are merged), with the columns below. A column entry is empty
%     when the annotation lacks the value.
%
%     - locus_tag : locus identifier (for example 'Cre01.g000001')
%     - old_locus_tag : earlier locus identifier
%     - GeneID : NCBI Gene ID (for example '5723799')
%     - gene_name : gene symbol (for example 'rbcL')
%     - GenBank_protein : protein accession (for example 'XP_001698190.2'),
%       as in the protein FASTA headers
%     - UniProt : UniProt accession from the Dbxref attribute of the gene
%
% Examples
% --------
%     % Fetch and parse gene data for Chlamydomonas reinhardtii
%     geneTable = getGeneData('GCF_000002595.2');
%
%     % Save the mapping table to a specific file
%     geneTable = getGeneData('GCF_000002595.2', 'chlamy.tsv');
%
%     % Save to a file and download to a specific directory
%     geneTable = getGeneData('GCF_000002595.2', 'chlamy.tsv', 'data/');
%
% Notes
% -----
% When outputFile is provided the table is written as UTF-8 encoded
% tab-delimited text with a header row, which readtable and
% renameModelGenes can read. A CDS whose gene has no locus_tag, or whose
% gene feature is missing, gives a row with an empty locus_tag.

	if nargin < 1 || isempty(accession)
        error('getGeneData:missingInput', ...
            'accession is required. Provide an NCBI assembly accession, e.g. ''GCF_000002595.2''.');
    end

    p = parseRAVENargs(varargin, {'outputFile', ''; 'downloadDir', pwd});
    outputFile = char(p.outputFile);
    downloadDir = char(p.downloadDir);

    % Resolve input → local GFF3 path. Decide whether input is a file path or an NCBI accession.
    if ischar(accession) || isstring(accession)
        inp = char(accession);
    else
        error('getGeneData:badInputType', 'Input must be a character vector or string.');
    end
    
    gffPath = '';
    % If input is an existing file, use it.
    if isfile(inp)
        gffPath = inp;
    else
        % Check if it looks like an NCBI assembly accession (GCF_... or GCA_...)
        if ~isempty(regexp(inp, '^GC[FA]_\d+(\.\d+)?$', 'once'))
            gffPath = downloadGenomeData(inp, downloadDir);
        else
            error('getGeneData:badInput', ...
                'Input "%s" is neither an existing file nor a recognised NCBI assembly accession.', inp);
        end
    end

    % Parse the GFF3 annotation into a gene-mapping table
    geneTable = parseGFF(gffPath);

    % Save result only if an output file was requested
    if ~isempty(outputFile)
        writetable(geneTable, outputFile, 'Delimiter', '\t', 'FileType', 'text');
        fprintf('Gene table saved to: %s  (%d rows)\n', outputFile, height(geneTable));
    end
end

%--------------------------------------------------------------------------
% Helper functions

function T = parseGFF(gffPath)
% Parse a GFF3 file extracting gene + CDS pairs. The protein accession is
% taken from each CDS protein_id, and the owning gene is resolved through
% the Parent chain (CDS → mRNA → gene for eukaryotes, CDS → gene for
% prokaryotes), so the GenBank_protein column matches the protein FASTA.
%
% Returns a table with columns:
%   locus_tag | old_locus_tag | GeneID | gene_name | GenBank_protein | UniProt

    batchSize = 10000;
    rows  = cell(batchSize, 6);
    nRows = 0;
    % Use a containers.Map to hold the last gene record keyed by gene ID,
    % and an mRNA → gene map to resolve eukaryote CDS parents.
    geneMap    = containers.Map('KeyType', 'char', 'ValueType', 'any');
    mrnaToGene = containers.Map('KeyType', 'char', 'ValueType', 'char');
    emptyGene = struct('locus_tag','','old_locus_tag','','geneID','','name','','uniProt','');

    fid = fopen(gffPath, 'r');
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || isempty(line) || line(1) == '#', continue, end
        fields = strsplit(line, '\t');
        if numel(fields) < 9, continue, end

        featureType = fields{3};
        attrs       = fields{9};

        if strcmp(featureType, 'gene')
            g = parseGeneAttrs(attrs);
            if ~isempty(g.id)
                geneMap(g.id) = g;
            end

        elseif strcmp(featureType, 'mRNA')
            id     = extractAttr(attrs, 'ID');
            parent = extractAttr(attrs, 'Parent');
            if ~isempty(id) && ~isempty(parent)
                mrnaToGene(id) = parent;
            end

        elseif strcmp(featureType, 'CDS')
            child     = parseCDSAttrs(attrs);
            proteinID = child.protein_id;
            if isempty(child.parent) || isempty(proteinID), continue, end

            % Resolve the owning gene through the Parent chain
            if isKey(mrnaToGene, child.parent)
                geneKey = mrnaToGene(child.parent);
            else
                geneKey = child.parent;  % prokaryote: CDS Parent is the gene
            end
            if isKey(geneMap, geneKey)
                g = geneMap(geneKey);
            else
                g = emptyGene;
            end

            nRows = nRows + 1;
            if nRows > size(rows, 1)
                rows = [rows; cell(batchSize, 6)];
            end
            rows(nRows, :) = {g.locus_tag, g.old_locus_tag, ...
                               g.geneID,    g.name,          proteinID, g.uniProt};
        end
    end
    fclose(fid);

    rows = rows(1:nRows, :);

    % Deduplicate rows (CDS repeats per exon).
    rowKeys = cellfun(@(a,b,c,d,e) strjoin({a,b,c,d,e}, '|'), ...
        rows(:,1), rows(:,2), rows(:,3), rows(:,4), rows(:,5), ...
        'UniformOutput', false);
    [~, ia] = unique(rowKeys, 'stable');
    rows = rows(ia, :);

    T = cell2table(rows, 'VariableNames', ...
        {'locus_tag','old_locus_tag','GeneID','gene_name','GenBank_protein','UniProt'});
end

% Attibute parsing utilities
function value = extractAttr(attrStr, key)
% Extract the value of a GFF3 attribute by key.
% Handles URL-encoded characters (%2C etc.) for robustness.
% Returns '' if key is not present.
    pattern = sprintf('(?:^|;)%s=([^;]+)', regexptranslate('escape', key));
    tok = regexp(attrStr, pattern, 'tokens', 'once');
    if isempty(tok)
        value = '';
    else
        value = gff3Decode(tok{1});
    end
end

function out = gff3Decode(str)
% Decode GFF3 percent-encoded characters (%XX). Unlike urldecode, '+' is
% left untouched, as GFF3 does not use it to encode spaces.
    out = regexprep(str, '%([0-9A-Fa-f]{2})', '${char(hex2dec($1))}');
end

function value = extractDbxrefField(dbxref, prefix)
% From a comma-separated Dbxref value like
%   "Phytozome:Cre16.g651050,GeneID:5723799,GenBank:XM_001698190.2"
% extract the part after "prefix:".
% Returns '' if not found.
    if isempty(dbxref)
        value = '';
        return
    end
    parts   = strsplit(dbxref, ',');
    pattern = sprintf('^%s:(.+)$', regexptranslate('escape', prefix));
    value   = '';
    for i = 1:numel(parts)
        tok = regexp(strtrim(parts{i}), pattern, 'tokens', 'once');
        if ~isempty(tok)
            value = tok{1};
            return
        end
    end
end

function g = parseGeneAttrs(attrs)
    g.id            = extractAttr(attrs, 'ID');
    g.name          = extractAttr(attrs, 'Name');
    g.locus_tag     = extractAttr(attrs, 'locus_tag');
    g.old_locus_tag = extractAttr(attrs, 'old_locus_tag');

    % Dbxref is a comma-separated list like:
    %   GeneID:5723799,GenBank:XM_001698190.2
    dbxref          = extractAttr(attrs, 'Dbxref');
    g.geneID        = extractDbxrefField(dbxref, 'GeneID');

    % UniProt accession appears as UniProtKB/Swiss-Prot:P12345 or UniProtKB:P12345
    g.uniProt       = extractDbxrefField(dbxref, 'UniProtKB/Swiss-Prot');
    if isempty(g.uniProt)
        g.uniProt   = extractDbxrefField(dbxref, 'UniProtKB');
    end
end

function c = parseCDSAttrs(attrs)
    c.parent     = extractAttr(attrs, 'Parent');
    c.protein_id = extractAttr(attrs, 'protein_id');
    if isempty(c.protein_id)
        % Some files put it in Dbxref as GenBank:WP_...
        dbxref       = extractAttr(attrs, 'Dbxref');
        c.protein_id = extractDbxrefField(dbxref, 'GenBank');
    end
end