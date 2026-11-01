function GSS = getUniProtScores(organism, varargin)
% getUniProtScores  Fetch curated subcellular locations from UniProt as a gene scoring structure.
%
% Queries the UniProtKB REST API for the "Subcellular location" annotations of an organism and
% parses them into the GSS consumed by predictLocalization. No manual export is needed. UniProt
% locations are qualitative, so each annotated compartment gets score 1.0 and a protein with
% several annotated locations scores 1.0 in each of them. The function requires an internet
% connection.
%
% Parameters
% ----------
% organism : double or char
%     a UniProt organism (taxonomy) id (for example 559292 for S. cerevisiae S288C).
%
% Name-Value Arguments
% --------------------
% idField : char
%     which identifier becomes the gene id: 'gene_oln' (ordered locus or ORF name, the default,
%     which matches yeast-GEM), 'accession', or 'gene_primary'.
% reviewed : logical
%     true to restrict the query to curated Swiss-Prot entries (default true).
% compartmentMap : containers.Map
%     maps UniProt location terms to the compartment ids of the model (default
%     defaultCompartmentMap). The keys are also the vocabulary searched for in the
%     "Subcellular location [CC]" free text of UniProt, matched as whole words and
%     case insensitively. The keys must be the noun forms used by UniProt
%     ('endoplasmic reticulum', not 'ER'), and a term that is not a key is not searched.
%     Several keys may share one compartment id, which merges them; since every hit
%     scores 1.0, a protein annotated to any of the merged terms scores 1.0 in that
%     compartment. One key cannot map to two compartment ids; parseScores describes
%     how to duplicate a column of the returned GSS instead.
% extraQuery : char
%     an additional UniProt query clause, combined with the query by AND (for example 'gene:CIT1').
%
% Returns
% -------
% GSS : struct
%     a gene scoring structure for predictLocalization.
%
% Examples
% --------
%     GSS = getUniProtScores(559292);                       % all reviewed S. cerevisiae proteins
%     GSS = getUniProtScores(559292, 'idField', 'accession');
%     map = defaultCompartmentMap; map('cell wall') = 'ce';  % extend the vocabulary
%     GSS = getUniProtScores(559292, 'compartmentMap', map);
%
% See also
% --------
% parseScores, predictLocalization, defaultCompartmentMap

p = parseRAVENargs(varargin, {'idField','gene_oln'; 'reviewed',true; ...
                              'compartmentMap',[]; 'extraQuery',[]});
idField = char(p.idField);
if ~any(strcmp(idField, {'gene_oln','accession','gene_primary'}))
    error('RAVEN:badInput', 'idField must be gene_oln, accession, or gene_primary');
end

query = ['organism_id:' char(string(organism))];
if p.reviewed; query = [query ' AND reviewed:true']; end
if ~isempty(p.extraQuery); query = [query ' AND (' char(p.extraQuery) ')']; end

% Put the chosen id field first so parseScores reads it as the gene id column.
fields = [idField ',cc_subcellular_location'];
outFile = [tempname '.tsv'];
cleanup = onCleanup(@() i_deleteIfExists(outFile)); %#ok<NASGU> — removes temp file on return/error
websave(outFile, 'https://rest.uniprot.org/uniprotkb/stream', ...
        'query', query, 'format', 'tsv', 'fields', fields, weboptions('Timeout', 120));

% The UniProt /stream endpoint gzip-compresses the body and MATLAB does not auto-decompress it,
% so detect the gzip magic bytes and inflate before parsing.
fid = fopen(outFile, 'r');
magic = fread(fid, 2, 'uint8')';
fclose(fid);
if isequal(magic, [31 139])
    gz = [outFile '.gz'];
    movefile(outFile, gz);
    gunzip(gz);                 % writes outFile (drops the .gz extension)
    delete(gz);
end

map = p.compartmentMap;
if isempty(map); map = defaultCompartmentMap(); end
GSS = parseScores(outFile, 'predictor', 'uniprot', 'compartmentMap', map);
end

function i_deleteIfExists(f)
if isfile(f); delete(f); end
end
