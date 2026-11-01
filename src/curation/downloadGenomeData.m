function [gffFile, faaFile] = downloadGenomeData(accession, varargin)
% downloadGenomeData  Download genome annotation files from NCBI.
%
% Downloads the GFF3 annotation and protein FASTA (.faa) files for an NCBI
% genome assembly accession using the NCBI Datasets v2 API, which returns
% both files in a single archive. The files are saved as
% <accession>_genomic.gff and <accession>_protein.faa in outputDir. They
% are input files for getGeneData.
%
% Parameters
% ----------
% accession : char
%     NCBI genome assembly accession, for example 'GCF_000002595.2'. The
%     prefix is 'GCF_' (RefSeq) or 'GCA_' (GenBank).
%
% Name-Value Arguments
% --------------------
% outputDir : char
%     Directory where downloaded files are saved (default: current working
%     directory).
% verbose : logical
%     Print download progress to the command window (default: true).
%
% Returns
% -------
% gffFile : char
%     Full path to the downloaded GFF3 annotation file.
% faaFile : char
%     Full path to the downloaded protein FASTA file.
%
% Examples
% --------
%     % Download to current directory
%     [gff, faa] = downloadGenomeData('GCF_000002595.2');
%
%     % Save to a specific directory
%     [gff, faa] = downloadGenomeData('GCF_000002595.2', 'data/');
%
%     % Suppress progress messages
%     [gff, faa] = downloadGenomeData('GCF_000002595.2', 'data/', false);
%
% Notes
% -----
% An internet connection is required. When both files are already present
% in outputDir, no download takes place; deleting them forces a new
% download. An error is raised when NCBI Datasets does not return both a
% GFF3 file and a protein FASTA file for the accession.

    if nargin < 1 || isempty(accession)
    error('downloadGenomeData:missingInput', ...
        ['accession is required. ', ...
         'Provide an NCBI assembly accession, e.g. ''GCF_000002595.2''.']);
    end
    accession = char(accession);
    
    if ~any(strncmpi(accession, {'GCF_', 'GCA_'}, 4))
        error('downloadGenomeData:invalidAccession', ...
        ['Accession must start with ''GCF_'' (RefSeq) or ''GCA_'' (GenBank). ', ...
         'Received: ''%s''.'], accession);
    end

    p = parseRAVENargs(varargin, {'outputDir', pwd; 'verbose', true});
    downloadDir = char(p.outputDir);
    verbose = p.verbose;

    if ~exist(downloadDir, 'dir')
        mkdir(downloadDir);
    end

    gffFile = fullfile(downloadDir, [accession '_genomic.gff']);
    faaFile = fullfile(downloadDir, [accession '_protein.faa']);

    % Skip the download if both files are already present
    if isfile(gffFile) && isfile(faaFile)
        if verbose; fprintf('Genome data for %s already present, skipping download.\n', accession); end
        return
    end

    % --- NCBI Datasets v2 API: download GFF3 and protein FASTA as one zip ---
    if verbose; fprintf('Downloading genome data for %s from NCBI Datasets...\n', accession); end
    url = ['https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/', ...
           accession, '/download?include_annotation_type=GENOME_GFF,PROT_FASTA'];
    zipFile    = [tempname '.zip'];
    extractDir = tempname;
    try
        websave(zipFile, url);
    catch ME
        error('downloadGenomeData:download', 'Download failed for %s: %s', accession, ME.message);
    end
    try
        unzip(zipFile, extractDir);
    catch ME
        delete(zipFile);
        error('downloadGenomeData:extract', 'Extraction failed for %s: %s', accession, ME.message);
    end
    delete(zipFile);

    % The archive stores the files under ncbi_dataset/data/<assembly>/
    dataParent = fullfile(extractDir, 'ncbi_dataset', 'data');
    sub = dir(dataParent);
    sub = sub([sub.isdir] & ~startsWith({sub.name}, '.'));
    if ~isempty(sub)
        assemblyDir = fullfile(dataParent, sub(1).name);
        srcGff = fullfile(assemblyDir, 'genomic.gff');
        srcFaa = fullfile(assemblyDir, 'protein.faa');
    else
        srcGff = '';
        srcFaa = '';
    end
    if ~isfile(srcGff) || ~isfile(srcFaa)
        rmdir(extractDir, 's');
        error('downloadGenomeData:noAnnotation', ...
            ['NCBI Datasets did not return both a GFF3 annotation and a protein ', ...
             'FASTA for %s. The assembly may lack annotation.'], accession);
    end
    movefile(srcGff, gffFile, 'f');
    movefile(srcFaa, faaFile, 'f');
    rmdir(extractDir, 's');
end