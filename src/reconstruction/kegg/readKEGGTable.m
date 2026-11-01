function tbl=readKEGGTable(gzFile)
% readKEGGTable  Read a gzipped KEGG relational table into a table of cellstr columns.
%
% Every column is a cellstr (RAVEN's cell-array-of-char convention), whatever the
% default text-import type of the caller.
%
% Parameters
% ----------
% gzFile : char
%     path to the gzipped TSV file, for example '.../kegg118_ko_reaction.tsv.gz'.
%
% Returns
% -------
% tbl : table
%     the table, with every column as a cellstr.

plainFile=gzFile(1:end-3); %strip the trailing '.gz'
if ~isfile(plainFile)
    gunzip(gzFile);
end
opts=detectImportOptions(plainFile,'FileType','text','Delimiter',char(9));
opts=setvartype(opts,opts.VariableNames,'char');
tbl=readtable(plainFile,opts);
end
