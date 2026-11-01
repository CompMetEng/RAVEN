function blastStructure=makeFakeBlastStructure(orthologList,sourceModelID,getModelFor)
% makeFakeBlastStructure  Make a fake blastStructure from an ortholog list.
%
% The structure has the format that getBlast generates. It feeds a
% predefined list of orthologs to getModelFromHomology. The orthologList
% is a cell array with at least 10 rows and two columns: the first column
% contains the genes from the source organism and the second column the
% genes from the target organism.
%
% Parameters
% ----------
% orthologList : cell
%     cell array of orthologous genes, where the first column contains the
%     genes from the source organism, while the second column contains the
%     genes from the target organism.
% sourceModelID : char
%     ID of the model that will be used as template, that contains the
%     genes in the first column of orthologList.
% getModelFor : char
%     the name of the organism to build a model for, identical to the
%     getModelFor parameter in the getModelFromHomology function.
%
% Returns
% -------
% blastStructure : struct
%     a fake blastStructure with one element per direction. For every pair,
%     evalue is 0, identity is 100, aligLen is 1000, bitscore is 100 and
%     ppos is 100. All orthologous pairs pass the filters of
%     getModelFromHomology.
%
% Examples
% --------
%     blastStructure = makeFakeBlastStructure(orthologList,sourceModelID,getModelFor);
%
% See also
% --------
% getModelFromHomology, getBlast, hitsFromSuppOrth

if nargin<3
    error('All three parameters should be set');
end
sourceModelID=char(sourceModelID);
getModelFor=char(getModelFor);

% Check whether orthologList is a cell array, with two columns, and more
% than 10 rows
if (~iscell(orthologList) || ~(size(orthologList,2)==2) || ~(size(orthologList,1)>=10))
    error('orthologList should be a cell array with two columns');
end

blastStructure(1).fromId=sourceModelID;
blastStructure(1).toId=getModelFor;
blastStructure(1).fromGenes=orthologList(:,1);
blastStructure(1).toGenes=orthologList(:,2);

blastStructure(2).fromId=getModelFor;
blastStructure(2).toId=sourceModelID;
blastStructure(2).fromGenes=orthologList(:,2);
blastStructure(2).toGenes=orthologList(:,1);

blastStructure(1).evalue=zeros(size(orthologList,1),1);
blastStructure(2).evalue=zeros(size(orthologList,1),1);
blastStructure(1).identity(1:size(orthologList,1),1)=100;
blastStructure(2).identity(1:size(orthologList,1),1)=100;
blastStructure(1).aligLen(1:size(orthologList,1),1)=1000;
blastStructure(2).aligLen(1:size(orthologList,1),1)=1000;
blastStructure(1).bitscore(1:size(orthologList,1),1)=100;
blastStructure(2).bitscore(1:size(orthologList,1),1)=100;
blastStructure(1).ppos(1:size(orthologList,1),1)=100;
blastStructure(2).ppos(1:size(orthologList,1),1)=100;
end
