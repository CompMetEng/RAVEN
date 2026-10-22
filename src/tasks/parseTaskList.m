function taskStruct=parseTaskList(inputFile)
% parseTaskList  Parse a task list file.
%
% Parameters
% ----------
% inputFile : char
%     a task list in either Excel (*.xlsx, with a sheet named TASKS with
%     all relevant content) or tab-delimited (*.txt) format. The old binary
%     *.xls format is not supported. The file may contain the following
%     column headers (note, all rows starting with a non-empty cell are
%     removed; the first row after that is considered the headers; the
%     header row itself must therefore start with an empty cell):
%
%     - ID : the only required header. Each task must have a unique id
%       (string or numeric). Tasks can span multiple rows, only the first
%       row in each task should have an id
%     - DESCRIPTION : description of the task
%     - IN : allowed input(s) for the task. Metabolite names should be on
%       the form "model.metName[model.comps]". Several inputs can be
%       delimited by ";". If so, then the same bounds are used for all
%       inputs. If that is not wanted, then use several rows for the task
%     - IN LB : lower bound for the uptake of the metabolites in the row
%       (optional, default 0 which corresponds to a minimal uptake of 0
%       units)
%     - IN UB : upper bound for the uptake of the metabolites in the row
%       (optional, default 1000 which corresponds to a maximal uptake of
%       1000 units)
%     - OUT : allowed output(s) for the task (see IN)
%     - OUT LB : lower bound for the production of the metabolites in the
%       row (optional, default 0 which corresponds to a minimal production
%       of 0 units)
%     - OUT UB : upper bound for the production of the metabolites in the
%       row (optional, default 1000 which corresponds to a maximal
%       production of 1000 units)
%     - EQU : equation to add. The equation should be on the form
%       "0.4 A + 2 B <=> (or =>) C" and the metabolites should be on the
%       form "model.metName[model.comps]" (optional)
%     - EQU LB : lower bound for the equation (optional, default -1000 for
%       reversible and 0 for irreversible)
%     - EQU UB : upper bound for the equation (optional, default 1000)
%     - CHANGED RXN : reaction ID for which to change the bounds for.
%       Several IDs can be delimited by ";". If so, then the same bounds
%       are used for all reactions. If that is not wanted, then use several
%       rows for the task
%     - CHANGED LB : lower bound for the reaction. If empty, the lower bound
%       of the reaction is left as it is
%     - CHANGED UB : upper bound for the reaction. If empty, the upper bound
%       of the reaction is left as it is
%     - SHOULD FAIL : true if the correct behavior of the model is to not
%       have a feasible solution given the constraints. Any non-empty value
%       except 0, false and no (any case) counts as true (optional, default
%       false)
%     - PRINT FLUX : true if the function should print the corresponding
%       flux distribution for a task. Can be useful for testing. Read in
%       the same way as SHOULD FAIL (optional, default false)
%
%     Headers are matched ignoring case and surrounding whitespace. In the
%     ";"-delimited cells, items are trimmed and empty items are dropped.
%
% Returns
% -------
% taskStruct : struct
%     array of structures with the following fields:
%
%     - id : the id of the task
%     - description : the description of the task
%     - shouldFail : true if the task should fail
%     - printFluxes : true if the fluxes should be printed
%     - comments : string with comments
%     - inputs : cell array with input metabolites (in the form metName[comps])
%     - LBin : array with lower bounds on inputs (default, 0)
%     - UBin : array with upper bounds on inputs (default, 1000)
%     - outputs : cell array with output metabolites (in the form metName[comps])
%     - LBout : array with lower bounds on outputs (default, 0)
%     - UBout : array with upper bounds on outputs (default, 1000)
%     - equations : cell array with equations (with mets in the form metName[comps])
%     - LBequ : array with lower bounds on equations (default, -1000 for
%       reversible and 0 for irreversible)
%     - UBequ : array with upper bounds on equations (default, 1000)
%     - changed : cell array with reactions to change bounds for
%     - LBrxn : array with lower bounds on changed reactions (NaN if empty,
%       meaning the bound is left unchanged)
%     - UBrxn : array with upper bounds on changed reactions (NaN if empty,
%       meaning the bound is left unchanged)
%
% Examples
% --------
%     taskStruct = parseTaskList(inputFile);
%
% Notes
% -----
% This function is used for defining a set of tasks for a model to perform.
% The tasks are defined by defining constraints on the model, and if the
% problem is feasible, then the task is considered successful. In general,
% each row can contain one constraint on uptakes, one constraint on
% outputs, one new equation, and one change of reaction bounds. If more
% bounds are needed to define the task, then several rows can be used for
% each task. To perform the task use checkTasks or fitTasks.
%
% The general metabolites "ALLMETS" and "ALLMETSIN[comps]" can be used as
% inputs or outputs in the similar manner to normal metabolites. This is a
% convenient way to, for example, allow excretion of all metabolites to
% check whether it is the synthesis of some metabolite that is limiting or
% whether it is the degradation of some byproduct. One important difference
% is that only the upper bounds are used for these general metabolites.
% That is, you can only say that uptake or excretion is allowed, not that
% it is required. This is to avoid conflicts where the constraints for the
% general metabolites overwrite those of the real ones.
%
% See also
% --------
% checkTasks, fitTasks

if ~isfile(inputFile)
    error('Task list %s cannot be found',string(inputFile));
end

%Load the tasks file
if strcmp(extractAfter(inputFile,strlength(inputFile) - 4), '.txt')
    %load from tab delimited text file
    fid = fopen(inputFile);
    %Need to read numeric columns as strings, this is converted further
    %down. If not, the titles would be lost.
    C = textscan(fid,'%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%*[^\n]', 'Delimiter', '\t');
    fclose(fid);
    raw = [C{:}];%unnest the cell array of cell arrays into a 2-dim cell array
else
    if endsWith(lower(char(inputFile)),'.xls')
        error('RAVEN:badInput', '%s', [inputFile ': the old binary .xls format is not supported. Save the task list as .xlsx (sheet ''TASKS'') or as a tab-delimited .txt file.']);
    end
    %Read the TASKS sheet from an Excel file
    try
        raw=readcell(inputFile,'Sheet','TASKS');
    catch
        error('RAVEN:badInput', '%s', ['Could not load sheet "TASKS" from ' inputFile]);
    end
    %Normalise blank cells (which readcell returns as "missing") to [], as
    %expected by the downstream cleanSheet/parsing.
    raw(cellfun(@(x) all(ismissing(x)),raw))={[]};
end

%Rows without any content carry no information (a header must not be
%mistaken for one of them). fileRow keeps the row numbers of the file for
%error messages.
if iscell(raw)
    emptyRow=all(cellfun(@(x) isempty(x) || (ischar(x) && isempty(strtrim(x))),raw),2);
    fileRow=find(~emptyRow);
    raw=raw(~emptyRow,:);
end

%Leading rows that start with a non-empty cell are comments; the first row
%with an empty first cell is the header row
if iscell(raw)
    if isempty(raw)
        error('RAVEN:badInput', '%s', 'The TASKS sheet must have a column named ID');
    end
    firstBlank=cellfun(@(x) isempty(x) || (ischar(x) && isempty(strtrim(x))),raw(:,1));
    first=find(firstBlank,1);
    if isempty(first)
        error('RAVEN:badInput', '%s', 'The TASKS sheet must have a column named ID');
    end
    raw=raw(first:end,:);
    fileRow=fileRow(first:end);
end

%Remove all lines starting with "#" (or actually any character) and all
%empty columns
[raw,keptRows]=cleanSheet(raw);
fileRow=fileRow(keptRows);

%Captions
columns={'ID';'DESCRIPTION';'IN';'IN LB';'IN UB';'OUT';'OUT LB';'OUT UB';'EQU';'EQU LB';'EQU UB';'CHANGED RXN';'CHANGED LB';'CHANGED UB';'SHOULD FAIL';'PRINT FLUX';'COMMENTS'};
%Match the columns. Headers are matched ignoring case and surrounding
%whitespace, so "id" and " IN LB " are accepted.
header=raw(1,:);
isText=cellfun(@ischar,header);
header(isText)=upper(strtrim(header(isText)));
[I, colI]=ismember(columns,header);

%Only ID is required. Columns that are absent get index 0 from ismember, which
%cannot be used to index raw, so point each of them at its own spare empty
%column instead: every optional field then reads as empty and its default
%applies.
absent=find(colI==0);
nOld=size(raw,2);
raw(:,nOld+(1:numel(absent)))={[]};
colI(absent)=nOld+(1:numel(absent));

%Text cells (every cell of a text file, and numbers typed as text in Excel)
%are converted to numbers; text that is not a number stays and is refused
%below
numericColumns = [0 0 0 1 1 0 1 1 0 1 1 0 1 1 0 0 0] == 1;
cols = colI(numericColumns);
numeric = cellfun(@isnumeric,raw(:,cols));
%trick to avoid messing up the title row:
numeric(1,:) = 1;
for colind = 1:numel(cols)
    col = cols(colind);
    for r = find(~numeric(:,colind))'
        x = raw{r,col};
        if ischar(x) && isempty(strtrim(x))
            raw{r,col} = [];
        elseif ischar(x)
            value = str2double(x);
            if ~isnan(value) || strcmpi(strtrim(x),'nan')
                raw{r,col} = value;
            end
        end
    end
end

%Every bound that is given must be a number
numericNames = columns([4 5 7 8 10 11 13 14]);
numericCols = colI([4 5 7 8 10 11 13 14]);
for k = 1:numel(numericCols)
    for r = 2:size(raw,1)
        x = raw{r,numericCols(k)};
        if ~isempty(x) && ~(isnumeric(x) && isscalar(x))
            error('RAVEN:badInput', 'The value in column %s of row %d of the task list is not a number', numericNames{k}, fileRow(r));
        end
    end
end

%Check that the ID field is present
if I(1)==0
    EM='The TASKS sheet must have a column named ID';
    error('RAVEN:badInput', '%s', EM);
end

%Rows before the first task id belong to no task
firstId=find(~cellfun(@isempty,raw(2:end,colI(1))),1)+1;
if isempty(firstId)
    firstId=size(raw,1)+1;
end
if firstId>2
    warning('RAVEN:warning', 'Skipping %d row(s) before the first task id', firstId-2);
    raw(2:firstId-1,:)=[];
    fileRow(2:firstId-1)=[];
end
if size(raw,1)<2
    taskStruct=[];
    return;
end

%Add default bounds where needed
for i=[4 5 7 8]
    I=cellfun(@isempty,raw(:,colI(i)));
    if i==5 || i==8
        raw(I,colI(i))={1000};
    else
        raw(I,colI(i))={0};
    end
end

%Create an empty task structure
eTask.id='';
eTask.description='';
eTask.shouldFail=false;
eTask.printFluxes=false;
eTask.comments='';
eTask.inputs={};
eTask.LBin=[];
eTask.UBin=[];
eTask.outputs={};
eTask.LBout=[];
eTask.UBout=[];
eTask.equations={};
eTask.LBequ=[];
eTask.UBequ=[];
eTask.changed={};
eTask.LBrxn=[];
eTask.UBrxn=[];

%Main loop
taskStruct=[];
task=eTask;
if isnumeric(raw{2,colI(1)})
    task.id=num2str(raw{2,colI(1)});
else
    task.id=raw{2,colI(1)};
end
task.description=raw{2,colI(2)};
if isTruthy(raw{2,colI(15)})
    task.shouldFail=true;
end
if isTruthy(raw{2,colI(16)})
    task.printFluxes=true;
end
if ~isempty(raw{2,colI(17)})
    task.comments=raw{2,colI(17)};
end

for i=2:size(raw,1)
    %Set the inputs
    if ischar(raw{i,colI(3)})
        inputs=splitList(raw{i,colI(3)});
        task.inputs=[task.inputs;inputs(:)];
        task.LBin=[task.LBin;ones(numel(inputs),1)*raw{i,colI(4)}];
        task.UBin=[task.UBin;ones(numel(inputs),1)*raw{i,colI(5)}];
    end
    %Set the outputs
    if ischar(raw{i,colI(6)})
        outputs=splitList(raw{i,colI(6)});
        task.outputs=[task.outputs;outputs(:)];
        task.LBout=[task.LBout;ones(numel(outputs),1)*raw{i,colI(7)}];
        task.UBout=[task.UBout;ones(numel(outputs),1)*raw{i,colI(8)}];
    end
    %Add new rxns
    if ischar(raw{i,colI(9)})
        task.equations=[task.equations;raw{i,colI(9)}];
        if ~isempty(raw{i,colI(10)})
            task.LBequ=[task.LBequ;raw{i,colI(10)}];
        else
            if any(strfind(raw{i,colI(9)},'<=>'))
                task.LBequ=[task.LBequ;-1000];
            else
                task.LBequ=[task.LBequ;0];
            end
        end
        if ~isempty(raw{i,colI(11)})
            task.UBequ=[task.UBequ;raw{i,colI(11)}];
        else
            task.UBequ=[task.UBequ;1000];
        end
    end
    %Add changed bounds
    if ischar(raw{i,colI(12)})
        changed=splitList(raw{i,colI(12)});
        task.changed=[task.changed;changed(:)];
        %An empty bound means "leave this bound as it is", stored as NaN
        task.LBrxn=[task.LBrxn;ones(numel(changed),1)*boundOrNaN(raw{i,colI(13)})];
        task.UBrxn=[task.UBrxn;ones(numel(changed),1)*boundOrNaN(raw{i,colI(14)})];
    end
    
    %Check if it should add more constraints
    if i<size(raw,1)
        if isempty(raw{i+1,colI(1)})
            continue;
        end
    end
    
    taskStruct=[taskStruct;task];
    task=eTask;
    if i<size(raw,1)
        if isnumeric(raw{i+1,colI(1)})
            task.id=num2str(raw{i+1,colI(1)});
        else
            task.id=raw{i+1,colI(1)};
        end
        task.description=raw{i+1,colI(2)};
        if isTruthy(raw{i+1,colI(15)})
            task.shouldFail=true;
        end
        if isTruthy(raw{i+1,colI(16)})
            task.printFluxes=true;
        end
        if ~isempty(raw{i+1,colI(17)})
            task.comments=raw{i+1,colI(17)};
        end
    end
end

%Should add more checks, such as unique IDs and missing headers

end

function tf = isTruthy(x)
% A flag cell is true unless it is empty, 0, or one of the texts "0", "false"
% and "no" (any case).
if isempty(x)
    tf = false;
elseif ischar(x) || isstring(x)
    tf = ~ismember(lower(strtrim(char(x))),{'','0','false','no'});
elseif isnumeric(x) || islogical(x)
    tf = any(x(:) ~= 0);
else
    tf = true;
end
end

function items = splitList(str)
% Split a ";"-delimited cell on ";", trim each item and drop empty ones, so
% that "A; B;" gives {'A';'B'}.
items = strtrim(regexp(str,';','split'));
items = items(~cellfun(@isempty,items));
end

function b = boundOrNaN(x)
% Bound for a changed reaction: NaN (leave unchanged) when the cell is empty.
if isempty(x)
    b = NaN;
else
    b = x;
end
end
