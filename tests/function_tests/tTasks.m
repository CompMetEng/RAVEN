classdef tTasks < RavenTestCase
% tTasks  Tests for the metabolic-task functions in tasks/.

    methods (Test)

        function checkTasksSucceedsOnFullModel(testCase)
            m = testCase.taskTestModel();
            task = testCase.taskTestStruct();
            [~, report] = evalc('checkTasks(m, [], true, false, false, task);');
            testCase.verifyTrue(report.ok == 1);
        end

        function checkTasksFailsWithoutKeyReaction(testCase)
            m = removeReactions(testCase.taskTestModel(), {'R2'});
            task = testCase.taskTestStruct();
            [~, report] = evalc('checkTasks(m, [], true, false, false, task);');
            testCase.verifyTrue(report.ok == 0);
        end

        function parseTaskListReturnsStruct(testCase)
            taskFile = fullfile(testCase.ravenRoot,'tests','function_tests', ...
                'test_data','test_tasks.txt');
            testCase.assumeDependency(exist(taskFile,'file')==2, 'test_tasks.txt');
            ts = parseTaskList(taskFile);
            testCase.verifyClass(ts, 'struct');
            testCase.verifyNotEmpty(ts);
        end

        function checkTasksWarnsWhenSeveralMetabolitesShareAName(testCase)
            m = testCase.taskTestModel();
            m.metNames{3} = 'a';   % bc is now also called a[s], like as
            m.metComps(3) = 1;
            task = testCase.taskTestStruct();
            lastwarn('');
            [~, report] = evalc('checkTasks(m, [], true, false, false, task);');
            testCase.verifySubstring(lastwarn, 'Several model metabolites match');
            testCase.verifyTrue(report.ok == 1);
        end

        function parseTaskListDefaultsMissingOptionalColumns(testCase)
            % Absent optional columns take their documented defaults: bounds
            % 0 (lower) and 1000 (upper), and empty comments.
            f = tTasks.writeTaskFile({
                {'' 'ID' 'IN' 'OUT'}
                {'' 'T1' 'A[c]' 'B[c]'}
                });
            ts = parseTaskList(f);
            testCase.verifyEqual(ts.LBin, 0);
            testCase.verifyEqual(ts.UBin, 1000);
            testCase.verifyEqual(ts.LBout, 0);
            testCase.verifyEqual(ts.UBout, 1000);
            testCase.verifyEqual(ts.comments, '');
        end

        function parseTaskListSkipsLeadingCommentRows(testCase)
            % Rows starting with a non-empty cell are comments, also when
            % they come before the header row.
            f = tTasks.writeTaskFile({
                {'# first comment'}
                {'# second comment'}
                {'' 'ID' 'IN'}
                {'' 'T1' 'A[c]'}
                });
            ts = parseTaskList(f);
            testCase.verifyEqual(numel(ts), 1);
            testCase.verifyEqual(ts.id, 'T1');
            testCase.verifyEqual(ts.inputs, {'A[c]'});
        end

        function parseTaskListWarnsAndSkipsRowsBeforeFirstId(testCase)
            f = tTasks.writeTaskFile({
                {'' 'ID' 'IN'}
                {'' '' 'X[c]'}
                {'' 'T1' 'A[c]'}
                });
            lastwarn('');
            ts = parseTaskList(f);
            testCase.verifySubstring(lastwarn, 'before the first task id');
            testCase.verifyEqual(numel(ts), 1);
            testCase.verifyEqual(ts.id, 'T1');
            testCase.verifyEqual(ts.inputs, {'A[c]'});
        end

        function parseTaskListRefusesNonNumericBound(testCase)
            f = tTasks.writeTaskFile({
                {'' 'ID' 'IN' 'IN LB'}
                {'' 'T1' 'A[c]' 'abc'}
                });
            testCase.verifyError(@() parseTaskList(f), 'RAVEN:badInput');
        end

        function parseTaskListRefusesXls(testCase)
            % The message names the supported formats, also for an upper
            % case extension.
            for ext = {'.xls','.XLS'}
                f = [tempname ext{1}];
                fid = fopen(f,'w'); fclose(fid);
                cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
                try
                    parseTaskList(f);
                    msg = '';
                    id = '';
                catch err
                    msg = err.message;
                    id = err.identifier;
                end
                testCase.verifyEqual(id, 'RAVEN:badInput');
                testCase.verifySubstring(msg, '.xlsx');
                testCase.verifySubstring(msg, '.txt');
            end
        end

        function parseTaskListIgnoresEmptyRows(testCase)
            % A row without content, also one of white space, is neither a
            % header nor a task row.
            f = tTasks.writeTaskFile({
                {'# comment' '' '' ''}
                {'' '' '' ''}
                {' ' '  ' '' ''}
                {'' 'ID' 'IN' 'IN LB'}
                {'' '' '' ''}
                {'' 'T1' 'A[c]' '2'}
                });
            ts = parseTaskList(f);
            testCase.verifyEqual(numel(ts), 1);
            testCase.verifyEqual(ts.id, 'T1');
            testCase.verifyEqual(ts.LBin, 2);

            x = [tempname '.xlsx'];
            writecell({'# comment' '' '' ''; '' '' '' ''; '' 'ID' 'IN' 'IN LB'; '' 'T1' 'A[c]' 2}, x, 'Sheet', 'TASKS');
            cleanupX = onCleanup(@() delete(x)); %#ok<NASGU>
            ts = parseTaskList(x);
            testCase.verifyEqual(numel(ts), 1);
            testCase.verifyEqual(ts.id, 'T1');
            testCase.verifyEqual(ts.LBin, 2);
        end

        function parseTaskListReportsDataRowOfBadNumber(testCase)
            % The row is counted after the header, comment rows and empty
            % rows are removed, and the message says so.
            f = tTasks.writeTaskFile({
                {'# c1' '' '' ''}
                {'' '' '' ''}
                {'' 'ID' 'IN' 'IN LB'}
                {'# c3' '' '' ''}
                {'' 'T1' 'A[c]' '1'}
                {'' 'T2' 'A[c]' 'abc'}
                });
            try
                parseTaskList(f);
                msg = '';
            catch err
                msg = err.message;
            end
            testCase.verifySubstring(msg, 'task list data row 2 is not a number');
            testCase.verifySubstring(msg, 'counted after');
        end

        function parseTaskListRefusesTextThatOnlyLooksLikeANumber(testCase)
            % str2double reads '1,5' as 15 and '2i' as a complex number.
            for bad = {'1,5','1,000','2i','1 2','1_000','0x10'}
                f = tTasks.writeTaskFile({
                    {'' 'ID' 'IN' 'IN LB'}
                    {'' 'T1' 'A[c]' bad{1}}
                    });
                testCase.verifyError(@() parseTaskList(f), 'RAVEN:badInput', bad{1});
            end
            x = [tempname '.xlsx'];
            writecell({'' 'ID' 'IN' 'IN LB'; '' 'T1' 'A[c]' '1,5'; '# note' '' '' ''}, x, 'Sheet', 'TASKS');
            cleanupX = onCleanup(@() delete(x)); %#ok<NASGU>
            testCase.verifyError(@() parseTaskList(x), 'RAVEN:badInput');
            f = tTasks.writeTaskFile({
                {'' 'ID' 'IN' 'IN LB' 'IN UB'}
                {'' 'T1' 'A[c]' '-1e3' 'Inf'}
                {'' 'T2' 'A[c]' ' 2.5 ' '.5'}
                });
            ts = parseTaskList(f);
            testCase.verifyEqual(ts(1).LBin, -1000);
            testCase.verifyEqual(ts(1).UBin, Inf);
            testCase.verifyEqual(ts(2).LBin, 2.5);
            testCase.verifyEqual(ts(2).UBin, 0.5);
        end

        function parseTaskListConvertsNumbersTypedAsText(testCase)
            % A bound that Excel holds as text is read as the number, as in
            % a tab-delimited file; other text is still refused.
            f = [tempname '.xlsx'];
            % a filled first column keeps the empty header cell in the sheet range
            writecell({'' 'ID' 'IN' 'IN LB'; '' 'T1' 'A[c]' '2'; '# note' '' '' ''}, f, 'Sheet', 'TASKS');
            cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
            ts = parseTaskList(f);
            testCase.verifyEqual(ts.LBin, 2);

            g = [tempname '.xlsx'];
            writecell({'' 'ID' 'IN' 'IN LB'; '' 'T1' 'A[c]' 'abc'; '# note' '' '' ''}, g, 'Sheet', 'TASKS');
            cleanupG = onCleanup(@() delete(g)); %#ok<NASGU>
            testCase.verifyError(@() parseTaskList(g), 'RAVEN:badInput');
        end

        function parseTaskListExcelMatchesText(testCase)
            % The same task definitions give the same structs from Excel and
            % from a tab-delimited file.
            txtFile = fullfile(testCase.ravenRoot,'tests','function_tests', ...
                'test_data','test_tasks.txt');
            testCase.assumeDependency(exist(txtFile,'file')==2, 'test_tasks.txt');
            cells = readcell(txtFile, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'char');
            cells(cellfun(@(x) all(ismissing(x)), cells)) = {[]};
            x = [tempname '.xlsx'];
            writecell(cells, x, 'Sheet', 'TASKS');
            cleanup = onCleanup(@() delete(x)); %#ok<NASGU>
            testCase.verifyEqual(parseTaskList(x), parseTaskList(txtFile));
        end

        function parseTaskListReadsXlsx(testCase)
            f = [tempname '.xlsx'];
            % a filled first column keeps the empty header cell in the sheet range
            writecell({'' 'ID' 'IN' 'IN LB'; '' 'T1' 'A[c];B[c]' 2; '' 'T2' 'C[c]' ''; '# note' '' '' ''}, f, 'Sheet', 'TASKS');
            cleanup = onCleanup(@() delete(f));
            ts = parseTaskList(f);
            testCase.verifyEqual(numel(ts), 2);
            testCase.verifyEqual(ts(1).inputs, {'A[c]';'B[c]'});
            testCase.verifyEqual(ts(1).LBin, [2;2]);
            testCase.verifyEqual(ts(2).LBin, 0);
        end

        function parseTaskListReadsFlagsListsAndHeadersLeniently(testCase)
            % Flags: "false"/"0"/"no" are false, other text is true. Lists
            % separated by ";" are trimmed and empty items dropped. Headers
            % are matched ignoring case and whitespace. An empty CHANGED
            % bound is NaN (leave the bound unchanged).
            rows = {
                {'' 'id' 'Description' 'Should Fail' 'IN' ' in lb ' 'IN UB' 'OUT' 'OUT LB' 'OUT UB' 'EQU' 'EQU LB' 'EQU UB' 'changed rxn' 'CHANGED LB' 'CHANGED UB' 'PRINT FLUX' 'COMMENTS'}
                {'' 'T1' 'first' 'false' 'A[c]; B[c];' '' '5' '' '' '' '' '' '' 'R1; R2' '' '10' '0' ''}
                {'' 'T2' 'second' 'yes' ';C[c]' '2' '' '' '' '' '' '' '' 'R3' '-5' '' 'No' ''}
                {'' 'T3' 'third' '0' '' '' '' 'D[c]' '' '' '' '' '' '' '' '' 'TRUE' ''}
                };
            nCols = 20;
            lines = cell(numel(rows),1);
            for k = 1:numel(rows)
                r = rows{k};
                r(end+1:nCols) = {''};
                lines{k} = strjoin(r, sprintf('\t'));
            end
            f = [tempname '.txt'];
            fid = fopen(f,'w');
            fprintf(fid,'%s\n',lines{:});
            fclose(fid);
            cleanup = onCleanup(@() delete(f));

            ts = parseTaskList(f);
            testCase.verifyEqual(numel(ts), 3);

            testCase.verifyFalse(ts(1).shouldFail);
            testCase.verifyFalse(ts(1).printFluxes);
            testCase.verifyEqual(ts(1).inputs, {'A[c]';'B[c]'});
            testCase.verifyEqual(ts(1).LBin, [0;0]);
            testCase.verifyEqual(ts(1).UBin, [5;5]);
            testCase.verifyEqual(ts(1).changed, {'R1';'R2'});
            testCase.verifyTrue(all(isnan(ts(1).LBrxn)));
            testCase.verifyEqual(ts(1).UBrxn, [10;10]);

            testCase.verifyTrue(ts(2).shouldFail);
            testCase.verifyFalse(ts(2).printFluxes);
            testCase.verifyEqual(ts(2).inputs, {'C[c]'});
            testCase.verifyEqual(ts(2).LBin, 2);
            testCase.verifyEqual(ts(2).LBrxn, -5);
            testCase.verifyTrue(isnan(ts(2).UBrxn));

            testCase.verifyFalse(ts(3).shouldFail);
            testCase.verifyTrue(ts(3).printFluxes);
        end

        function checkTasksLeavesNaNChangedBoundsAlone(testCase)
            m = testCase.taskTestModel();
            task = testCase.taskTestStruct();
            task.changed = {'R2'};
            task.LBrxn = NaN;
            task.UBrxn = 1000;
            [~, report] = evalc('checkTasks(m, [], true, false, false, task);');
            testCase.verifyTrue(report.ok == 1);

            % R2 can no longer carry flux from a[s] to a[c] with UB 0, and
            % its lower bound (-1000) is kept
            task.UBrxn = 0;
            [~, report] = evalc('checkTasks(m, [], true, false, false, task);');
            testCase.verifyTrue(report.ok == 0);
        end

        function parseTaskListMatchesExpected(testCase)
            % The text task definitions must parse to structs with exactly the
            % expected fields and values.
            dataDir   = fullfile(testCase.ravenRoot,'tests','function_tests','test_data');
            txtFile   = fullfile(dataDir,'test_tasks.txt');
            testCase.assumeDependency(exist(txtFile,'file')==2, 'test_tasks.txt');

            taskStruct      = parseTaskList(txtFile);

            expTask1.id          = 'ER';
            expTask1.description = 'Aerobic rephosphorylation of ATP from glucose';
            expTask1.shouldFail  = true;
            expTask1.printFluxes = true;
            expTask1.comments    = 'Messed up reaction';
            expTask1.inputs      = {'O2[s]';'glucose[s]'};
            expTask1.LBin        = [23.6;23.6];
            expTask1.UBin        = [23.8;23.8];
            expTask1.outputs     = {'H2O[s]';'CO2[s]'};
            expTask1.LBout       = [26.1;26.1];
            expTask1.UBout       = [26.2;26.2];
            expTask1.equations   = {'ATP[c] + H2O[c] => ADP[c] + Pi[c] + H+[c]'};
            expTask1.LBequ       = 30.2;
            expTask1.UBequ       = 30.6;
            expTask1.changed     = {'ATP[a] + H2O[a] => ADP[a] + Pi[a] + H+[a]'};
            expTask1.LBrxn       = 56.2;
            % incomplete struct must not match
            testCase.verifyNotEqual(taskStruct(1), expTask1);
            expTask1.UBrxn = 60;
            testCase.verifyEqual(taskStruct(1), expTask1);

            testCase.verifyEqual(length(taskStruct), 2);

            expTask2.id          = 'BS';
            expTask2.description = 'ATP de novo synthesis';
            expTask2.shouldFail  = false;
            expTask2.printFluxes = false;
            expTask2.comments    = '';
            expTask2.inputs      = {'O2[s]';'glucose[s]';'NH3[s]';'Pi[s]'};
            expTask2.LBin        = [0;0;0;0];
            expTask2.UBin        = [1000;1000;1000;1000];
            expTask2.outputs     = {'H2O[s]';'CO2[s]';'ATP[c]'};
            expTask2.LBout       = [0;0;1];
            expTask2.UBout       = [1.1;1.1;1.3];
            expTask2.equations   = {'ATP[c] + H2O[c] <=> ADP[c] + Pi[c] + H+[c]'};
            expTask2.LBequ       = -1000;
            expTask2.UBequ       = 1000;
            expTask2.changed     = {};
            expTask2.LBrxn       = [];
            expTask2.UBrxn       = [];
            testCase.verifyEqual(taskStruct(2), expTask2);
        end

    end

    methods (Static, Access = private)
        function f = writeTaskFile(rows)
            % Write rows (cell arrays of cells) as a tab-delimited task list.
            % The file is deleted when the calling test ends.
            lines = cell(numel(rows),1);
            for k = 1:numel(rows)
                lines{k} = strjoin(rows{k}, sprintf('\t'));
            end
            f = [tempname '.txt'];
            fid = fopen(f,'w');
            fprintf(fid,'%s\n',lines{:});
            fclose(fid);
            evalin('caller', sprintf('cleanup = onCleanup(@() delete(''%s''));', f));
        end
    end
end
