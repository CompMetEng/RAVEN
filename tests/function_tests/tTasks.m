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
            % The text and Excel task definitions must parse to identical
            % structs with exactly the expected fields and values.
            dataDir   = fullfile(testCase.ravenRoot,'tests','function_tests','test_data');
            txtFile   = fullfile(dataDir,'test_tasks.txt');
            xlsFile   = fullfile(dataDir,'test_tasks.xls');
            testCase.assumeDependency(exist(txtFile,'file')==2, 'test_tasks.txt');
            testCase.assumeDependency(exist(xlsFile,'file')==2, 'test_tasks.xls');

            taskStruct      = parseTaskList(txtFile);
            taskStructExcel = parseTaskList(xlsFile);

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
            testCase.verifyEqual(taskStructExcel(1), expTask1);

            testCase.verifyEqual(length(taskStruct), 2);
            testCase.verifyEqual(length(taskStructExcel), 2);

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
            testCase.verifyEqual(taskStructExcel(2), expTask2);

            % the two formats must produce identical results
            testCase.verifyEqual(taskStruct, taskStructExcel);
        end

    end
end
