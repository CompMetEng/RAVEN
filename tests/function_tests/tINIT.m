classdef tINIT < RavenTestCase
% tINIT  Tests for the ftINIT context-specific modelling functions in INIT/.
%
%   The MILP-based pipeline is exercised end-to-end on small synthetic models
%   (built by the local getTstModel* helpers, illustrated in tINITtestInfo/):
%   ftINITPipelineRuns / ftINITWithTaskRuns / ftINITMetabolomicsRuns /
%   ftINITFullVsThreeStepRuns cover prepINITModel + ftINIT and the internal
%   ftINITInternalAlg / groupRxnScores helpers, fitTasksPreMerged* covers the
%   gap-filling step (fitTasks with gapFillMode 'preMerged', ftINIT's mode),
%   all guarded on a MILP solver (Gurobi/SCIP). The self-contained helpers
%   (mergeLinear, groupRxnScores, reverseRxns, rescaleModelForINIT,
%   scoreModel, getExprForRxnScore, parseHPA, parseHPArna) are tested
%   directly against known results.

    methods (Test)

        function getINITStepsReturnsSteps(testCase)
            steps = getINITSteps([], '1+1');
            testCase.verifyNotEmpty(steps);
        end

        function prepINITModelRefusesTaskIdsThatAreNotText(testCase)
            % Non-text ids and a missing id field are reported as RAVEN:badInput.
            bad = {struct('id',{"a","a"}), struct('id',{1,1}), struct('name',{'a','b'})};
            for i = 1:numel(bad)
                err = [];
                try
                    prepINITModel(struct(), bad{i});
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier, 'RAVEN:badInput');
            end
        end

        function prepINITModelAcceptsRepeatedAndEmptyTaskIds(testCase)
            % Ids may repeat and may be empty: a task is one row group of the
            % task list. A task that cannot be performed is removed by its
            % position, so the task with the same id that can be performed stays.
            testCase.assumeMILPSolver();
            testModel = getTstModel();
            base = getTstModelTasks();
            second = base;  second.id = 'ER';  second.description = 'second';
            first = base;   first.id = 'BS';   first.description = 'first';
            failing = base; failing.id = 'BS'; failing.description = 'cannot be performed';
            failing.inputs = {}; failing.LBin = []; failing.UBin = [];
            third = base;   third.id = '';     third.description = 'empty id';
            fourth = base;  fourth.id = '';    fourth.description = 'empty id again';
            tasks = [first, second, failing, third, fourth];
            evalc('prepData = prepINITModel(testModel, tasks, {}, false, {}, ''s'');');
            testCase.verifyEqual({prepData.taskStruct.description}, ...
                {'first','second','empty id','empty id again'});
            testCase.verifyEqual({prepData.taskStruct.id}, {'BS','ER','',''});
        end

        function prepINITModelKeepsTasksOfARepeatedIdListLikeHumanGEM(testCase)
            % A task list in Human-GEM's layout: category codes as ids, several
            % tasks per code, comment rows above the header.
            testCase.assumeMILPSolver();
            taskFile = fullfile(fileparts(mfilename('fullpath')), ...
                'test_data','test_tasks_repeated_ids.txt');
            tasks = parseTaskList(taskFile);
            testCase.verifyEqual({tasks.id}, {'ER','BS','ER','BS'});
            evalc('prepData = prepINITModel(getTstModel(), tasks, {}, false, {}, ''s'');');
            testCase.verifyEqual({prepData.taskStruct.description}, ...
                {'Produce e from a','Produce e from a, second','Produce e from a, third'});
            testCase.verifyEqual({prepData.taskStruct.id}, {'ER','BS','BS'});
        end

        function getExprForRxnScoreRuns(testCase)
            expr = getExprForRxnScore(rand(10, 1), 1);
            testCase.verifyNotEmpty(expr);
        end

        function removeLowScoreGenesRuns(testCase)
            geneScores = randn(numel(testCase.model.genes), 1);
            evalc('m2 = removeLowScoreGenes(testCase.model, geneScores);');
            testCase.verifyClass(m2, 'struct');
        end

        function removeLowScoreGenesKeepsFieldsAligned(testCase)
            % getGenesFromGrRules returns a sorted gene list, but the
            % annotation fields are trimmed with a mask in the original gene
            % order. If newModel.genes took the sorted order, every annotation
            % field would shift relative to it. Use deliberately unsorted genes
            % so sorted != original, and drop one isozyme gene.
            m = struct();
            m.id = 'test';
            m.rxns = {'R1'}; m.rxnNames = {'R1'};
            m.mets = {'A';'B'}; m.metNames = {'A';'B'}; m.metComps = [1;1];
            m.comps = {'c'}; m.compNames = {'c'};
            m.S = sparse([-1;1]); m.lb = 0; m.ub = 1000; m.rev = 0; m.c = 0;
            m.b = [0;0];
            m.genes = {'G3';'G1';'G2'};                      % unsorted
            m.geneShortNames = {'short3';'short1';'short2'}; % aligned to genes
            m.grRules = {'G3 or G1 or G2'};
            m.rxnGeneMat = sparse([1 1 1]);
            scores = [1; -1; 1];   % G1 (negative) is dropped from the isozyme
            evalc('[mm, removed] = removeLowScoreGenes(m, scores);');

            testCase.verifyEqual(removed, {'G1'});
            % Each surviving gene keeps its own short name.
            for k = 1:numel(mm.genes)
                want = strrep(mm.genes{k}, 'G', 'short');
                testCase.verifyEqual(mm.geneShortNames{k}, want);
            end
            % rxnGeneMat columns must correspond to mm.genes, and R1 still uses
            % both surviving genes.
            testCase.verifyEqual(numel(mm.genes), size(mm.rxnGeneMat, 2));
            testCase.verifyEqual(full(mm.rxnGeneMat), ones(1, numel(mm.genes)));
        end

        function removeLowScoreGenesTieKeepsFirstGeneInRule(testCase)
            % With every isozyme negative and tied, the gene listed first in
            % the rule is kept (not the alphabetically first one), and the
            % order of the surviving genes follows the rule.
            m = struct();
            m.id = 'test';
            m.rxns = {'R1';'R2'}; m.rxnNames = {'R1';'R2'};
            m.mets = {'A';'B'}; m.metNames = {'A';'B'}; m.metComps = [1;1];
            m.comps = {'c'}; m.compNames = {'c'};
            m.S = sparse([-1 -1;1 1]); m.lb = [0;0]; m.ub = [1000;1000];
            m.rev = [0;0]; m.c = [0;0]; m.b = [0;0];
            m.genes = {'G1';'G2';'G3';'G4'};
            m.grRules = {'G2 or G1';'G4 or G3 or G2'};
            m.rxnGeneMat = sparse([1 1 0 0;0 1 1 1]);
            evalc('[mm, removed] = removeLowScoreGenes(m, [-1;-1;1;-1]);');
            testCase.verifyEqual(mm.grRules{1}, 'G2');
            testCase.verifyEqual(mm.grRules{2}, 'G3');
            testCase.verifyEqual(sort(removed), {'G1';'G4'});

            evalc('mm = removeLowScoreGenes(m, [1;1;1;1]);');
            testCase.verifyEqual(mm.grRules, m.grRules);
        end

        function reverseRxnsRuns(testCase)
            % R1 = '=> a[s]'; R3 = 'a[c] <=> b[c] + c[c]'
            testModel = getTstModel();
            tmpModel  = reverseRxns(testModel, {'R1';'R3'});
            res       = constructEquations(tmpModel, {'R1';'R3'});
            expRes    = {'a[s] => ';'b[c] + c[c] <=> a[c]'};
            testCase.verifyTrue(all(strcmp(res, expRes)));
        end

        function rescaleModelForINITRuns(testCase)
            miniModel      = struct();
            miniModel.S    = [1,1000;-1,-40];
            miniModel.rxns = {'1';'2'};
            miniModel.mets = {'1';'2'};
            res = rescaleModelForINIT(miniModel, 10);
            testCase.verifyTrue(abs(res.S(1,2) - res.S(2,2)*-10) < 10^-6);
            testCase.verifyTrue(abs((abs(res.S(1,2)) + abs(res.S(2,2)))/2) - 1 < 10^-6);
        end

        function mergeLinearGroupsAndScores(testCase)
            % mergeLinear merges linear reaction chains; groupRxnScores then
            % aggregates the per-reaction scores onto the merged groups.
            testModel     = getTstModel();
            testRxnScores = getTstModelRxnScores();
            evalc('[reducedModel,origRxnIds,groupIds,reversedRxns] = mergeLinear(testModel, {});');
            % {R1,R2}, {R3,R5}, {R4,R6}, {R7,R8}, {R9,R10} are merged
            testCase.verifyTrue(all(groupIds == [1;1;2;3;2;3;4;4;5;5]));
            % R1, R3, R4, R7 irreversible, R9 reversible
            testCase.verifyTrue(all(reducedModel.rev == [0;0;0;0;1]));
            testCase.verifyTrue(all(reducedModel.lb == [0;0;0;0;-1000]));

            newRxnScores = groupRxnScores(reducedModel, testRxnScores, ...
                origRxnIds, groupIds, ismember(origRxnIds, {'R1';'R2';'R8'}));
            testCase.verifyTrue(all(newRxnScores == [0;-0.5;7.5;-1;0.5]));

            % testModel4 has reactions that are not all linearly merged
            testModel4 = getTstModel4();
            evalc('[reducedModel,origRxnIds,groupIds,reversedRxns] = mergeLinear(testModel4, {});');
            % {R5,R6}, {R7,R8}, {R9,R10} are merged
            testCase.verifyTrue(all(groupIds == [0;0;0;0;1;1;2;2;3;3;0]));
            testCase.verifyTrue(all(reducedModel.rev == [1;0;1;0;1;1;0;0]));
            % some reactions flip direction when turned irreversible
            testCase.verifyTrue(strcmp(constructEquations(reducedModel, 'R9'), 'g[s] => '));
            testCase.verifyTrue(all(find(reversedRxns) == [6;9]));
        end

        function mergeLinearNumbersGroupsByFirstRxn(testCase)
            % Group numbers follow the first reaction of each group in
            % model.rxns, not the merge order. The metabolites are visited
            % in model.mets order, so c (R5, R6) merges first.
            model = struct();
            model.id = 'numbering';
            model.rxns = {'R1';'R2';'R3';'R4';'R5';'R6'};
            model.rxnNames = model.rxns;
            model.mets = {'z';'c';'w';'x';'a';'m';'b';'y'};
            model.metNames = model.mets;
            model.comps = {'c'};
            model.compNames = {'c'};
            model.metComps = ones(8,1);
            subs = {'x';'a';'m';'b';'z';'c'};
            prods = {'a';'m';'b';'y';'c';'w'};
            model.S = sparse(8,6);
            for i = 1:6
                model.S(strcmp(model.mets,subs{i}),i) = -1;
                model.S(strcmp(model.mets,prods{i}),i) = 1;
            end
            model.lb = zeros(6,1);
            model.ub = 1000*ones(6,1);
            model.rev = zeros(6,1);
            model.c = zeros(6,1);
            model.b = zeros(8,1);
            model.genes = cell(0,1);
            model.grRules = repmat({''},6,1);
            model.rxnGeneMat = sparse(6,0);
            evalc('[reducedModel,~,groupIds] = mergeLinear(model, {});');
            testCase.verifyEqual(groupIds(:), [1;1;1;1;2;2]);
            testCase.verifyEqual(sort(reducedModel.rxns), {'R1';'R5'});
        end

        function mergeLinearFlipsRxnInLaterPass(testCase)
            % R1 and R2 merge through A in the first pass and R2 is removed,
            % so R4 sits at index 3 of the reduced model. The second pass
            % merges R1 and R4 through Y and flips R4, which must be the
            % reaction marked as reversed (not R3, the third original one).
            model = struct();
            model.id = 'lateflip';
            model.rxns = {'R1';'R2';'R3';'R4';'R5';'R6'};
            model.rxnNames = model.rxns;
            model.mets = {'A';'Y';'Q';'W'};
            model.metNames = model.mets;
            model.comps = {'c'};
            model.compNames = {'c'};
            model.metComps = ones(4,1);
            S = zeros(4,6);
            S([1 2],1) = [1;1];   % R1: -> A + Y
            S([1 2],2) = [-1;1];  % R2: A -> Y
            S(3,3) = -1;          % R3: Q ->
            S([4 2],4) = [-1;1];  % R4: W <=> Y
            S(4,5) = 1;           % R5: -> W
            S(4,6) = -1;          % R6: W ->
            model.S = sparse(S);
            model.lb = [0;0;0;-1000;0;0];
            model.ub = 1000*ones(6,1);
            model.rev = double(model.lb < 0);
            model.c = zeros(6,1);
            model.b = zeros(4,1);
            model.genes = cell(0,1);
            model.grRules = repmat({''},6,1);
            model.rxnGeneMat = sparse(6,0);
            evalc('[reducedModel,origRxnIds,groupIds,reversedRxns] = mergeLinear(model, {});');
            testCase.verifyEqual(groupIds(:), [1;1;0;1;0;0]);
            testCase.verifyEqual(origRxnIds(reversedRxns), {'R4'});
            testCase.verifyEqual(sort(reducedModel.rxns), {'R1';'R3';'R5';'R6'});
        end

        function scoreModelRuns(testCase)
            arrayData.genes   = testCase.model.genes;
            arrayData.tissues = {'t1'};
            arrayData.levels  = abs(randn(numel(testCase.model.genes), 1)) + 1;
            arrayData.threshold = 1;
            evalc('rxnScores = scoreModel(testCase.model, [], arrayData, ''t1'', []);');
            testCase.verifyNumElements(rxnScores, numel(testCase.model.rxns));
        end

        function scoreModelHpaUnmeasuredCellTypeIsNotZero(testCase)
            % A gene not detected in one cell type of a tissue and not
            % measured in the other must keep its 'Not detected' score. An
            % unmeasured cell type is not a measurement of zero, and zero
            % ranks above 'Not detected' (-8), so scoring it as such would
            % also keep the gene from ever being pruned.
            m = testCase.model;
            hpaData.genes      = m.genes(1);
            hpaData.tissues    = {'liver';'liver'};
            hpaData.celltypes  = {'hepatocyte';'kupffer'};
            hpaData.levels     = {'Not detected'};
            hpaData.gene2Level = sparse(1,2);
            hpaData.gene2Level(1,1) = 1;   % not detected in hepatocyte
                                           % not measured in kupffer

            evalc(['[~,~,hpaScores] = scoreModel(m, hpaData, [], ' ...
                '''liver'', ''multipleCellScoring'', ''max'');']);
            testCase.verifyEqual(hpaScores(1), -8, 'AbsTol', 1e-9);

            % The average is likewise taken over the measurements that exist
            evalc(['[~,~,hpaScores] = scoreModel(m, hpaData, [], ' ...
                '''liver'', ''multipleCellScoring'', ''average'');']);
            testCase.verifyEqual(hpaScores(1), -8, 'AbsTol', 1e-9);
        end

        function scoreModelExactScores(testCase)
            % scoreModel must reproduce the known reaction scores for
            % the reference model prepared by prepINITModel.
            testCase.assumeMILPSolver();
            testModel = getTstModel();
            evalc('prepData = prepINITModel(testModel, {}, {}, false, {}, ''s'');');
            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;
            rxnScores = scoreModel(prepData.refModel, [], arrayData, arrayData.tissues{1}, []);
            testCase.verifyTrue(all(abs(rxnScores - getTstModelRxnScores()) < 10^-10));
        end

        function parseHPAparsesSyntheticDump(testCase)
            % parseHPA reads a tab-separated HPA dump with the columns
            % Gene / Gene name / Tissue / Cell type / Level / Reliability.
            % Synthesise a tiny one so the parser is exercised without the
            % multi-gigabyte external download.
            g = testCase.model.genes(1:3);
            f = [tempname '.tsv'];
            testCase.addTeardown(@() delete(f));
            fid = fopen(f, 'w');
            fprintf(fid, 'Gene\tGene name\tTissue\tCell type\tLevel\tReliability\n');
            fprintf(fid, '%s\t%s\tliver\thepatocyte\tHigh\tApproved\n', g{1}, g{1});
            fprintf(fid, '%s\t%s\tliver\thepatocyte\tMedium\tApproved\n', g{2}, g{2});
            fprintf(fid, '%s\t%s\tbrain\tneuron\tLow\tSupported\n', g{3}, g{3});
            fclose(fid);
            evalc('hpa = parseHPA(f);');
            testCase.verifyClass(hpa, 'struct');
            testCase.verifyEqual(sort(hpa.genes), sort(g));
            testCase.verifyNumElements(hpa.tissues, 2);   % liver+brain cell types
        end

        function parseHPArnaParsesSyntheticDump(testCase)
            % parseHPArna (version 19) reads Gene / Gene name / Tissue /
            % TPM / pTPM / NX and builds a genes-by-tissues TPM matrix.
            g = testCase.model.genes(1:3);
            f = [tempname '.tsv'];
            testCase.addTeardown(@() delete(f));
            fid = fopen(f, 'w');
            fprintf(fid, 'Gene\tGene name\tTissue\tTPM\tpTPM\tNX\n');
            fprintf(fid, '%s\t%s\tliver\t12.5\t10.1\t8.0\n', g{1}, g{1});
            fprintf(fid, '%s\t%s\tbrain\t3.2\t2.9\t1.5\n', g{2}, g{2});
            fprintf(fid, '%s\t%s\tliver\t0.0\t0.0\t0.0\n', g{3}, g{3});
            fclose(fid);
            evalc('rna = parseHPArna(f);');
            testCase.verifyClass(rna, 'struct');
            testCase.verifyNumElements(rna.genes, 3);
            testCase.verifyEqual(size(rna.levels), [3 2]);   % 3 genes x 2 tissues
        end

        function getINITModelRuns(testCase)
            % The tINIT path is legacy but supported, and until now nothing
            % exercised it: a refactor of the machinery it shares with ftINIT
            % (checkTasks, simplifyModel, solveLP, removeReactions) could
            % break it without any test noticing.
            testCase.assumeMILPSolver();
            refModel = getTstModel();
            % getINITModel wants the closed form, which this fixture is not
            % built in; its own documentation gives this as the way to add it
            refModel.unconstrained = false(numel(refModel.mets),1);

            arrayData.genes     = refModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;

            evalc(['m = getINITModel(refModel, arrayData.tissues{1}, ' ...
                '''arrayData'', arrayData, ''printReport'', false);']);
            testCase.verifyClass(m, 'struct');
            testCase.verifyNotEmpty(m.rxns);
            testCase.verifyTrue(all(ismember(m.rxns, refModel.rxns)));
            % R4 and R10 carry the highest scores, so neither should be cut
            testCase.verifyTrue(all(ismember({'R4';'R10'}, m.rxns)));
        end

        function ftINITPipelineRuns(testCase)
            % prepINITModel + ftINIT end-to-end on testModel without tasks.
            testCase.assumeMILPSolver();
            testModel  = getTstModel();
            testParams = struct();
            evalc('prepData = prepINITModel(testModel, {}, {}, false, {}, ''s'');');
            % R1 and R8 are exchange reactions without a GPR
            testCase.verifyTrue(all(strcmp( ...
                prepData.refModel.rxns(prepData.toIgnoreExch), {'R1';'R8'})));

            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;

            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,testParams,false);']);
            % R1 and R8 are added (no GPR, exchange); R2 (transport, no GPR) is
            % removed by the standard third step; R7 is not added (it has a GPR).
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R4';'R6';'R8';'R9';'R10'})));

            % make R7 and R10 spontaneous
            evalc('prepData = prepINITModel(testModel, {}, {''R7'';''R10''}, false, {}, ''s'');');
            testCase.verifyTrue(all(strcmp( ...
                prepData.refModel.rxns(prepData.toIgnoreExch | prepData.toIgnoreSpont), ...
                {'R1';'R7';'R8';'R10'})));
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,testParams,false);']);
            % the model now takes the "correct" path (incl. R2), skipping R9/R10
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R2';'R4';'R6';'R7';'R8'})));
        end

        function ftINITWithoutRemoveGenesDropsGenesOfRemovedReactions(testCase)
            % Genes that no remaining reaction uses are not in the output,
            % also when removeGenes is false.
            testCase.assumeMILPSolver();
            testModel = getTstModel();
            evalc('prepData = prepINITModel(testModel, {}, {}, false, {}, ''s'');');
            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),false,true,struct(),false);']);
            testCase.verifyLessThan(numel(resModel.genes), numel(testModel.genes));
            testCase.verifyTrue(all(full(any(resModel.rxnGeneMat,1))));
        end

        function ftINITResolveTiesAndProveAbsGapStillRunCorrectly(testCase)
            % resolveTies and proveAbsGap are threaded through ftINIT ->
            % ftINITInternalAlg (main extraction) and, for resolveTies, also
            % into fitTasks' gap-filling (raven-gecko-parity#104). This
            % fixture (same as ftINITPipelineRuns above) is too small to carry
            % a genuine tie or to expose the escalation's near-zero-objective
            % gap defect -- both are covered at genome scale in the PR
            % description -- so this only checks that turning either or both
            % on does not change the correct result or throw.
            % proveAbsGap is deliberately set to 0.05 here, not the
            % genome-scale-recommended 1.0: this fixture's reaction scores
            % (getTstModelRxnScores, as small as 0.5) make 1.0 genuinely loose
            % enough to accept a worse, different incumbent by design -- a
            % real property of a fixed absolute gap, not a bug -- so a value
            % below the model's own score granularity is what actually proves
            % the true optimum here.
            testCase.assumeMILPSolver();
            testModel  = getTstModel();
            testParams = struct();
            evalc('prepData = prepINITModel(testModel, {}, {}, false, {}, ''s'');');

            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;

            expected = {'R1';'R4';'R6';'R8';'R9';'R10'};

            evalc(['resTies = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,testParams,false,' ...
                '''resolveTies'',true);']);
            testCase.verifyTrue(all(strcmp(resTies.rxns, expected)));

            evalc(['resGap = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,testParams,false,' ...
                '''proveAbsGap'',0.05);']);
            testCase.verifyTrue(all(strcmp(resGap.rxns, expected)));

            evalc(['resBoth = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,testParams,false,' ...
                '''resolveTies'',true,''proveAbsGap'',0.05);']);
            testCase.verifyTrue(all(strcmp(resBoth.rxns, expected)));
        end

        function ftINITResolveTiesRanksMergedRxnBySmallestOriginalId(testCase)
            % a[c] reaches m[c] by path P (Rm, Rc) or path Q (Rn, Rb), each
            % path one linear merge group with equal negative scores. Model
            % order and group numbers put P first, and the merged reactions
            % are named by their first members Rm and Rn, but the smallest
            % original id of Q (Rb) sorts before that of P (Rc). resolveTies
            % keeps Q; ranking by the merged id, the model order or the group
            % number would keep P.
            testCase.assumeMILPSolver();
            testModel = getTstModelTwoPaths();
            evalc('prepData = prepINITModel(testModel, getTstModelTasks(), {}, false, {}, ''s'');');
            g = @(r) prepData.groupIds(strcmp(prepData.refModel.rxns, r));
            testCase.verifyEqual(g('Rm'), g('Rc'));
            testCase.verifyEqual(g('Rn'), g('Rb'));
            testCase.verifyLessThan(g('Rm'), g('Rn'));
            testCase.verifyTrue(all(ismember({'Rm';'Rn'}, prepData.minModel.rxns)));

            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore([-1;-1;-1;-1;5;5]);
            arrayData.threshold = 1;
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,struct(),false,' ...
                '''resolveTies'',true);']);
            testCase.verifyTrue(all(ismember({'Rn';'Rb'}, resModel.rxns)));
            testCase.verifyFalse(any(ismember({'Rm';'Rc'}, resModel.rxns)));
        end

        function ftINITInternalAlgResolveTiesGivesMinimalFluxes(testCase)
            % With resolveTies the fluxes come from the least-total-flux
            % solution for the tie-broken on/off state. The positive
            % reversible R3 (m <=> q) can be on in either direction: backward
            % (import q, feed the essential E) carries 0.4 in total, forward
            % (a via the free R1 into R3 and E) 0.8. The free loop L1/L2
            % carries nothing and the negative-score dead end R2 stays off.
            testCase.assumeMILPSolver();
            m.rxns = {'R1';'R2';'R3';'E';'L1';'L2';'EX_a';'EX_p';'EX_q'};
            m.mets = {'a';'m';'p';'q';'n';'z'};
            %         R1 R2 R3  E L1 L2 EXa EXp EXq
            m.S = sparse([-1 -1  0  0  0  0 -1   0   0;  % a
                           1  0 -1 -1 -1  1  0   0   0;  % m
                           0  0  0  1  0  0  0  -1   0;  % p
                           0  0  1  0  0  0  0   0  -1;  % q
                           0  0  0  0  1 -1  0   0   0;  % n
                           0  1  0  0  0  0  0   0   0]); % z
            m.lb = [0;0;-1000;0;0;0;-1000;-1000;-1000];
            m.ub = 1000*ones(9,1);
            m.rev = double(m.lb < 0);
            scores = [0;-1;1;0;0;0;0;0;0];
            fluxIn = 0.1*ones(9,1);
            [deleted,~,~,turnedOn,fluxes] = ftINITInternalAlg(m,scores,[],{'E'}, ...
                [],false,false,struct(),[],fluxIn,false,0.1,true);
            testCase.verifyTrue(turnedOn(3));
            testCase.verifyTrue(deleted(2));
            % E is forced at min(0.99*0.1,0.1) and needs the 0.1 of R3.
            expected = [0;0;-0.1;0.1;0;0;0;0.1;-0.1];
            testCase.verifyEqual(fluxes, expected, 'AbsTol', 1e-6);
        end

        function ftINITInternalAlgResolveTiesFixesEqualFluxDirections(testCase)
            % The essential E consumes y, which the reversible R1 (y <=> x)
            % or R2 (y <=> z) can supply when run backward from the imports
            % IN_x or IN_z. Both routes carry the same total flux; with
            % resolveTies the higher-ranked reaction carries the backward
            % flux: R2 when ranked by id, R1 when rankKeys rank it last. N is
            % a negative-score dead end.
            testCase.assumeMILPSolver();
            m.rxns = {'IN_x';'IN_z';'R1';'R2';'E';'OUT_p';'N'};
            m.mets = {'x';'y';'z';'p';'d'};
            %        INx INz R1 R2  E OUTp  N
            m.S = sparse([1   0  1  0  0   0 -1;  % x
                          0   0 -1 -1 -1   0  0;  % y
                          0   1  0  1  0   0  0;  % z
                          0   0  0  0  1  -1  0;  % p
                          0   0  0  0  0   0  1]); % d
            m.lb = [0;0;-1000;-1000;0;0;0];
            m.ub = 1000*ones(7,1);
            m.rev = double(m.lb < 0);
            scores = [0;0;0;0;0;0;-1];
            fluxIn = ones(7,1);
            [~,~,~,~,fluxes] = ftINITInternalAlg(m,scores,[],{'E'}, ...
                [],false,false,struct(),[],fluxIn,false,0.1,true);
            testCase.verifyEqual(fluxes(3:4), [0;-0.1], 'AbsTol', 1e-6);
            keys = {'IN_x';'IN_z';'R9';'R2';'E';'OUT_p';'N'};
            [~,~,~,~,fluxes] = ftINITInternalAlg(m,scores,[],{'E'}, ...
                [],false,false,struct(),[],fluxIn,false,0.1,true,keys);
            testCase.verifyEqual(fluxes(3:4), [-0.1;0], 'AbsTol', 1e-6);
        end

        function ftINITWithTaskRuns(testCase)
            % A task requiring e[s] from a[s] forces R2 and R7 to be essential.
            testCase.assumeMILPSolver();
            testModel      = getTstModel();
            testModelTasks = getTstModelTasks();
            testParams     = struct();
            evalc('prepData = prepINITModel(testModel, testModelTasks, {}, false, {}, ''s'');');
            % essential rxns are reported on the linearly merged model (R1, R7)
            testCase.verifyTrue(all(strcmp(prepData.essentialRxns, {'R1';'R7'})));

            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],getINITSteps(),true,true,testParams,false);']);
            % R2 and R7 essential => all rxns on except R3 and R5 (negative
            % score, not needed for the task)
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})));
        end

        function fitTasksPreMergedRuns(testCase)
            % Remove exchange rxns (required for gap filling) and create a gap
            % by removing R7; fitTasks with gapFillMode 'preMerged' must add
            % R7 back.
            testCase.assumeMILPSolver();
            testModel      = getTstModel();
            testModelTasks = getTstModelTasks();
            testRxnScores  = getTstModelRxnScores();

            mTempRef = closeModel(testModel);
            mTempRef = removeReactions(mTempRef, {'R1';'R8'});
            mTemp    = removeReactions(mTempRef, {'R7'});
            mTemp.id = 'tmp';
            tmpRxnScores = testRxnScores([2;3;4;5;6;7;9;10]);
            evalc(['[~,addedRxnMat] = fitTasks(mTemp,mTempRef,[],' ...
                '''printOutput'',false,''rxnScores'',min(tmpRxnScores,-0.1),' ...
                '''taskStructure'',testModelTasks,''gapFillMode'',''preMerged'',' ...
                '''params'',struct(),''verbose'',false);']);
            % An equality rather than all(strcmp(...)): an empty selection
            % makes strcmp return an empty logical, which all() reports as
            % true, so a run that added nothing would pass silently.
            testCase.verifyEqual(mTempRef.rxns(any(addedRxnMat,2)), {'R7'});
        end

        function fitTasksPreMergedResolveTiesStillFillsCorrectly(testCase)
            % resolveTies is threaded through fitTasks (gapFillMode preMerged)
            % -> ftINITFillGaps -> getMinNrFluxes (raven-gecko-parity#104); this
            % exercises that whole chain end to end. Same fixture as
            % fitTasksPreMergedRuns above -- correctness (R7 still gets added
            % back) is what this checks, since that fixture has a unique fill
            % (no tie to actually break); getMinNrFluxesResolveTiesIsDeterministic
            % in tAnalysis.m covers the tie-breaking property itself.
            testCase.assumeMILPSolver();
            testModel      = getTstModel();
            testModelTasks = getTstModelTasks();
            testRxnScores  = getTstModelRxnScores();

            mTempRef = closeModel(testModel);
            mTempRef = removeReactions(mTempRef, {'R1';'R8'});
            mTemp    = removeReactions(mTempRef, {'R7'});
            mTemp.id = 'tmp';
            tmpRxnScores = testRxnScores([2;3;4;5;6;7;9;10]);
            evalc(['[~,addedRxnMat] = fitTasks(mTemp,mTempRef,[],' ...
                '''printOutput'',false,''rxnScores'',min(tmpRxnScores,-0.1),' ...
                '''taskStructure'',testModelTasks,''gapFillMode'',''preMerged'',' ...
                '''params'',struct(),''verbose'',false,''resolveTies'',true);']);
            testCase.verifyEqual(mTempRef.rxns(any(addedRxnMat,2)), {'R7'});
        end

        function fitTasksPreMergedReportsUnfillableTask(testCase)
            % R7 is the only producer of e[s], so removing it from the
            % reference model as well leaves the task unfillable. R8 is kept
            % so that e[s] still takes part in a reaction, and R9 is withheld
            % from the model so the reference has something to offer: the
            % MILP is then the thing that has to report the task as
            % unfillable rather than as "Added 0 reaction(s)".
            testCase.assumeMILPSolver();
            testModel      = getTstModel();
            testModelTasks = getTstModelTasks();
            testRxnScores  = getTstModelRxnScores();

            mTempRef = closeModel(testModel);
            mTempRef = removeReactions(mTempRef, {'R1';'R7'});
            mTemp    = removeReactions(mTempRef, {'R9'});
            mTemp.id = 'tmp';
            tmpRxnScores = testRxnScores([2;3;4;5;6;8;9;10]);
            lastwarn('');
            evalc(['[~,addedRxnMat,failedTasks] = fitTasks(mTemp,mTempRef,[],' ...
                '''printOutput'',false,''rxnScores'',min(tmpRxnScores,-0.1),' ...
                '''taskStructure'',testModelTasks,''gapFillMode'',''preMerged'');']);
            testCase.verifyTrue(failedTasks(1));
            testCase.verifyFalse(any(addedRxnMat(:)));
            % Every way of failing sets failedTasks, so the reason has to be
            % checked too: without this the test would also pass if the
            % gap-filling attempt threw before it reached the solver.
            testCase.verifySubstring(lastwarn, 'no feasible solution exists');
        end

        function fitTasksPreMergedHandlesEmptyReferenceSet(testCase)
            % The reference model holds nothing the model does not already
            % have, so no reaction can be added. An empty reaction list
            % means "all reactions" to the MILP, which would otherwise
            % report the model itself as a solution to an unfillable task.
            testCase.assumeMILPSolver();
            testModel      = getTstModel();
            testModelTasks = getTstModelTasks();
            testRxnScores  = getTstModelRxnScores();

            mTempRef = closeModel(testModel);
            mTempRef = removeReactions(mTempRef, {'R1';'R7';'R8'});
            mTemp    = mTempRef;
            mTemp.id = 'tmp';
            tmpRxnScores = testRxnScores([2;3;4;5;6;9;10]);
            lastwarn('');
            evalc(['[~,addedRxnMat,failedTasks] = fitTasks(mTemp,mTempRef,[],' ...
                '''printOutput'',false,''rxnScores'',min(tmpRxnScores,-0.1),' ...
                '''taskStructure'',testModelTasks,''gapFillMode'',''preMerged'',' ...
                '''params'',struct(),''verbose'',false);']);
            testCase.verifyTrue(failedTasks(1));
            testCase.verifyFalse(any(addedRxnMat(:)));
            testCase.verifySubstring(lastwarn, 'no feasible solution exists');
        end

        function fitTasksPreMergedWarningHasRealNewlineAndPercent(testCase)
            % The "could not be gap-filled" warning embeds the task
            % id/description; a literal "\n" must become a real newline
            % (not print as the two characters backslash-n), and a "%" in
            % the task id must survive intact.
            testCase.assumeMILPSolver();
            testModel      = getTstModel();
            testModelTasks = getTstModelTasks();
            testModelTasks.id = 'Task 50% test';
            testModelTasks.description = testModelTasks.id;
            testRxnScores  = getTstModelRxnScores();

            mTempRef = closeModel(testModel);
            mTempRef = removeReactions(mTempRef, {'R1';'R7';'R8'});
            mTemp    = mTempRef;
            mTemp.id = 'tmp';
            tmpRxnScores = testRxnScores([2;3;4;5;6;9;10]);
            lastwarn('');
            evalc(['fitTasks(mTemp,mTempRef,[],' ...
                '''printOutput'',false,''rxnScores'',min(tmpRxnScores,-0.1),' ...
                '''taskStructure'',testModelTasks,''gapFillMode'',''preMerged'',' ...
                '''params'',struct(),''verbose'',false);']);
            msg = lastwarn();
            testCase.verifySubstring(msg, 'Task 50% test');
            testCase.verifyFalse(contains(msg, '\n'));
        end

        function ftINITFillGapsReportsTaskNeedingAnOrphanMet(testCase)
            % e[s] takes part only in R7 and R8. With both gone from the
            % reference model no reaction touches e[s] at all, so a task
            % requiring net production of it cannot be filled from that
            % reference. That requirement is carried in b, and simplifyModel
            % must leave the metabolite in place: without its row the MILP
            % solves a problem that no longer asks for e[s], and reports the
            % task as satisfiable.
            testCase.assumeMILPSolver();
            testModel     = getTstModel();
            testRxnScores = getTstModelRxnScores();

            mTempRef = closeModel(testModel);
            mTempRef = removeReactions(mTempRef, {'R1';'R7';'R8'});
            % R5 is dropped from the model only so that the reference offers
            % a candidate reaction to consider; adding it back does not make
            % e[s] reachable.
            mTemp    = removeReactions(mTempRef, {'R5'});
            mTemp.id = 'tmp';
            tmpRxnScores = min(testRxnScores([2;3;4;5;6;9;10]), -0.1);

            tModel = setTaskBounds(mTemp);
            tRef   = setTaskBounds(mTempRef);
            sol = solveLP(tModel);
            testCase.assertEmpty(sol.x, 'the task must start out infeasible');

            evalc(['[addedRxns,~,exitFlag] = ftINITFillGaps(tModel,mTemp,tRef,' ...
                'false,true,tmpRxnScores,struct(),false);']);
            testCase.verifyEqual(exitFlag, -1);
            testCase.verifyEmpty(addedRxns);
        end

        function ftINITFillGapsReportsNoCandidateReactions(testCase)
            % The reference model holds no reaction the model does not
            % already have, so nothing can be added to satisfy the task. An
            % empty reaction list means "all reactions" to the MILP, and its
            % solution is then sized for every reaction rather than for the
            % empty set of candidates it is indexed against.
            testCase.assumeMILPSolver();
            testModel     = getTstModel();
            testRxnScores = getTstModelRxnScores();

            m = closeModel(testModel);
            m = removeReactions(m, {'R1'});   % R7/R8 stay, so e[s] is reachable
            mTemp = m;
            mTemp.id = 'tmp';
            tmpRxnScores = min(testRxnScores(2:10), -0.1);

            % Both sides hold the same reactions, so there is no candidate
            % to add. Only the reference may take up a[s], so only it can
            % satisfy the task.
            tRef   = setTaskBounds(m);
            tModel = setTaskBounds(mTemp);
            tModel.b(strcmp(tModel.mets, 'as'), 1) = 0;
            sol = solveLP(tModel);
            testCase.assertEmpty(sol.x, 'the task must start out infeasible');

            evalc(['[addedRxns,~,exitFlag] = ftINITFillGaps(tModel,mTemp,tRef,' ...
                'false,true,tmpRxnScores,struct(),false);']);
            testCase.verifyEqual(exitFlag, -1);
            testCase.verifyEmpty(addedRxns);
        end

        function ftINITMetabolomicsRuns(testCase)
            % Detected metabolites steer ftINIT towards alternative pathways.
            testCase.assumeMILPSolver();
            testModel  = getTstModel();
            testParams = struct();
            evalc('prepData = prepINITModel(testModel, {}, {''R7'';''R10''}, false, {}, ''s'');');

            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;

            % adding metabolite f favours R9/R10 again (presence of R2 is random)
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,{''f''},getINITSteps(),true,true,testParams,false);']);
            if length(resModel.rxns) == 7
                testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                    {'R1';'R4';'R6';'R7';'R8';'R9';'R10'})));
            else
                testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                    {'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})));
            end

            % metabolites a, e, f give the same result
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,{''f'';''a'';''e''},getINITSteps(),true,true,testParams,false);']);
            if length(resModel.rxns) == 7
                testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                    {'R1';'R4';'R6';'R7';'R8';'R9';'R10'})));
            else
                testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                    {'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})));
            end

            % metabolite b turns on R2 and R3/R5 and turns off R9
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,{''b''},getINITSteps(),true,true,testParams,false);']);
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8'})));

            % testModel5 contains reactions that are not linearly merged
            testModel5 = getTstModel5();
            arrayData.genes  = testModel5.genes;
            arrayData.levels = getExprForRxnScore(getTstModel5RxnScores());

            evalc('prepData5 = prepINITModel(testModel5, {}, {''R7'';''R10''}, false, {}, ''s'');');
            evalc(['resModel = ftINIT(prepData5,arrayData.tissues{1},[],[],' ...
                'arrayData,{},getINITSteps(),true,true,testParams,false);']);
            % a->g->e via R11/R13 (score -2) ties the R2 path (score -2), so
            % which of the two is returned is a tie-break that differs between
            % solvers. Neither route uses R9/R10. isequal rather than strcmp,
            % which errors on two lists of different length.
            testCase.verifyTrue( ...
                isequal(resModel.rxns, {'R1';'R4';'R6';'R7';'R8';'R11';'R13'}) || ...
                isequal(resModel.rxns, {'R1';'R2';'R4';'R6';'R7';'R8'}));

            % adding metabolite g drops R7
            evalc('prepData5 = prepINITModel(testModel5, {}, {''R10''}, false, {}, ''s'');');
            arrayData.levels(7) = getExprForRxnScore(-1.1); % avoid randomness
            evalc(['resModel = ftINIT(prepData5,arrayData.tissues{1},[],[],' ...
                'arrayData,{''g''},getINITSteps(),true,true,testParams,false);']);
            testCase.verifyTrue(all(strcmp(resModel.rxns, ...
                {'R1';'R4';'R6';'R8';'R11';'R13'})));
        end

        function ftINITFullVsThreeStepRuns(testCase)
            % On a larger model the three-step and "full" variants must give
            % consistent results.
            testCase.assumeMILPSolver();
            testModelL          = getTstModelL();
            testModelLGeneScores = getTstModelLGeneScores();
            testParams          = struct();

            arrayDataL.genes     = testModelL.genes;
            arrayDataL.tissues   = {'t1'};
            arrayDataL.levels    = getExprForRxnScore(testModelLGeneScores, 1);
            arrayDataL.threshold = 1;

            evalc('prepDataL = prepINITModel(testModelL, [], {}, false, {}, ''s'');');
            evalc(['mres = ftINIT(prepDataL,arrayDataL.tissues{1},[],[],' ...
                'arrayDataL,[],getINITSteps(),true,true,testParams,false);']);
            evalc(['mres2 = ftINIT(prepDataL,arrayDataL.tissues{1},[],[],' ...
                'arrayDataL,[],getINITSteps([], ''full''),true,true,testParams,false);']);

            expResult = {'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'E1';'E2'; ...
                'E3';'E4';'E5';'E6';'E8';'E9';'R1';'R2';'R3';'R4';'R5';'R6';'R7'; ...
                'R8';'R9';'R12';'R13';'R14';'R15'};
            testCase.verifyTrue(all(contains(mres.rxns,  expResult)));
            testCase.verifyTrue(all(contains(mres2.rxns, expResult)));
        end

        function ftINITSeriesVariantsRun(testCase)
            % Exercises the 2-step series from the paper, including the
            % allowExcretion constraint it relies on.
            testCase.assumeMILPSolver();
            testModel  = getTstModel();
            testParams = struct();
            evalc('prepData = prepINITModel(testModel, {}, {}, false, {}, ''s'');');
            arrayData.genes     = testModel.genes;
            arrayData.tissues   = {'a'};
            arrayData.levels    = getExprForRxnScore(getTstModelRxnScores());
            arrayData.threshold = 1;

            steps = getINITSteps([], '2+1');
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],steps,true,true,testParams,false);']);
            % Same answer as the '1+1' series in ftINITPipelineRuns.
            testCase.verifyEqual(resModel.rxns, {'R1';'R4';'R6';'R8';'R9';'R10'});

            steps = getINITSteps([], '2+0');
            evalc(['resModel = ftINIT(prepData,arrayData.tissues{1},[],[],' ...
                'arrayData,[],steps,true,true,testParams,false);']);
            % '2+0' skips step 3, so the GPR-less transport R2 survives.
            testCase.verifyEqual(resModel.rxns, {'R1';'R2';'R4';'R6';'R8';'R9';'R10'});
        end

        function scoreModelDataPrecedenceReactionPrefersHpaWholeRule(testCase)
            % dataPrecedence 'reaction': any gene of a reaction having HPA
            % data means array data is ignored for the whole reaction, not
            % just that gene. R7 is "G1 or G4", G1 has (low) HPA data and G4
            % only array data; per-gene scoring would let G4's higher array
            % score win, per-reaction scoring must not.
            m = scoringTestModel();
            a = scoringArrayData();
            h = scoringHpaDataLowG1();
            perGene = scoreModel(m, h, a, 't1', ...
                'isozymeScoring', 'max', 'complexScoring', 'max');
            perRxn  = scoreModel(m, h, a, 't1', ...
                'isozymeScoring', 'max', 'complexScoring', 'max', ...
                'dataPrecedence', 'reaction');
            r7 = strcmp(m.rxns, 'R7');

            testCase.verifyEqual(perRxn(r7), -8, 'AbsTol', 1e-12);
            % Per gene, G4's array score is seen too, and it is higher
            testCase.verifyEqual(perGene(r7), 5*log(3), 'AbsTol', 1e-12);
            testCase.verifyNotEqual(perRxn(r7), perGene(r7));
        end

        function scoreModelAverageReducesDownTheRule(testCase)
            % multipleGeneScoring 'average' averages down the grRule, so a
            % complex counts once against its isozymes rather than once per
            % subunit. R5 is "(G1 and G2) or G3".
            m = scoringTestModel();
            a = scoringArrayData();
            scores = scoreModel(m, [], a, 't1', ...
                'isozymeScoring', 'average', 'complexScoring', 'average');
            r5 = strcmp(m.rxns, 'R5');

            g = 5*log([1.5; 2; 2.5]);          % scores of G1, G2, G3
            testCase.verifyEqual(scores(r5), mean([mean(g(1:2)); g(3)]), ...
                'AbsTol', 1e-12);
            testCase.verifyNotEqual(scores(r5), mean(g));   % not a flat mean
        end

        function scoreModelReadsGrRulesNotRxnGeneMat(testCase)
            % A model carrying its gene associations only in grRules (no
            % rxnGeneMat) must still score every reaction correctly.
            m = scoringTestModel();
            a = scoringArrayData();
            noMat = m;
            noMat.rxnGeneMat = [];
            withMat = scoreModel(m, [], a, 't1');
            withoutMat = scoreModel(noMat, [], a, 't1');
            testCase.verifyEqual(withoutMat, withMat, 'AbsTol', 1e-12);
        end

        function scoreModelRejectsUnknownDataPrecedence(testCase)
            m = scoringTestModel();
            a = scoringArrayData();
            testCase.verifyError( ...
                @() scoreModel(m, [], a, 't1', 'dataPrecedence', 'model'), ...
                'RAVEN:badInput');
        end

    end
end

%==========================================================================
% Synthetic test models, illustrated in tINITtestInfo/06.svg and 07.svg.
%==========================================================================

function testModel = getTstModel()
testModel = struct();
testModel.id = 'testModel';
testModel.rxns = {};
testModel.S = [];
testModel.rev = [];
testModel.mets = {'as';'ac';'bc';'cc';'dc';'ec';'es';'fc'};
testModel.metNames = {'a';'a';'b';'c';'d';'e';'e';'f'};
testModel.comps = {'s';'c'};
testModel.compNames = testModel.comps;
testModel.metComps = [1;2;2;2;2;2;1;2];
testModel.genes = {'G1';'G2';'G3';'G4';'G5';'G6';'G7';'G8';'G9';'G10'};
testModel.grRules = {};
testModel.rxnGeneMat = [];

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R10'};
rxnsToAdd.equations = {'=> a[s]';...
    'a[s] <=> a[c]';...
    'a[c] <=> b[c] + c[c]';...
    'a[c] <=> 2 d[c]';...
    'b[c] + c[c] => e[c]';...
    '2 d[c] => e[c]';...
    'e[c] => e[s]';...
    'e[s] =>';...
    'a[c] <=> f[c]';...
    'f[c] <=> e[c]'};
rxnsToAdd.grRules = {'';'';'G3';'G4';'G5';'G6';'G7';'';'G9';'G10'};
evalc('testModel = addRxns(testModel,rxnsToAdd, 3);');
testModel.c = [0;0;0;0;0;0;0;1;0;0];
testModel.ub = repmat(1000,10,1);
testModel.lb = [0;-1000;-1000;-1000;0;0;0;0;-1000;-1000];
testModel.rxnNames = testModel.rxns;
testModel.b = repmat(0,8,1);
end

function testModel = getTstModelTwoPaths()
%a[c] reaches m[c] by Rm, Rc or by Rn, Rb; E, R7 and R8 export e[s]. The
%paths are reversible: simplifyModel in prepINITModel fails on a model
%without a reversible reaction after the first simplification.
testModel = struct();
testModel.id = 'twoPaths';
testModel.rxns = {};
testModel.S = [];
testModel.rev = [];
testModel.mets = {'as';'ac';'xc';'yc';'mc';'ec';'es'};
testModel.metNames = {'a';'a';'x';'y';'m';'e';'e'};
testModel.comps = {'s';'c'};
testModel.compNames = testModel.comps;
testModel.metComps = [1;2;2;2;2;2;1];
testModel.genes = {'Gm';'Gc';'Gn';'Gb';'GE';'G7'};
testModel.grRules = {};
testModel.rxnGeneMat = [];

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R1';'R2';'Rm';'Rc';'Rn';'Rb';'E';'R7';'R8'};
rxnsToAdd.equations = {'=> a[s]';...
    'a[s] <=> a[c]';...
    'a[c] <=> x[c]';...
    'x[c] <=> m[c]';...
    'a[c] <=> y[c]';...
    'y[c] <=> m[c]';...
    'm[c] => e[c]';...
    'e[c] => e[s]';...
    'e[s] =>'};
rxnsToAdd.grRules = {'';'';'Gm';'Gc';'Gn';'Gb';'GE';'G7';''};
evalc('testModel = addRxns(testModel,rxnsToAdd, 3);');
testModel.c = [0;0;0;0;0;0;0;0;1];
testModel.ub = repmat(1000,9,1);
testModel.lb = [0;-1000;-1000;-1000;-1000;-1000;0;0;0];
testModel.rxnNames = testModel.rxns;
testModel.b = zeros(7,1);
end

function testModelRxnScores = getTstModelRxnScores()
testModelRxnScores = [-2;-2;-1;7;0.5;0.5;-1;-2;-3;3.5];
end

function model = setTaskBounds(model)
% The b bounds that fitTasks (gapFillMode 'preMerged') derives from getTstModelTasks:
% a[s] may be taken up freely, e[s] has to leave at exactly one unit.
model.b = zeros(numel(model.mets), 2);
model.b(strcmp(model.mets, 'as'), :) = [-inf 0];
model.b(strcmp(model.mets, 'es'), :) = [1 1];
end

function testModelTasks = getTstModelTasks()
testModelTasks = struct();
testModelTasks.id = 'Gen e[s] from a[s]';
testModelTasks.description = 'Gen e[s] from a[s]';
testModelTasks.shouldFail = false;
testModelTasks.printFluxes = false;
testModelTasks.comments = '';
testModelTasks.inputs = {'a[s]'};
testModelTasks.LBin = 0;
testModelTasks.UBin = inf;
testModelTasks.outputs = {'e[s]'};
testModelTasks.LBout = 1;
testModelTasks.UBout = 1;
testModelTasks.equations = {};
testModelTasks.LBequ = [];
testModelTasks.UBequ = [];
testModelTasks.changed = {};
testModelTasks.LBrxn = {};
testModelTasks.UBrxn = {};
end

function testModel2 = getTstModel2()
testModel2 = struct();
testModel2.id = 'testModel2';
testModel2.rxns = {};
testModel2.S = [];
testModel2.rev = [];
testModel2.mets = {'a';'b'};
testModel2.metNames = {'a';'b'};
testModel2.comps = {'s'};
testModel2.compNames = testModel2.comps;
testModel2.metComps = [1;1];
testModel2.genes = {'G1';'G2';'G3';'G4'};
testModel2.grRules = {};
testModel2.rxnGeneMat = [];

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R1';'R2';'R3';'R4'};
rxnsToAdd.equations = {'a[s] <=>';...
    'a[s] => b[s]';...
    'a[s] <=> b[s]';...
    'b[s] =>'};
rxnsToAdd.grRules = testModel2.genes;
evalc('testModel2 = addRxns(testModel2,rxnsToAdd,3,'''',true,true);');
testModel2.c = [0;0;0;1];
testModel2.ub = repmat(1000,4,1);
testModel2.lb = [-1000;0;-1000;0];
testModel2.rxnNames = testModel2.rxns;
testModel2.b = zeros(2,1);
end

function testModel4 = getTstModel4()
testModel4 = getTstModel2();

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R5';'R6';'R7';'R8';'R9';'R10';'R11'};
rxnsToAdd.equations = {'5 a[s] <=> 5 d[s]';...
    'e[s] <=> d[s]';
    'f[s] + g[s] <=> e[s]';
    'b[s] <=> f[s]';
    'h[s] <=> g[s]';
    'h[s] =>';
    'e[s] => g[s]'};
rxnsToAdd.grRules = {'G5';'G6';'G7';'G8';'G9';'G10';'G11'};
[~,testModel4] = evalc("addRxns(testModel4,rxnsToAdd, 3, [], true, true);");
end

function testModel5 = getTstModel5()
testModel = getTstModel();

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R11';'R12';'R13';'R14'};
rxnsToAdd.equations = {'a[c] <=> g[c]';...
    'a[c] <=> g[c]';...
    'g[c] <=> e[c]';...
    'g[c] <=> e[c]'};
rxnsToAdd.grRules = {'G11';'G12';'G13';'G14'};
[~,testModel5] = evalc("addRxns(testModel,rxnsToAdd, 3, [], true, true);");
end

function testModel5RxnScores = getTstModel5RxnScores()
testModel5RxnScores = [getTstModelRxnScores();-1;-1.5;-1;-1.5];
end

function testModelL = getTstModelL()
testModelL = struct();
testModelL.id = 'testModel';
testModelL.rxns = {};
testModelL.S = [];
testModelL.rev = [];
testModelL.metNames = {'e1';'e2';'e3';'e4';'e5';'e6';'e7';'e8';'e9';'e1';'e2';'e3';'e4';'e5';'e6';'e7';'e8';'e9';'x1';'x2';'x3';'x4';'x5';'x6';'x7';'x8';'x9';'x10';'x11'};
testModelL.comps = {'s';'c'};
testModelL.compNames = testModelL.comps;
testModelL.metComps = [1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2];
testModelL.mets = strcat(testModelL.metNames, testModelL.comps(testModelL.metComps));

testModelL.grRules = {};
testModelL.rxnGeneMat = [];

testModelL.genes = {'Ge1';'Ge2';'Ge4';'Ge5';'Ge7';'Ge9'; 'Gr1';'Gr2';'Gr3';'Gr5';'Gr6';'Gr7';'Gr8';'Gr9';'Gr10';'Gr11';'Gr12';'Gr14';'Gr15'};

testModelL.ub = [];
testModelL.lb = [];

rxnsToAdd = struct();
rxnsToAdd.rxns = {  'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'E1';'E2';'E2b';'E3';'E4';'E5';'E6';'E7';'E8';'E9';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R10';'R11';'R12';'R13';'R14';'R15'};
rxnsToAdd.grRules = {'';  '';  '';  '';  '';  '';  '';  '';  ''; 'Ge1';'Ge2';'';'';'Ge4';'Ge5';'';'Ge7';'';'Ge9'; 'Gr1';'Gr2';'Gr3';'';'Gr5';'Gr6';'Gr7';'Gr8';'Gr9';'Gr10';'Gr11';'Gr12';'';'Gr14';'Gr15'};
rxnsToAdd.equations = {'e1[s] <=>';...
    'e2[s] <=>';...
    'e3[s] <=>';...
    'e4[s] <=>';...
    'e5[s] <=>';...
    'e6[s] <=>';...
    'e7[s] <=>';...
    'e8[s] <=>';...
    'e9[s] <=>';...
    'e1[s] <=> e1[c]';...
    'e2[s] <=> e2[c]';...
    'e2[s] <=> e2[c]';... %b variant
    'e3[s] <=> e3[c]';...
    'e4[s] <=> e4[c]';...
    'e5[s] <=> e5[c]';...
    'e6[s] <=> e6[c]';...
    'e7[s] <=> e7[c]';...
    'e8[s] <=> e8[c]';...
    'e9[s] <=> e9[c]';...
    'e1[c] + e2[c] <=> x1[c]';... %R1
    'e1[c] + e3[c] => x2[c] + x3[c]';... %R2
    'e4[c] + x3[c] => x4[c] + x5[c]';... %R3
    'e5[c] + e6[c] + x4[c] => 2 x2[c] + x6[c]';... %R4
    'x1[c] + x2[c] <=> x7[c] + 2 x8[c]';... %R5
    'x2[c] + x8[c] => x3[c] + x9[c]';... %R6
    'x4[c] <=> x9[c]';... %R7
    'x5[c] <=> x9[c]';... %R8
    'x6[c] <=> x10[c]';... %R9
    'x6[c] <=> x11[c]';... %R10
    'x10[c] + 2 x11[c] => e7[c]';... %R11
    'x9[c] + x10[c] <=> e8[c]';... %R12
    'x7[c] + x8[c] + x9[c] => e9[c]';... %R13
    'x6[c] => x9[c]';... %R14
    'x3[c] => x9[c]'... %R15
    };
evalc('testModelL = addRxns(testModelL,rxnsToAdd, 3);');
testModelL.c = [0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
testModelL.rxnNames = testModelL.rxns;
testModelL.b = repmat(0,length(testModelL.mets),1);
end

function testModelLGeneScores = getTstModelLGeneScores()
%testModelL.genes = {'Ge1';'Ge2';'Ge4';'Ge5';'Ge7';'Ge9'; 'Gr1';'Gr2';'Gr3';'Gr5';'Gr6';'Gr7';'Gr8';'Gr9';'Gr10';'Gr11';'Gr12';'Gr14';'Gr15'};
testModelLGeneScores = [3; -1;   8;    6;    -5;    5;     4;    5;    2;    3;    6;    1;    3;    1;    -3;    1;     3;      1;    2];
end

%==========================================================================
% Fixtures for the scoreModel* tests: a small model whose grRules cover
% every rule shape the scoring treats differently (none, single, OR, AND,
% nested, mixed-data-source OR), and matching array/HPA data.
%==========================================================================

function m = scoringTestModel()
m = struct();
m.id = 'scoringTest';
m.rxns = {}; m.S = []; m.rev = [];
m.mets     = {'ac'; 'bc'; 'cc'; 'dc'; 'ec'};
m.metNames = {'a'; 'b'; 'c'; 'd'; 'e'};
m.comps = {'c'}; m.compNames = m.comps;
m.metComps = [1; 1; 1; 1; 1];
m.genes = {'G1'; 'G2'; 'G3'; 'G4'};
m.grRules = {}; m.rxnGeneMat = [];
r = struct();
r.rxns = {'R1'; 'R2'; 'R3'; 'R4'; 'R5'; 'R6'; 'R7'};
r.equations = {'=> a[c]'; 'a[c] => b[c]'; 'b[c] => c[c]'; ...
    'c[c] => d[c]'; 'a[c] => d[c]'; 'd[c] => e[c]'; 'e[c] =>'};
r.grRules = {''; 'G1'; 'G1 or G2'; 'G1 and G2'; ...
    '(G1 and G2) or G3'; 'G4'; 'G1 or G4'};
evalc('m = addRxns(m, r, 3);');
m.c = zeros(7, 1);
m.lb = zeros(7, 1);
m.ub = repmat(1000, 7, 1);
m.rxnNames = m.rxns;
m.b = zeros(5, 1);
evalc('[m.grRules, m.rxnGeneMat] = standardizeGrRules(m, true);');
end

function a = scoringArrayData()
% Distinct, uncapped scores: with a threshold of 1 the score of each gene
% is 5*log(level), i.e. 2.03, 3.47, 4.58 and 5.49.
a = struct();
a.genes     = {'G1'; 'G2'; 'G3'; 'G4'};
a.tissues   = {'t1'; 't2'};
a.celltypes = {'ct1'; 'ct2'};
a.levels    = [1.5 1; 2 1; 2.5 1; 3 1];
a.threshold = ones(4, 1);
end

function h = scoringHpaDataLowG1()
% Only G1 has HPA data, and its level is the lowest one, so any array
% score for the other genes outranks it.
h = struct();
h.genes      = {'G1'};
h.tissues    = {'t1'};
h.celltypes  = {'ct1'};
h.levels     = {'High', 'Medium', 'Low', 'None'};
h.gene2Level = 4;
end
