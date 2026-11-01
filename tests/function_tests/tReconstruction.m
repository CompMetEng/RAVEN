classdef tReconstruction < RavenTestCase
% tReconstruction  Tests for the de-novo reconstruction functions in reconstruction/.
%
%   Homology and KEGG reconstruction depend on external binaries (BLAST+,
%   DIAMOND) and on a local KEGG dump / aligners. Those are guarded and skip
%   when unavailable; the self-contained functions are tested directly.

    methods (Test)

        function guessCompositionRuns(testCase)
            evalc('[m2, guessed] = guessComposition(testCase.model);');
            testCase.verifyClass(m2, 'struct');
        end

        function guessCompositionAccountsForStoichiometricCoefficient(testCase)
            % A -> 2 B, with A's formula known (CH4). B's own coefficient is
            % 2, so its atom counts must be half of A's (C0.5H4/2 = C0.5H2),
            % not identical to A's.
            model = tReconstruction.toyModel();
            model.metFormulas = {'CH4'; ''};
            evalc('[m2, guessedFor, couldNotGuess] = guessComposition(model);');
            testCase.verifyEqual(guessedFor, {'B'});
            testCase.verifyEmpty(couldNotGuess);
            testCase.verifyEqual(m2.metFormulas{2}, 'C0.5H2');
        end

        function makeFakeBlastStructureReturnsStruct(testCase)
            % makeFakeBlastStructure requires at least 10 ortholog pairs.
            ol = [testCase.model.genes(1:10), strcat('t_', testCase.model.genes(1:10))];
            bs = makeFakeBlastStructure(ol, 'srcModel', 'tgtOrg');
            testCase.verifyClass(bs, 'struct');
        end

        function getModelFromHomologyBuildsDraft(testCase)
            src = testCase.model; src.id = 'srcModel';
            n = min(20, numel(src.genes));
            ol = [src.genes(1:n), strcat('t_', src.genes(1:n))];
            bs = makeFakeBlastStructure(ol, 'srcModel', 'tgtOrg');
            evalc('draft = getModelFromHomology({src}, bs, ''tgtOrg'');');
            testCase.verifyClass(draft, 'struct');
        end

        function getModelFromHomologyComplexPolicyAndOptions(testCase)
            % A small template: a complex, an isozyme pair, a gene-free
            % reaction, and two template genes that map to one new gene.
            t = struct();
            t.id = 'tmpl'; t.name = 'tmpl';
            t.mets = {'a_c';'b_c';'d_c'}; t.metNames = {'A';'B';'D'};
            t.comps = {'c'}; t.compNames = {'cytosol'}; t.metComps = [1;1;1];
            t.rxns = {'Rcplx';'Riso';'Rfree';'Rdup'};
            t.rxnNames = t.rxns;
            t.S = sparse([-1 0 0 -1; 0 -1 0 1; 1 1 -1 0]);
            t.lb = [0;0;0;0]; t.ub = [1000;1000;1000;1000]; t.c = [0;0;0;0]; t.b = [0;0;0];
            t.rev = [0;0;0;0];
            t.grRules = {'tg1 and tg2'; 'tg3 or tg4'; ''; 'tg3 or tg4'};
            t.genes = {'tg1';'tg2';'tg3';'tg4'};
            t.rxnGeneMat = sparse([1 1 0 0; 0 0 1 1; 0 0 0 0; 0 0 1 1]);
            t.rxnNotes = {'n1';'n2';'n3';'n4'};
            t.rxnConfidenceScores = [4;4;4;4];
            t.version = '9.9.9';
            % makeFakeBlastStructure wants at least 10 pairs; the extra ones
            % are genes the template does not have
            pairs = [{'tg1','ng1'; 'tg3','ng3'; 'tg4','ng3'}; ...
                     [strcat('x', string(1:7))', strcat('nx', string(1:7))']];
            pairs = cellstr(pairs);
            bs = makeFakeBlastStructure(pairs, 'tmpl', 'new');
            evalc('flag = getModelFromHomology({t}, bs, ''new'', ''minLen'', 0, ''minIde'', 0);');
            evalc('keep = getModelFromHomology({t}, bs, ''new'', ''minLen'', 0, ''minIde'', 0, ''complexPolicy'', ''keep'', ''keepGeneFree'', true, ''preserveNotes'', true);');
            evalc('drop = getModelFromHomology({t}, bs, ''new'', ''minLen'', 0, ''minIde'', 0, ''complexPolicy'', ''drop'');');
            % flag: unmapped subunit kept as OLD_; gene-free reaction dropped
            testCase.verifyTrue(contains(flag.grRules{strcmp(flag.rxns,'Rcplx')}, 'OLD_tmpl_tg2'));
            testCase.verifyFalse(ismember('Rfree', flag.rxns));
            % two template genes mapping to one new gene appear once
            testCase.verifyEqual(flag.grRules{strcmp(flag.rxns,'Riso')}, 'ng3');
            % keep: subunit dropped, gene-free kept, template notes kept
            testCase.verifyEqual(keep.grRules{strcmp(keep.rxns,'Rcplx')}, 'ng1');
            testCase.verifyTrue(ismember('Rfree', keep.rxns));
            testCase.verifyEqual(keep.rxnNotes{strcmp(keep.rxns,'Rcplx')}, 'n1');
            testCase.verifyEqual(keep.rxnConfidenceScores(strcmp(keep.rxns,'Rcplx')), 4);
            % drop: the incomplete complex is removed
            testCase.verifyFalse(ismember('Rcplx', drop.rxns));
            testCase.verifyTrue(ismember('Riso', drop.rxns));
            % the template version is not inherited
            testCase.verifyFalse(isfield(flag, 'version') && strcmp(flag.version, '9.9.9'));
            testCase.verifyEqual(flag.rxnNotes{1}, 'Included by getModelFromHomology');
        end

        function getModelFromHomologyComplexPolicyOnNestedRules(testCase)
            % An isozyme group whose members all lack an ortholog leaves one flagged
            % subunit; "drop" removes every reaction that keeps a flagged gene.
            % Columns: rule, mapped template genes, policy, genes in the rule ({} = reaction removed).
            cases = {'(A or B) and C', {'C'}, 'flag', {'OLD_T_B', 'nC'}; ...
                     '(A or B) and C', {'C'}, 'keep', {'nC'}; ...
                     '(A or B) and C', {'C'}, 'drop', {}; ...
                     'A or (B and C)', {'A', 'C'}, 'drop', {}; ...
                     '(A and B) or C', {'C'}, 'drop', {}};
            for k = 1:size(cases, 1)
                t = tReconstruction.homologyTemplate('T', {'R1'}, cases(k, 1));
                pairs = cell(numel(cases{k, 2}), 2);
                for q = 1:numel(cases{k, 2})
                    pairs(q, :) = {cases{k, 2}{q}, ['n' cases{k, 2}{q}]};
                end
                draft = tReconstruction.homologyDraft({t}, {pairs}, cases{k, 3}, {});
                if isempty(cases{k, 4})
                    testCase.verifyFalse(ismember('R1', draft.rxns), cases{k, 1});
                else
                    testCase.verifyEqual(sort(draft.genes), sort(cases{k, 4}(:)), cases{k, 1});
                    testCase.verifyFalse(contains(draft.grRules{strcmp(draft.rxns, 'R1')}, ' or '), cases{k, 1});
                end
            end
            % the flagged subunit stays next to its isozyme in the other branch
            t = tReconstruction.homologyTemplate('T', {'R1'}, {'A or (B and C)'});
            draft = tReconstruction.homologyDraft({t}, {{'A', 'nA'; 'C', 'nC'}}, 'flag', {});
            testCase.verifyTrue(contains(draft.grRules{1}, 'OLD_T_B'));
            testCase.verifyTrue(contains(draft.grRules{1}, 'nA'));
            draft = tReconstruction.homologyDraft({t}, {{'A', 'nA'; 'C', 'nC'}}, 'keep', {});
            testCase.verifyTrue(contains(draft.grRules{1}, 'nA'));
            testCase.verifyTrue(contains(draft.grRules{1}, 'nC'));
            testCase.verifyFalse(contains(draft.grRules{1}, 'OLD_'));
        end

        function getModelFromHomologyPreferredOrderDropsSharedGenes(testCase)
            % A template gene that maps to a gene already taken by an earlier model is
            % removed from the later model altogether, including its other mappings.
            t1 = tReconstruction.homologyTemplate('T1', {'R_T1'}, {'Y'});
            t2 = tReconstruction.homologyTemplate('T2', {'R_T2'}, {'X'});
            pairs = {{'Y', 'N1'}, {'X', 'N1'; 'X', 'N2'}};
            draft = tReconstruction.homologyDraft({t1, t2}, pairs, 'flag', {'T1', 'T2'});
            testCase.verifyEqual(draft.rxns, {'R_T1'});
            draft = tReconstruction.homologyDraft({t1, t2}, pairs, 'flag', {'T2', 'T1'});
            testCase.verifyEqual(draft.rxns, {'R_T2'});
            testCase.verifyEqual(sort(draft.genes), {'N1'; 'N2'});
        end

        function getModelFromHomologyPreferredOrderWithoutSharedGenes(testCase)
            % Templates with the same reaction id and no shared gene: the order does not
            % change the result, and no empty gene list reaches removeGenes.
            t1 = tReconstruction.homologyTemplate('T1', {'R1'}, {'Y'});
            t2 = tReconstruction.homologyTemplate('T2', {'R1'}, {'X'});
            pairs = {{'Y', 'N1'}, {'X', 'N2'}};
            none = tReconstruction.homologyDraft({t1, t2}, pairs, 'flag', {});
            ordered = tReconstruction.homologyDraft({t1, t2}, pairs, 'flag', {'T2', 'T1'});
            testCase.verifyEqual(sort(ordered.rxns), sort(none.rxns));
            testCase.verifyEqual(numel(ordered.rxns), 2);
            testCase.verifyEqual(ordered.grRules{strcmp(ordered.rxns, 'R1')}, 'N1');
        end

        function getModelFromHomologyRejectsGeneIdsThatSplit(testCase)
            src = testCase.model; src.id = 'srcModel';
            newGenes = cellstr(strcat('sp|P', string(1:10)', '|X_BUG'));
            bs = makeFakeBlastStructure([src.genes(1:10), newGenes], 'srcModel', 'tgtOrg');
            testCase.verifyError(@() getModelFromHomology({src}, bs, 'tgtOrg'), 'RAVEN:badInput');
        end

        function getModelFromHomologyRejectsUnknownComplexPolicy(testCase)
            src = testCase.model; src.id = 'srcModel';
            bs = makeFakeBlastStructure([src.genes(1:10), strcat('t_', src.genes(1:10))], 'srcModel', 'tgtOrg');
            testCase.verifyError(@() getModelFromHomology({src}, bs, 'tgtOrg', 'complexPolicy', 'maybe'), ?MException);
        end

        function hitsFromSuppOrthBothDirectionsAndPlaceholders(testCase)
            file = tReconstruction.suppOrthTable();
            bs = hitsFromSuppOrth(file, 'bug');
            ids = strcat({bs.fromId}, '>', {bs.toId});
            testCase.verifyEqual(sort(ids), sort({'other>bug','templateGEM>bug','bug>other','bug>templateGEM'}));
            fwd = bs(strcmp(ids, 'bug>templateGEM')); rev = bs(strcmp(ids, 'templateGEM>bug'));
            testCase.verifyEqual(sort(strcat(fwd.fromGenes, '-', fwd.toGenes)), ...
                sort(strcat(rev.toGenes, '-', rev.fromGenes)));
            testCase.verifyTrue(all(fwd.evalue == 0 & fwd.identity == 100 & fwd.aligLen == 1000 & fwd.ppos == 100));
            % ids such as NA, 0123 and UniProt headers stay text
            other = bs(strcmp(ids, 'bug>other'));
            testCase.verifyEqual(sort(other.fromGenes), {'0123';'ng1'});
            testCase.verifyTrue(all(ismember({'NA';'sp|P1|NG4_BUG'}, fwd.fromGenes)));
            % no alignment hit: bitscore 0
            testCase.verifyEqual(fwd.bitscore(strcmp(fwd.toGenes, 'tg3')), 0);
            testCase.verifyEqual(fwd.bitscore(strcmp(fwd.toGenes, 'tg1')), 512);
        end

        function hitsFromSuppOrthFilters(testCase)
            file = tReconstruction.suppOrthTable();
            pairs = @(bs) sort(strcat(bs(strcmp({bs.fromId}, 'bug') & strcmp({bs.toId}, 'templateGEM')).fromGenes, '-', ...
                bs(strcmp({bs.fromId}, 'bug') & strcmp({bs.toId}, 'templateGEM')).toGenes));
            testCase.verifyEqual(pairs(hitsFromSuppOrth(file, 'bug', 'minSupport', 2)), {'ng1-tg1';'ng2-tg2'});
            testCase.verifyEqual(pairs(hitsFromSuppOrth(file, 'bug', 'requireTools', {'broccoli'})), {'ng1-tg1';'ng2-tg3'});
            testCase.verifyEqual(pairs(hitsFromSuppOrth(file, 'bug', 'maxEvalue', 1e-85)), {'NA-tg5';'ng1-tg1'});
            p = pairs(hitsFromSuppOrth(file, 'bug', 'minIdentity', 50));
            testCase.verifyFalse(any(ismember({'ng2-tg3';'sp|P1|NG4_BUG-tg4'}, p)));
            testCase.verifyTrue(ismember('NA-tg5', p));
            testCase.verifyError(@() hitsFromSuppOrth(file, 'bug', 'requireTools', {'oma'}), 'RAVEN:badInput');
            testCase.verifyError(@() hitsFromSuppOrth(file, 'bug', 'minSupport', 4), 'RAVEN:badInput');
            testCase.verifyError(@() hitsFromSuppOrth(file, 'yeast'), 'RAVEN:badInput');
        end

        function hitsFromSuppOrthSpeciesToModel(testCase)
            file = tReconstruction.suppOrthTable();
            bs = hitsFromSuppOrth(file, 'bug', 'speciesToModel', {'templateGEM','iTemplate';'bug','bugGEM';'other',''});
            testCase.verifyEqual(unique({bs.fromId}), {'bugGEM','iTemplate'});
            testCase.verifyError(@() hitsFromSuppOrth(file, 'bug', 'speciesToModel', {'templateGEM','t';'bug','b'}), 'RAVEN:badInput');
        end

        function hitsFromSuppOrthProteinToGene(testCase)
            file = tReconstruction.suppOrthTable({ ...
                'templateGEM	bug	tg1	ng1.1	1	F	1	0	0	60.0	1e-60	250.0'
                'templateGEM	bug	tg1	ng1.2	3	F,S,B	1	1	1	75.0	1e-90	350.0'
                'templateGEM	bug	tg2	ng2.1	2	S,B	0	1	1	50.0	1e-40	150.0'});
            map = {'ng1.1','ng1';'ng1.2','ng1';'ng2.1','ng2'};
            testCase.verifyWarning(@() hitsFromSuppOrth(file, 'bug', 'proteinToGene', map), ...
                'RAVEN:hitsFromSuppOrth:unmappedProteins');
            tsv = [tempname '.tsv'];
            fid = fopen(tsv, 'w'); fprintf(fid, '# protein\tgene\n'); m2 = map'; fprintf(fid, '%s\t%s\n', m2{:}); fclose(fid);
            json = [tempname '.json'];
            fid = fopen(json, 'w'); fprintf(fid, '{"ng1.1": "ng1", "ng1.2": "ng1", "ng2.1": "ng2"}'); fclose(fid);
            for m = {map, tsv, json}
                w = warning('off', 'RAVEN:hitsFromSuppOrth:unmappedProteins');
                [bs, info] = hitsFromSuppOrth(file, 'bug', 'proteinToGene', m{1});
                warning(w);
                fwd = bs(strcmp({bs.fromId}, 'bug'));
                testCase.verifyEqual(sort(strcat(fwd.fromGenes, '-', fwd.toGenes)), {'ng1-tg1';'ng2-tg2'});
                row = info(strcmp(info.newGene, 'ng1'), :);
                testCase.verifyEqual([row.support row.identity row.evalue row.bitscore], [3 75 1e-90 350]);
                testCase.verifyEqual(row.supportedBy{1}, 'F,S,B');
            end
            delete(tsv, json);
        end

        function hitsFromSuppOrthThroughGetModelFromHomology(testCase)
            file = tReconstruction.suppOrthTable();
            bs = hitsFromSuppOrth(file, 'bug', 'minSupport', 2, 'speciesToModel', ...
                {'templateGEM','templateGEM';'bug','bug';'other',''});
            t = tReconstruction.suppOrthTemplate();
            for s = 1:3
                evalc('draft = getModelFromHomology({t}, bs, ''bug'', ''strictness'', s);');
                testCase.verifyEqual(sort(draft.rxns), {'R_iso';'R_single'});
                testCase.verifyEqual(draft.grRules(strcmp(draft.rxns, 'R_iso')), {'ng2'});
            end
        end

        function runSuppOrthWithStub(testCase)
            testCase.assumeFalse(ispc, 'SuppOrth runs on Linux and macOS only');
            dir = tempname; mkdir(dir);
            table = tReconstruction.suppOrthTable();
            % A stub that checks its arguments and writes the fixture table
            stub = fullfile(dir, 'suppOrth');
            fid = fopen(stub, 'w');
            fprintf(fid, '#!/bin/sh\necho "$@" > "%s/args.txt"\nwhile [ "$1" != "-o" ]; do shift; done\ncp "%s" "$2/suppOrth.tsv"\necho done\n', dir, table);
            fclose(fid); fileattrib(stub, '+x');
            fa = fullfile(dir, 'bug.fa'); fid = fopen(fa, 'w'); fprintf(fid, '>ng1 x\nMA\n'); fclose(fid);
            fb = fullfile(dir, 'tmpl.fa'); fid = fopen(fb, 'w'); fprintf(fid, '>tg1 y\nMA\n'); fclose(fid);
            out = fullfile(dir, 'out');
            bs = runSuppOrth('bug', fa, {'templateGEM'}, {fb}, out, 'suppOrth', stub, 'threads', 2, ...
                'tools', 'orthofinder,broccoli', 'speciesToModel', {'templateGEM','templateGEM';'bug','bug';'other',''});
            args = strtrim(fileread(fullfile(dir, 'args.txt')));
            testCase.verifySubstring(args, '--labels bug,templateGEM');
            testCase.verifySubstring(args, '--supportLv 2 --identity_filter 80 --alignment diamond -t 2');
            testCase.verifyTrue(isfile(fullfile(out, 'run_supporth.log')));
            testCase.verifyEqual(unique({bs.fromId}), {'bug','templateGEM'});
            % shared ids get a label prefix that is removed again
            fid = fopen(fb, 'w'); fprintf(fid, '>ng1 y\nMA\n'); fclose(fid);
            fid = fopen(stub, 'w');
            fprintf(fid, '#!/bin/sh\nprintf "Species1\\tSpecies2\\tProtein1\\tProtein2\\tLv_support\\tSupportedBy\\n" > "%s/suppOrth.tsv"\nprintf "templateGEM\\tbug\\ttemplateGEM__ng1\\tbug__ng1\\t2\\tF,B\\n" >> "%s/suppOrth.tsv"\n', out, out);
            fclose(fid);
            bs = runSuppOrth('bug', fa, {'templateGEM'}, {fb}, out, 'suppOrth', stub, 'tools', 'orthofinder,broccoli');
            testCase.verifyEqual(bs(1).fromGenes, {'ng1'}); testCase.verifyEqual(bs(1).toGenes, {'ng1'});
            testCase.verifyTrue(isfile(fullfile(out, 'prefixed_inputs', 'templateGEM.faa')));
            % a failing run raises with the end of its output
            fid = fopen(stub, 'w'); fprintf(fid, '#!/bin/sh\necho "orthofinder FAILED"\nexit 3\n'); fclose(fid);
            testCase.verifyError(@() runSuppOrth('bug', fa, {'templateGEM'}, {fb}, out, 'suppOrth', stub, ...
                'tools', 'orthofinder,broccoli'), 'RAVEN:runSuppOrth:failed');
            rmdir(dir, 's');
        end

        function getBlastRunsWhenAvailable(testCase)
            fa1 = fullfile(testCase.ravenRoot,'tests','function_tests','test_data','human_galactosidases.fa');
            fa2 = fullfile(testCase.ravenRoot,'tests','function_tests','test_data','yeast_galactosidases.fa');
            testCase.assumeDependency(exist(fa1,'file')==2, 'test FASTA files');
            ok = true;
            try, evalc('bs = getBlast(''human'', fa1, {''yeast''}, {fa2});'); catch, ok = false; end
            testCase.assumeTrue(ok, 'BLAST+ binary not functional in this environment');
            testCase.verifyClass(bs, 'struct');
        end

        function getDiamondRunsWhenAvailable(testCase)
            fa1 = fullfile(testCase.ravenRoot,'tests','function_tests','test_data','human_galactosidases.fa');
            fa2 = fullfile(testCase.ravenRoot,'tests','function_tests','test_data','yeast_galactosidases.fa');
            testCase.assumeDependency(exist(fa1,'file')==2, 'test FASTA files');
            ok = true;
            try, evalc('bs = getDiamond(''human'', fa1, {''yeast''}, {fa2});'); catch, ok = false; end
            testCase.assumeTrue(ok, 'DIAMOND binary not functional in this environment');
            testCase.verifyClass(bs, 'struct');
        end

        function getKEGGModelForOrganismNeedsData(testCase)
            testCase.assumeFail('Requires keggModel.mat and an HMM library from toolbox-assets.');
        end

        function getModelFromKEGGNeedsData(testCase)
            matFile = fullfile(testCase.ravenRoot,'src','reconstruction', 'kegg', 'keggModel.mat');
            testCase.assumeFalse(exist(matFile, 'file') == 2, ...
                'keggModel.mat is present; skipping the build-from-artefacts path.');
            testCase.assumeFail('Downloads and assembles the full KEGG artefact set from toolbox-assets; not run automatically.');
        end

        function buildGlobalGPRJoinsGenesThroughSharedKO(testCase)
            % Offline unit test of the join at the core of buildGlobalKEGGModel:
            % two reactions sharing a KO, a KO used by two organisms' genes, and
            % a reaction with no KO at all (should end up with no genes).
            rxns = {'R1'; 'R2'; 'R3'};
            koReaction = table({'K1'; 'K1'; 'K2'}, {'R1'; 'R2'; 'R2'}, ...
                'VariableNames', {'ko', 'reaction'});
            organismGeneKO = table({'a'; 'a'; 'b'}, {'g1'; 'g2'; 'g1'}, {'K1'; 'K2'; 'K1'}, ...
                'VariableNames', {'organism', 'gene', 'ko'});
            [genes, rxnGeneMat] = buildGlobalGPR(rxns, koReaction, organismGeneKO);
            testCase.verifyEqual(genes, {'a:g1'; 'a:g2'; 'b:g1'});
            expected = [1 0 1; 1 1 1; 0 0 0];
            testCase.verifyEqual(full(rxnGeneMat), expected);
        end

        function getPhylDistNeedsData(testCase)
            distFile = fullfile(testCase.ravenRoot,'src','reconstruction', 'kegg', 'keggPhylDist.mat');
            testCase.assumeFalse(exist(distFile, 'file') == 2, ...
                'keggPhylDist.mat is present; skipping error-path test.');
            testCase.verifyError(@() getPhylDist(), 'getPhylDist:noData');
        end

    end

    methods (Static)

        function t = homologyTemplate(id, rxns, rules)
            % One-reaction-per-rule template converting a_c to b_c; genes come from the rules.
            n = numel(rxns);
            tokens = cellfun(@(r) regexp(r, '[A-Za-z0-9_]+', 'match'), rules, 'UniformOutput', false);
            genes = unique(setdiff([tokens{:}], {'and', 'or'}), 'stable')';
            t = struct();
            t.id = id; t.name = id;
            t.mets = {'a_c'; 'b_c'}; t.metNames = {'a'; 'b'};
            t.comps = {'c'}; t.compNames = {'cytosol'}; t.metComps = [1; 1];
            t.rxns = rxns(:); t.rxnNames = rxns(:);
            t.S = sparse(repmat([-1; 1], 1, n));
            t.lb = zeros(n, 1); t.ub = 1000*ones(n, 1); t.c = zeros(n, 1);
            t.b = zeros(2, 1); t.rev = zeros(n, 1);
            t.grRules = rules(:);
            t.genes = genes;
            t.rxnGeneMat = sparse(n, numel(genes));
            for i = 1:n
                t.rxnGeneMat(i, ismember(genes, tokens{i})) = 1;
            end
        end

        function draft = homologyDraft(templates, pairsPerTemplate, policy, order)
            % Draft for organism "N" from fixed ortholog pairs, one list per template.
            blast = [];
            for k = 1:numel(templates)
                list = pairsPerTemplate{k};
                % makeFakeBlastStructure wants at least 10 pairs; the padding genes are in no template
                for pad = size(list, 1)+1:10
                    list(pad, :) = {sprintf('pad_%s_%d', templates{k}.id, pad), sprintf('npad_%s_%d', templates{k}.id, pad)};
                end
                fake = makeFakeBlastStructure(list, templates{k}.id, 'N');
                if isempty(blast)
                    blast = fake;
                else
                    blast = [blast, fake]; %#ok<AGROW>
                end
            end
            args = {'complexPolicy', policy, 'minLen', 0, 'minIde', 0};
            if ~isempty(order)
                args = [args, {'preferredOrder', order}];
            end
            evalc('draft = getModelFromHomology(templates, blast, ''N'', args{:});');
        end
        function file = suppOrthTable(rows)
            % SuppOrth table with input order templateGEM, bug, other
            if nargin < 1
                rows = {
                    'templateGEM	bug	tg1	ng1	3	F,S,B	1	1	1	85.2	1e-150	512.0'
                    'templateGEM	bug	tg2	ng2	2	F,S	1	1	0	55.0	1e-80	300.0'
                    'templateGEM	bug	tg3	ng2	1	B	0	0	1			'
                    'templateGEM	bug	tg4	sp|P1|NG4_BUG	1	F	1	0	0	30.0	1e-20	120.0'
                    'templateGEM	bug	tg5	NA	1	S	0	1	0	82.0	1e-90	280.0'
                    'templateGEM	other	tg1	o1	3	F,S,B	1	1	1	90.0	0.0	600.0'
                    'bug	other	ng1	o1	2	F,B	1	0	1	70.0	1e-100	400.0'
                    'bug	other	0123	o2	1	F	1	0	0	81.0	1e-50	200.0'};
            end
            file = [tempname '.tsv'];
            fid = fopen(file, 'w');
            fprintf(fid, 'Species1\tSpecies2\tProtein1\tProtein2\tLv_support\tSupportedBy\torthofinder\tsonicparanoid\tbroccoli\tIdentity\te_val\tBitscore\n');
            fprintf(fid, '%s\n', rows{:});
            fclose(fid);
        end

        function t = suppOrthTemplate()
            t = struct();
            t.id = 'templateGEM'; t.name = 'templateGEM';
            t.mets = {'a_c';'b_c';'d_c'}; t.metNames = {'A';'B';'D'};
            t.comps = {'c'}; t.compNames = {'cytosol'}; t.metComps = [1;1;1];
            t.rxns = {'R_single';'R_iso';'R_cplx'}; t.rxnNames = t.rxns;
            t.S = sparse([-1 0 -1; 1 -1 0; 0 1 1]);
            t.lb = [0;0;0]; t.ub = [1000;1000;1000]; t.c = [0;0;0]; t.b = [0;0;0]; t.rev = [0;0;0];
            t.grRules = {'tg1'; 'tg2 or tg3'; 'tg4 and tg5'};
            t.genes = {'tg1';'tg2';'tg3';'tg4';'tg5'};
            t.rxnGeneMat = sparse([1 0 0 0 0; 0 1 1 0 0; 0 0 0 1 1]);
        end

        function model = toyModel()
            % Single reaction: A -> 2 B, one compartment.
            model = struct();
            model.id          = 'toy';
            model.name        = 'toy';
            model.mets        = {'m1'; 'm2'};
            model.metNames    = {'A'; 'B'};
            model.comps       = {'c'};
            model.compNames   = {'cytosol'};
            model.metComps    = [1; 1];
            model.metFormulas = {''; ''};
            model.rxns        = {'R1'};
            model.rxnNames    = {'R1'};
            model.S           = sparse([-1; 2]);
            model.lb          = 0;
            model.ub          = 1000;
            model.rev         = 0;
            model.c           = 0;
            model.b           = zeros(2, 1);
        end
    end
end
