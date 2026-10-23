classdef tCuration < RavenTestCase
% tCuration  Tests for the table-driven curation functions in curation/.
%
%   curateModelFromTables reads its mets/genes/rxns input from tab-delimited
%   files; the tests below write those files to a temporary directory.

    methods (Test)

        function curateModelFromTablesKeepsSupersetReaction(testCase)
            % A model reaction that uses the candidate's metabolites *and*
            % others is a different reaction, not a duplicate: adding
            % A + B -> C to a model holding A + B -> C + H must leave the
            % existing reaction alone and add the candidate separately.
            model = tCuration.toyModel();
            [coeffsF, infoF] = tCuration.writeRxnTables(testCase, ...
                {'A','B','C'}, [-1 -1 1]);
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',infoF);']);

            testCase.verifyNumElements(out.rxns, 2);
            % The existing reaction keeps its name and its fourth metabolite
            testCase.verifyEqual(out.rxnNames{1}, 'ATP hydrolysis with proton');
            testCase.verifyNumElements(find(out.S(:,1)), 4);
            % The candidate is added as a new reaction
            testCase.verifyEqual(out.rxnNames{2}, 'newRxn');
            testCase.verifyNumElements(find(out.S(:,2)), 3);
        end

        function curateModelFromTablesOverwritesExactDuplicate(testCase)
            % An exact stoichiometric match is still treated as the same
            % reaction and has its annotation overwritten.
            model = tCuration.toyModel();
            [coeffsF, infoF] = tCuration.writeRxnTables(testCase, ...
                {'A','B','C','H'}, [-1 -1 1 1]);
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',infoF);']);

            testCase.verifyNumElements(out.rxns, 1);
            testCase.verifyEqual(out.rxnNames{1}, 'newRxn');
        end

        function curateModelFromTablesReadsMetColumnsByHeader(testCase)
            % Columns are found by header, so a table without an inchi
            % column and in another order must not write its notes over
            % the InChIs; notes and charges reach existing metabolites,
            % creating the fields the model lacks.
            model = tCuration.toyModel();
            model.inchis = {'InChI=A';'InChI=B';'InChI=C';'InChI=H'};
            metsF = tCuration.writeTsv(testCase, { ...
                'comps\tmetNames\tcharge\tmetNotes\tformula', ...
                'c\tA\t-2\tnote A\tC1'});
            evalc('out = curateModelFromTables(model,metsF);');

            testCase.verifyNumElements(out.mets, 4);
            testCase.verifyEqual(out.inchis, model.inchis);
            testCase.verifyEqual(out.metFormulas{1}, 'C1');
            testCase.verifyEqual(out.metNotes, {'note A';'';'';''});
            testCase.verifyEqual(out.metCharges(1), -2);
            testCase.verifyTrue(all(isnan(out.metCharges(2:4))));
        end

        function curateModelFromTablesWritesExistingRxnFields(testCase)
            % Notes, references, confidence scores and MIRIAM columns
            % reach an existing reaction, and grRules are matched to the
            % right reaction when existing and new reactions are mixed.
            % The yeast-GEM tables name the index column "index".
            model = tCuration.toyModel();
            model.genes      = {'g1'};
            model.rxnGeneMat = sparse(1,1);
            coeffsF = tCuration.writeTsv(testCase, { ...
                'index\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', ...
                '2\tnewRxn\tA\tc\t-1', '2\tnewRxn\tB\tc\t-1', ...
                '2\tnewRxn\tC\tc\t1'});
            infoF = tCuration.writeTsv(testCase, { ...
                ['index\trxnNames\tgrRules\tlb\tub\trev\tsubSystems\t' ...
                 'eccodes\trxnNotes\trxnReferences\trxnConfidenceScores\tkegg.reaction'], ...
                '1\tdupRxn\t\t0\t1000\t0\t\t\ta note\tPMID:1\t3\tR00001', ...
                '2\tnewRxn\tg1\t0\t1000\t0\t\t\t\t\t\t'});
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',infoF);']);

            testCase.verifyNumElements(out.rxns, 2);
            testCase.verifyEqual(out.rxnNames{1}, 'dupRxn');
            testCase.verifyEqual(out.rxnNotes{1}, 'a note');
            testCase.verifyEqual(out.rxnReferences{1}, 'PMID:1');
            testCase.verifyEqual(out.rxnConfidenceScores(1), 3);
            testCase.verifyEqual(out.rxnMiriams{1}.name, {'kegg.reaction'});
            testCase.verifyEqual(out.rxnMiriams{1}.value, {'R00001'});
            testCase.verifyEqual(out.grRules{1}, '');
            testCase.verifyEqual(out.grRules{2}, 'g1');
        end

        function curateModelFromTablesWritesExistingGeneFields(testCase)
            % yeast-GEM gene tables use the header geneShortName; it must
            % be read as the short name, not as a MIRIAM namespace.
            model = tCuration.toyModel();
            model.genes      = {'g1'};
            model.rxnGeneMat = sparse(1,1);
            genesF = tCuration.writeTsv(testCase, { ...
                'genes\tgeneShortName\tuniprot', 'g1\tGENE1\tP12345'});
            evalc('out = curateModelFromTables(model,''none'',''genesInfo'',genesF);');

            testCase.verifyEqual(out.geneShortNames, {'GENE1'});
            testCase.verifyEqual(out.geneMiriams{1}.name, {'uniprot'});
            testCase.verifyEqual(out.geneMiriams{1}.value, {'P12345'});
        end

        function curateModelFromTablesOverwritesExistingWithEmptyCells(testCase)
            % By default an empty cell empties the field of an existing
            % metabolite, gene or reaction; an empty charge becomes 0.
            model = tCuration.existingFieldsModel();
            [metsF, genesF, coeffsF, infoF] = tCuration.emptyCellTables(testCase);
            evalc(['out = curateModelFromTables(model,metsF,''genesInfo'',genesF,' ...
                '''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF);']);

            testCase.verifyEqual(out.metFormulas{1}, '');
            testCase.verifyEqual(out.metFormulas{2}, 'F2');
            testCase.verifyEqual(out.metCharges(1:2), [0;5]);
            testCase.verifyEqual(out.geneShortNames, {''});
            testCase.verifyEqual(out.rxnNotes, {''});
            testCase.verifyEqual(out.eccodes, {''});
        end

        function curateModelFromTablesKeepsExistingWithEmptyCellsWhenAsked(testCase)
            % overwriteWithEmpty=false writes only the cells with a value.
            model = tCuration.existingFieldsModel();
            [metsF, genesF, coeffsF, infoF] = tCuration.emptyCellTables(testCase);
            evalc(['out = curateModelFromTables(model,metsF,''genesInfo'',genesF,' ...
                '''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF,''overwriteWithEmpty'',false);']);

            testCase.verifyEqual(out.metFormulas{1}, 'CH4');
            testCase.verifyEqual(out.metFormulas{2}, 'F2');
            testCase.verifyEqual(out.metCharges(1:2), [1;5]);
            testCase.verifyEqual(out.geneShortNames, {'OLD'});
            testCase.verifyEqual(out.rxnNotes, {'old note'});
            testCase.verifyEqual(out.rxnConfidenceScores, 3);
            testCase.verifyEqual(out.eccodes, {'1.1.1.1'});
            testCase.verifyEqual(out.rxnNames, {'ATP hydrolysis with proton'});
            testCase.verifyEqual(out.rxnMiriams{1}.value, {'R00001'});
            testCase.verifyEqual(out.metMiriams{1}.value, {'C00001'});
        end

        function curateModelFromTablesRejectsNonLogicalOverwriteWithEmpty(testCase)
            model = tCuration.toyModel();
            metsF = tCuration.writeTsv(testCase, {'metNames\tcomps', 'A\tc'});
            err = [];
            try
                curateModelFromTables(model,metsF,'overwriteWithEmpty','yes');
            catch err
            end
            testCase.verifyNotEmpty(err);
            testCase.verifyEqual(err.identifier, 'RAVEN:badInput');
        end

        function curateModelFromTablesTreatsEmptyBoundsAsNotGiven(testCase)
            % An empty lb or ub cell leaves the bound of an existing
            % reaction unchanged. A new reaction gets -Inf for a missing lb
            % and the default upper bound for a missing ub.
            model = tCuration.toyModel();
            model.lb = -5;
            model.ub = 7;
            model.annotation.defaultUB = 500;
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', ...
                '2\tnewRxn\tA\tc\t-1', '2\tnewRxn\tB\tc\t-1', ...
                '2\tnewRxn\tC\tc\t1', ...
                '3\tboundRxn\tA\tc\t-1', '3\tboundRxn\tC\tc\t1'});
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tlb\tub\trev', ...
                '1\tdupRxn\t\t\t', '2\tnewRxn\t\t\t', '3\tboundRxn\t-2\t3\t1'});
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',infoF);']);

            testCase.verifyEqual(out.lb, [-5; -Inf; -2]);
            testCase.verifyEqual(out.ub, [7; 500; 3]);
        end

        function curateModelFromTablesRefusesTextAsBound(testCase)
            model = tCuration.toyModel();
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tnewRxn\tA\tc\t-1', '1\tnewRxn\tC\tc\t1'});
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tlb\tub', '1\tnewRxn\tlow\t1000'});
            err = [];
            try
                curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF);
            catch err
            end
            testCase.verifyNotEmpty(err);
            testCase.verifyEqual(err.identifier, 'RAVEN:badInput');
            testCase.verifySubstring(err.message, 'low');
        end

        function curateModelFromTablesReadsTrueFalseRev(testCase)
            % A rev cell can be written as true or false.
            model = tCuration.toyModel();
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', ...
                '2\tnewRxn\tA\tc\t-1', '2\tnewRxn\tC\tc\t1'});
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tlb\tub\trev', ...
                '1\tdupRxn\t\t\ttrue', '2\tnewRxn\t-3\t4\tfalse'});
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',infoF);']);
            testCase.verifyNumElements(out.rxns, 2);
            testCase.verifyEqual(out.lb(2), -3);
            testCase.verifyEqual(out.ub(2), 4);
        end

        function curateModelFromTablesRefusesComplexBoundAndKeepsOnNaN(testCase)
            % A complex bound is refused; a NaN bound counts as not given.
            model = tCuration.toyModel();
            model.lb = -5;
            model.ub = 7;
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1'});
            complexF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tlb', '1\tdupRxn\t1i'});
            err = [];
            try
                curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',complexF);
            catch err
            end
            testCase.verifyNotEmpty(err);
            testCase.verifyEqual(err.identifier, 'RAVEN:badInput');
            nanF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tlb', '1\tdupRxn\tNaN'});
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',nanF);']);
            testCase.verifyEqual(out.lb, -5);
        end

        function curateModelFromTablesRejectsOverwriteWithEmptyThatIsNotZeroOrOne(testCase)
            model = tCuration.toyModel();
            metsF = tCuration.writeTsv(testCase, {'metNames\tcomps', 'A\tc'});
            for bad = {2, NaN, [true false], {true}}
                err = [];
                try
                    curateModelFromTables(model,metsF,'overwriteWithEmpty',bad{1});
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier, 'RAVEN:badInput');
            end
        end

        function curateModelFromTablesBlankMiriamAndGrRulesCountAsEmpty(testCase)
            % A cell with only blanks leaves the MIRIAM annotation of an
            % existing metabolite when overwriteWithEmpty is false, and
            % never replaces a gene rule.
            model = tCuration.existingFieldsModel();
            model.grRules = {'g1'};
            model.rxnGeneMat = sparse(1,1,1);
            metsF = tCuration.writeTsv(testCase, { ...
                'metNames\tcomps\tkegg.compound', 'A\tc\t '});
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\t\tA\tc\t-1', '1\t\tB\tc\t-1', ...
                '1\t\tC\tc\t1',  '1\t\tH\tc\t1'});
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tgrRules', '1\t\t   '});
            evalc(['out = curateModelFromTables(model,metsF,''rxnsCoeffs'',coeffsF,' ...
                '''rxnsInfo'',infoF,''overwriteWithEmpty'',false);']);
            testCase.verifyEqual(out.metMiriams{1}.value, {'C00001'});
            testCase.verifyEqual(out.grRules, {'g1'});
        end

        function curateModelFromTablesRequiresKeyColumns(testCase)
            model = tCuration.toyModel();
            metsF = tCuration.writeTsv(testCase, {'metNames\tformula', 'A\tC1'});
            testCase.verifyError(@() curateModelFromTables(model,metsF), ?MException);
        end

    end

    methods (Static, Access = private)

        function model = toyModel()
            % Single reaction: A + B -> C + H, all in one compartment.
            model = struct();
            model.id        = 'toy';
            model.rxns      = {'R1'};
            model.rxnNames  = {'ATP hydrolysis with proton'};
            model.mets      = {'m1';'m2';'m3';'m4'};
            model.metNames  = {'A';'B';'C';'H'};
            model.comps     = {'c'};
            model.compNames = {'cytosol'};
            model.metComps  = [1;1;1;1];
            model.S         = sparse([-1;-1;1;1]);
            model.lb        = 0;
            model.ub        = 1000;
            model.rev       = 0;
            model.c         = 0;
            model.b         = zeros(4,1);
            model.genes     = {};
            model.grRules   = {''};
            model.rxnGeneMat= sparse(1,0);
            model.metFormulas = {'';'';'';''};
            model.subSystems  = {{''}};
            model.eccodes     = {''};
            model.rxnNotes    = {''};
            model.rxnReferences = {''};
            model.rxnConfidenceScores = 0;
        end

        function model = existingFieldsModel()
            % toyModel with every curated field filled in.
            model = tCuration.toyModel();
            model.metFormulas = {'CH4';'F2';'';''};
            model.metCharges  = [1;2;3;4];
            model.genes       = {'g1'};
            model.geneShortNames = {'OLD'};
            model.rxnGeneMat  = sparse(1,1);
            model.eccodes     = {'1.1.1.1'};
            model.rxnNotes    = {'old note'};
            model.rxnConfidenceScores = 3;
            model.metMiriams  = {struct('name',{{'kegg.compound'}},'value',{{'C00001'}});[];[];[]};
            model.rxnMiriams  = {struct('name',{{'kegg.reaction'}},'value',{{'R00001'}})};
        end

        function [metsF, genesF, coeffsF, infoF] = emptyCellTables(testCase)
            % Rows for existing entities: A has no formula, charge or MIRIAM
            % value (B has a formula and a charge); the gene has no short
            % name; the reaction has no name, EC code, note, confidence
            % score or MIRIAM value.
            metsF = tCuration.writeTsv(testCase, { ...
                'metNames\tcomps\tformula\tcharge\tkegg.compound', ...
                'A\tc\t\t\t', 'B\tc\tF2\t5\t'});
            genesF = tCuration.writeTsv(testCase, { ...
                'genes\tgeneShortNames', 'g1\t'});
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\t\tA\tc\t-1', '1\t\tB\tc\t-1', ...
                '1\t\tC\tc\t1',  '1\t\tH\tc\t1'});
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\teccodes\trxnNotes\trxnConfidenceScores\tkegg.reaction', ...
                '1\t\t\t\t\t'});
        end

        function [coeffsF, infoF] = writeRxnTables(testCase, metNames, coeffs)
            % Write a one-reaction rxnsCoeffs/rxnsInfo pair for 'newRxn'.
            d = tempname; mkdir(d);
            testCase.addTeardown(@() rmdir(d,'s'));
            coeffsF = fullfile(d,'coeffs.tsv');
            infoF   = fullfile(d,'info.tsv');

            fid = fopen(coeffsF,'w');
            fprintf(fid,'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient\n');
            for i=1:numel(metNames)
                fprintf(fid,'1\tnewRxn\t%s\tc\t%d\n', metNames{i}, coeffs(i));
            end
            fclose(fid);

            fid = fopen(infoF,'w');
            fprintf(fid,['rxnIdx\trxnNames\tgrRules\tlb\tub\trev\tsubSystems\t' ...
                'eccodes\trxnNotes\trxnReferences\trxnConfidenceScores\n']);
            fprintf(fid,'1\tnewRxn\t\t0\t1000\t0\t\t\t\t\t\n');
            fclose(fid);
        end

        function f = writeTsv(testCase, lines)
            % Write tab-delimited lines (with \t escapes) to a temporary file.
            d = tempname; mkdir(d);
            testCase.addTeardown(@() rmdir(d,'s'));
            f = fullfile(d,'table.tsv');
            fid = fopen(f,'w');
            for i=1:numel(lines)
                fprintf(fid,[lines{i} '\n']);
            end
            fclose(fid);
        end

    end
end
