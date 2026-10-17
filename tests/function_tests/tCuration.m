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
