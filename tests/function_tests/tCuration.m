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

        function curateModelFromTablesRefusesGenesListedTwice(testCase)
            % A gene id listed twice in genesInfo raises RAVEN:badInput, for a
            % new gene and for a gene the model has; distinct ids are added.
            model = tCuration.toyModel();
            model.genes      = {'g1'};
            model.rxnGeneMat = sparse(1,1);
            for rows = {{'g2\ta','g2\tb'},{'g1\ta','g1\tb'},{'g2\ta','g3\tb','g2\t'}}
                genesF = tCuration.writeTsv(testCase, [{'genes\tgeneShortNames'}, rows{1}]);
                err = [];
                try
                    curateModelFromTables(model,'none','genesInfo',genesF);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,'listed more than once');
            end
            genesF = tCuration.writeTsv(testCase, {'genes\tgeneShortNames','g2\ta','g3\tb'});
            evalc('out = curateModelFromTables(model,''none'',''genesInfo'',genesF);');
            testCase.verifyEqual(sort(out.genes),{'g1';'g2';'g3'});
        end

        function curateModelFromTablesRefusesRepeatedRxnIdx(testCase)
            % A rxnIdx used by two rows of rxnsInfo raises RAVEN:badInput and
            % names the index and the rows; 1, 1.0 and 1e0 are one index.
            model = tCuration.toyModel();
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tnewA\tA\tc\t-1', '1\tnewA\tC\tc\t1', ...
                '1\tnewB\tA\tc\t-1', '1\tnewB\tC\tc\t1'});
            for second = {'1','1.0','1e0'}
                infoF = tCuration.writeTsv(testCase, { ...
                    'rxnIdx\trxnNames', '1\tnewA', [second{1} '\tnewB']});
                err = [];
                try
                    curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,'rxnIdx 1 is used by data rows 1, 2');
            end
        end

        function curateModelFromTablesTrimsPaddedQuotedHeaders(testCase)
            % Spaces inside the quotes of a header cell are ignored, so
            % " comps " and "  metNames" name the core columns.
            model = tCuration.toyModel();
            metsF = tCuration.writeTsv(testCase, { ...
                '"  metNames"\t" comps "\t formula ', 'newMet\tc\tH2O'});
            evalc('out = curateModelFromTables(model,''none'',''metsInfo'',metsF);');
            testCase.verifyNumElements(out.mets, numel(model.mets)+1);
            testCase.verifyEqual(out.metFormulas{end}, 'H2O');
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

        function curateModelFromTablesTreatsNaNChargeAsNotGiven(testCase)
            % A NaN charge counts like an empty charge cell: 0 when other
            % charges are given and overwriteWithEmpty is true, unchanged
            % when it is false; a new metabolite gets 0.
            model = tCuration.existingFieldsModel();
            metsF = tCuration.writeTsv(testCase, { ...
                'metNames\tcomps\tcharge', 'A\tc\tNaN', 'B\tc\t5', 'Z\tc\tnan'});
            evalc('out = curateModelFromTables(model,metsF);');
            idxA = find(strcmp(out.metNames,'A'));
            idxZ = find(strcmp(out.metNames,'Z'));
            testCase.verifyEqual(out.metCharges(idxA), 0);
            testCase.verifyEqual(out.metCharges(idxZ), 0);
            testCase.verifyFalse(any(isnan(out.metCharges)));

            evalc('out = curateModelFromTables(model,metsF,''overwriteWithEmpty'',false);');
            testCase.verifyEqual(out.metCharges(idxA), 1);
            testCase.verifyEqual(out.metCharges(strcmp(out.metNames,'B')), 5);
            testCase.verifyEqual(out.metCharges(idxZ), 0);
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
            % A rev cell can be written as true or false in any case. It sets
            % the reversibility of an existing reaction; for a new reaction
            % the reversibility follows from lb.
            model = tCuration.toyModel();
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', ...
                '2\tnewRxn\tA\tc\t-1', '2\tnewRxn\tC\tc\t1'});
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tlb\tub\trev', ...
                '1\tdupRxn\t\t\tTRUE', '2\tnewRxn\t-3\t4\tfalse'});
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',infoF);']);
            testCase.verifyNumElements(out.rxns, 2);
            testCase.verifyEqual(out.rev(1), 1);
            testCase.verifyEqual(out.lb(2), -3);
            testCase.verifyEqual(out.ub(2), 4);
            testCase.verifyEqual(out.rev(2), 1);
            model.rev = 1;
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\trev', '1\tdupRxn\tFalse'});
            coeffs1F = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1'});
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffs1F,''rxnsInfo'',infoF);']);
            testCase.verifyEqual(out.rev(1), 0);
        end

        function curateModelFromTablesRefusesCommasHexAndOtherNumberText(testCase)
            % Only plain decimal numbers, Inf, NaN and true/false are read;
            % a decimal comma, hexadecimal and similar text raise
            % RAVEN:badInput.
            model = tCuration.toyModel();
            coeffsF = tCuration.dupRxnCoeffs(testCase);
            for txt = {'1,5','-0,5','1,000','0x10','1d3','1_0','Infinity','1+2i','1i','--5','abc'}
                infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames\tlb','1\tdupRxn\t' txt{1}});
                testCase.verifyError(@() curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF),'RAVEN:badInput');
                metsF = tCuration.writeTsv(testCase, {'metNames\tcomps\tcharge',['A\tc\t' txt{1}]});
                testCase.verifyError(@() curateModelFromTables(model,metsF),'RAVEN:badInput');
            end
        end

        function curateModelFromTablesReadsPlainDecimalNumbers(testCase)
            model = tCuration.toyModel();
            coeffsF = tCuration.dupRxnCoeffs(testCase);
            texts = {'+5','.5','5.','1e3',' -2.5E-1 ','TRUE','false'};
            expected = [5 0.5 5 1000 -0.25 1 0];
            for i = 1:numel(texts)
                infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames\tlb',['1\tdupRxn\t' texts{i}]});
                evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF);');
                testCase.verifyEqual(out.lb(1), expected(i));
            end
            infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames\tlb\tub','1\tdupRxn\tNaN\t-INF'});
            evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF);');
            testCase.verifyEqual(out.lb(1), 0);
            testCase.verifyEqual(out.ub(1), -Inf);
        end

        function curateModelFromTablesRefusesFractionalAndInfiniteCharge(testCase)
            % A charge is a finite whole number; empty and NaN are not given.
            model = tCuration.toyModel();
            for txt = {'0.5','2.5','1e-3','Inf','-Inf'}
                metsF = tCuration.writeTsv(testCase, {'metNames\tcomps\tcharge',['A\tc\t' txt{1}]});
                err = [];
                try
                    curateModelFromTables(model,metsF);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,'not a finite whole number');
            end
            for txt = {'2','-1','3.0','1e1',' 2 ','0'}
                metsF = tCuration.writeTsv(testCase, {'metNames\tcomps\tcharge',['Anew\tc\t' txt{1}]});
                evalc('out = curateModelFromTables(model,metsF);');
                testCase.verifyEqual(out.metCharges(end), round(str2double(txt{1})));
            end
            metsF = tCuration.writeTsv(testCase, {'metNames\tcomps\tcharge','Anew\tc\tNaN','Bnew\tc\t4'});
            evalc('out = curateModelFromTables(model,metsF);');
            testCase.verifyEqual(out.metCharges(end-1:end), [0; 4]);
        end

        function curateModelFromTablesRefusesEmptyCoefficientInRowWithContent(testCase)
            % A row of the coefficient table that names a reaction or a
            % metabolite needs a coefficient; a row without any content is
            % ignored.
            model = tCuration.toyModel();
            infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn'});
            for coefficient = {'','NaN'}
                coeffsF = tCuration.writeTsv(testCase, { ...
                    'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                    '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                    '1\tdupRxn\tC\tc\t1',  ['1\tdupRxn\tH\tc\t' coefficient{1}]});
                err = [];
                try
                    curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,'coefficient');
            end
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', '\t\t\t\t'});
            evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF);');
            testCase.verifyEqual(numel(out.rxns), 1);
        end

        function curateModelFromTablesRefusesMalformedGeneRules(testCase)
            % A malformed grRules cell raises RAVEN:badInput for existing and
            % new reactions; well-formed rules and empty cells are accepted.
            model = tCuration.toyModel();
            model.genes      = {'G1';'G2';'G3'};
            model.rxnGeneMat = sparse(1,3);
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', ...
                '2\tnewRxn\tA\tc\t-1', '2\tnewRxn\tC\tc\t1'});
            malformed = {'G1 and and (','(G1 or G2','G1 G2','and G1','G1 or','()', ...
                'G1 and (G2 or )','G1)','G1 and or G2','(G1) (G2)','G1 xor G2'};
            for k = 1:numel(malformed)
                for row = {['1\tdupRxn\t' malformed{k}],['2\tnewRxn\t' malformed{k}]}
                    lines = {'rxnIdx\trxnNames\tgrRules','1\tdupRxn\t','2\tnewRxn\t'};
                    lines{2 + strncmp(row{1},'2',1)} = row{1};
                    infoF = tCuration.writeTsv(testCase, lines);
                    err = [];
                    try
                        curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF);
                    catch err
                    end
                    testCase.verifyNotEmpty(err, malformed{k});
                    testCase.verifyEqual(err.identifier,'RAVEN:badInput', malformed{k});
                    testCase.verifySubstring(err.message,'is not a valid gene rule');
                    testCase.verifySubstring(err.message,'data row');
                end
            end
            valid = {'G1','(G1 and G2) or G3','G1 AND G2','(G1)','G1 And G2','G1 & G2','G1 && G2', ...
                'G1 | G2','G1 || G2','G1 and(G2 or G3)','((G1 and G2) or G3)','  G1  ',''};
            for k = 1:numel(valid)
                infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames\tgrRules',['1\tdupRxn\t' valid{k}],['2\tnewRxn\t' valid{k}]});
                evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF);');
                testCase.verifyEqual(numel(out.rxns), 2, valid{k});
            end
        end

        function curateModelFromTablesRefusesRevOtherThanZeroOrOne(testCase)
            % A rev cell is 0, 1, true or false (any case); empty and NaN
            % are not given. Other values raise an error for existing and
            % new reactions.
            model = tCuration.toyModel();
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', ...
                '2\tnewRxn\tA\tc\t-1', '2\tnewRxn\tC\tc\t1'});
            for txt = {'2','0.5','Inf','-1','-Inf'}
                for row = {['1\tdupRxn\t' txt{1}],['2\tnewRxn\t' txt{1}]}
                    lines = {'rxnIdx\trxnNames\trev','1\tdupRxn\t','2\tnewRxn\t'};
                    lines{2 + strncmp(row{1},'2',1)} = row{1};
                    infoF = tCuration.writeTsv(testCase, lines);
                    err = [];
                    try
                        curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF);
                    catch err
                    end
                    testCase.verifyNotEmpty(err);
                    testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                    testCase.verifySubstring(err.message,'not 0, 1, true or false');
                end
            end
            for txt = {'0','1','true','FALSE','NaN',''}
                infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames\trev',['1\tdupRxn\t' txt{1}],['2\tnewRxn\t' txt{1}]});
                evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF);');
                testCase.verifyEqual(numel(out.rxns), 2);
            end
        end

        function curateModelFromTablesRefusesRxnIdxThatIsNotAWholeNumberFromOne(testCase)
            % A reaction index is a whole number of at least 1.
            model = tCuration.toyModel();
            coeffsGood = tCuration.dupRxnCoeffs(testCase);
            infoGood = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn'});
            for txt = {'1.5','Inf','0','-1','-Inf'}
                coeffsF = tCuration.writeTsv(testCase, { ...
                    'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                    '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                    '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', ...
                    [txt{1} '\tother\tA\tc\t-1']});
                err = [];
                try
                    curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoGood);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,'whole number of at least 1');
                infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames',['1\tdupRxn'],[txt{1} '\tother']});
                err = [];
                try
                    curateModelFromTables(model,'none','rxnsCoeffs',coeffsGood,'rxnsInfo',infoF);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,'whole number of at least 1');
            end
            evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsGood,''rxnsInfo'',infoGood);');
            testCase.verifyEqual(numel(out.rxns), 1);
        end

        function curateModelFromTablesRefusesTextInIndexAndCoefficientColumns(testCase)
            % A coefficient or index that is not a number raises
            % RAVEN:badInput.
            model = tCuration.toyModel();
            infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn'});
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\tabc', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1'});
            testCase.verifyError(@() curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF),'RAVEN:badInput');
            coeffsF = tCuration.dupRxnCoeffs(testCase);
            infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','x1\tdupRxn'});
            testCase.verifyError(@() curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF),'RAVEN:badInput');
        end

        function curateModelFromTablesTreatsNaNScoreAsNotGiven(testCase)
            % A NaN confidence score counts like an empty cell, so the
            % existing score is not replaced by NaN.
            model = tCuration.existingFieldsModel();
            coeffsF = tCuration.dupRxnCoeffs(testCase);
            infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames\trxnConfidenceScores','1\tdupRxn\tNaN'});
            for flag = [true false]
                evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF,''overwriteWithEmpty'',flag);');
                testCase.verifyEqual(out.rxnConfidenceScores(1), 3);
            end
        end

        function curateModelFromTablesReadsNumbersWithoutEvaluatingText(testCase)
            % Cell text is read as a number or true/false; no expression is
            % evaluated, and other text raises RAVEN:badInput.
            model = tCuration.toyModel();
            d = tempname; mkdir(d);
            testCase.addTeardown(@() rmdir(d,'s'));
            marker = fullfile(d,'executed');
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1'});
            cells = {sprintf('system(''touch %s'')',marker), 'pi', '1/2', '[3]'};
            for i = 1:numel(cells)
                infoF = tCuration.writeTsv(testCase, { ...
                    'rxnIdx\trxnNames\tlb', ['1\tdupRxn\t' strrep(cells{i},'\','\\')]});
                err = [];
                try
                    curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier, 'RAVEN:badInput');
            end
            testCase.verifyEqual(exist(marker,'file'), 0);
            metsF = tCuration.writeTsv(testCase, {'metNames\tcomps\tcharge', 'A\tc\tpi'});
            err = [];
            try
                curateModelFromTables(model,metsF);
            catch err
            end
            testCase.verifyEqual(err.identifier, 'RAVEN:badInput');
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\trxnConfidenceScores', '1\tdupRxn\t2+1'});
            err = [];
            try
                curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF);
            catch err
            end
            testCase.verifyEqual(err.identifier, 'RAVEN:badInput');
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tlb\tub\trxnConfidenceScores', '1\tdupRxn\t -1e1 \t-Inf\t3'});
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',infoF);']);
            testCase.verifyEqual(out.lb(1), -10);
            testCase.verifyEqual(out.ub(1), -Inf);
            testCase.verifyEqual(out.rxnConfidenceScores(1), 3);
        end

        function curateModelFromTablesBlankCellsCountAsEmptyWhenOverwriting(testCase)
            % A cell with only blanks is an empty cell: it removes the
            % MIRIAM namespace and empties the formula of an existing
            % metabolite.
            model = tCuration.existingFieldsModel();
            metsF = tCuration.writeTsv(testCase, { ...
                'metNames\tcomps\tformula\tkegg.compound', ...
                'A\tc\t" "\t" "', 'B\tc\tF2\tC00002'});
            evalc('out = curateModelFromTables(model,metsF);');
            testCase.verifyEqual(out.metFormulas{1}, '');
            testCase.verifyTrue(isempty(out.metMiriams{1}) || ~isfield(out.metMiriams{1},'name') ...
                || ~any(strcmp(out.metMiriams{1}.name,'kegg.compound')));
            testCase.verifyEqual(out.metMiriams{2}.value, {'C00002'});
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

        function curateModelFromTablesRejectsOverwriteWithEmptyNotZeroOrOne(testCase)
            model = tCuration.toyModel();
            metsF = tCuration.writeTsv(testCase, {'metNames\tcomps', 'A\tc'});
            for bad = {2, NaN, [true false], {true}, complex(1,0)}
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
                'metNames\tcomps\tkegg.compound', 'A\tc\t" "'});
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\t\tA\tc\t-1', '1\t\tB\tc\t-1', ...
                '1\t\tC\tc\t1',  '1\t\tH\tc\t1'});
            infoF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tgrRules', '1\t\t"   "'});
            evalc(['out = curateModelFromTables(model,metsF,''rxnsCoeffs'',coeffsF,' ...
                '''rxnsInfo'',infoF,''overwriteWithEmpty'',false);']);
            testCase.verifyEqual(out.metMiriams{1}.value, {'C00001'});
            testCase.verifyEqual(out.grRules, {'g1'});
        end

        function curateModelFromTablesRefusesRowsWithValuesButNoRxnIdx(testCase)
            % A row of the rxnsInfo or rxnsCoeffs table that has values needs
            % a rxnIdx; a row without any value is ignored.
            model = tCuration.toyModel();
            coeffsGood = tCuration.dupRxnCoeffs(testCase);
            infoGood = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn'});
            for idx = {'','NaN'}
                infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn',[idx{1} '\tskipped']});
                err = [];
                try
                    curateModelFromTables(model,'none','rxnsCoeffs',coeffsGood,'rxnsInfo',infoF);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,'rxnIdx');
                coeffsF = tCuration.writeTsv(testCase, { ...
                    'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                    '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                    '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', ...
                    [idx{1} '\tother\tA\tc\t-1']});
                err = [];
                try
                    curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoGood);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,'rxnIdx');
            end
            infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn','\t'});
            coeffsF = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', '\t\t\t\t7'});
            evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF);');
            testCase.verifyEqual(numel(out.rxns), 1);
        end

        function curateModelFromTablesRequiresMatchingReactionTables(testCase)
            % The reactions of the rxnsInfo and rxnsCoeffs tables match
            % exactly; every mismatch is reported with RAVEN:badInput.
            model = tCuration.toyModel();
            coeffsGood = tCuration.dupRxnCoeffs(testCase);
            infoGood = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn'});
            coeffsExtra = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1', ...
                '2\tnewRxn\tA\tc\t-1', '2\tnewRxn\tC\tc\t1'});
            infoExtra = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn','2\tnewRxn'});
            infoRenamed = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tother'});
            cases = {coeffsExtra, infoGood, 'newRxn'; coeffsGood, infoExtra, 'newRxn'; ...
                coeffsGood, infoRenamed, 'other'};
            for i = 1:size(cases,1)
                err = [];
                try
                    curateModelFromTables(model,'none','rxnsCoeffs',cases{i,1},'rxnsInfo',cases{i,2});
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,cases{i,3});
            end
            evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsExtra,''rxnsInfo'',infoExtra);');
            testCase.verifyEqual(numel(out.rxns), 2);
        end

        function curateModelFromTablesMatchesRxnIdxByValue(testCase)
            % The indices 1, 1.0 and 1e0 of the two reaction tables are the
            % same index.
            model = tCuration.toyModel();
            coeffsF = tCuration.dupRxnCoeffs(testCase);
            for txt = {'1.0','1e0','01','+1'}
                infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames\tlb',['' txt{1} '\tdupRxn\t5']});
                evalc('out = curateModelFromTables(model,''none'',''rxnsCoeffs'',coeffsF,''rxnsInfo'',infoF);');
                testCase.verifyEqual(numel(out.rxns), 1);
                testCase.verifyEqual(out.lb(1), 5);
            end
        end

        function curateModelFromTablesAcceptsTrueFalseOnlyInBoundsAndRev(testCase)
            % The words true and false are numbers in the lb, ub and rev
            % cells only.
            model = tCuration.toyModel();
            coeffsF = tCuration.dupRxnCoeffs(testCase);
            infoGood = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn'});
            for txt = {'true','FALSE'}
                metsF = tCuration.writeTsv(testCase, {'metNames\tcomps\tcharge',['A\tc\t' txt{1}]});
                testCase.verifyError(@() curateModelFromTables(model,metsF),'RAVEN:badInput');
                infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames\trxnConfidenceScores',['1\tdupRxn\t' txt{1}]});
                testCase.verifyError(@() curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF),'RAVEN:badInput');
                infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames',[txt{1} '\tdupRxn']});
                testCase.verifyError(@() curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF),'RAVEN:badInput');
                coeffsBad = tCuration.writeTsv(testCase, { ...
                    'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                    '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                    '1\tdupRxn\tC\tc\t1',  ['1\tdupRxn\tH\tc\t' txt{1}]});
                testCase.verifyError(@() curateModelFromTables(model,'none','rxnsCoeffs',coeffsBad,'rxnsInfo',infoGood),'RAVEN:badInput');
            end
        end

        function curateModelFromTablesRefusesInfiniteAndZeroCoefficient(testCase)
            % A coefficient is a finite number other than 0.
            model = tCuration.toyModel();
            infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn'});
            for txt = {'Inf','-Inf','0','0.0'}
                coeffsF = tCuration.writeTsv(testCase, { ...
                    'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                    '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                    '1\tdupRxn\tC\tc\t1',  ['1\tdupRxn\tH\tc\t' txt{1}]});
                err = [];
                try
                    curateModelFromTables(model,'none','rxnsCoeffs',coeffsF,'rxnsInfo',infoF);
                catch err
                end
                testCase.verifyNotEmpty(err);
                testCase.verifyEqual(err.identifier,'RAVEN:badInput');
                testCase.verifySubstring(err.message,'finite non-zero');
            end
        end

        function curateModelFromTablesTreatsEmptyOverwriteWithEmptyAsDefault(testCase)
            % An empty overwriteWithEmpty means the default, true.
            model = tCuration.existingFieldsModel();
            metsF = tCuration.writeTsv(testCase, {'metNames\tcomps\tformula','A\tc\t'});
            evalc('out = curateModelFromTables(model,metsF,''overwriteWithEmpty'',[]);');
            testCase.verifyEqual(out.metFormulas{1}, '');
        end

        function curateModelFromTablesRaisesBadInputForUnreadableFiles(testCase)
            % Problems with the table files are reported with RAVEN:badInput.
            model = tCuration.toyModel();
            testCase.verifyError(@() curateModelFromTables(model,fullfile(tempdir,'does_not_exist_curate.tsv')),'RAVEN:badInput');
            emptyF = tCuration.writeTsv(testCase, {});
            testCase.verifyError(@() curateModelFromTables(model,emptyF),'RAVEN:badInput');
            metsF = tCuration.writeTsv(testCase, {'metNames\tformula', 'A\tC1'});
            testCase.verifyError(@() curateModelFromTables(model,metsF),'RAVEN:badInput');
            coeffsF = tCuration.dupRxnCoeffs(testCase);
            testCase.verifyError(@() curateModelFromTables(model,'none','rxnsCoeffs',coeffsF),'RAVEN:badInput');
            metsOK = tCuration.writeTsv(testCase, {'metNames\tcomps','A\tc'});
            missingF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames\tmetNames\tcomps','1\tdupRxn\tA\tc'});
            infoF = tCuration.writeTsv(testCase, {'rxnIdx\trxnNames','1\tdupRxn'});
            testCase.verifyError(@() curateModelFromTables(model,metsOK,'rxnsCoeffs',missingF,'rxnsInfo',infoF),'RAVEN:badInput');
            coeffsMetMissing = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient','1\tdupRxn\tNotThere\tc\t-1'});
            testCase.verifyError(@() curateModelFromTables(model,'none','rxnsCoeffs',coeffsMetMissing,'rxnsInfo',infoF),'RAVEN:badInput');
        end

        function curateModelFromTablesRequiresKeyColumns(testCase)
            model = tCuration.toyModel();
            metsF = tCuration.writeTsv(testCase, {'metNames\tformula', 'A\tC1'});
            testCase.verifyError(@() curateModelFromTables(model,metsF), ?MException);
        end

    end

    methods (Static, Access = private)

        function f = dupRxnCoeffs(testCase)
            % Coefficients of A + B -> C + H as reaction 1, "dupRxn".
            f = tCuration.writeTsv(testCase, { ...
                'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient', ...
                '1\tdupRxn\tA\tc\t-1', '1\tdupRxn\tB\tc\t-1', ...
                '1\tdupRxn\tC\tc\t1',  '1\tdupRxn\tH\tc\t1'});
        end

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
