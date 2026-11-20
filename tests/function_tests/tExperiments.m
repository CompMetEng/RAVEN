classdef tExperiments < RavenTestCase
% tExperiments  Tests for the functions in experiments/.
%
%   A toy model with two carbon sources feeds every simulation, so that each
%   expected value follows from the stoichiometry. Glucose and ethanol feed the
%   pool metabolite A_c, which the reaction R turns into biomass and into a
%   byproduct (0.1 per unit of R flux). R needs the enzyme prot_E, whose usage
%   is capped by a protein pool of 1 mg/gDCW, so R is capped at
%   3600*kcat/MW = 3600/100 = 36.

    properties (Constant)
        enzymeCap = 36
        byproductYield = 0.1
        biomassCarbon = 41
    end

    methods (Test)

        %% loadFluxData

        function loadFluxDataMissingFile(testCase)
            testCase.verifyError(@() loadFluxData(tempname), 'loadFluxData:fileNotFound');
        end

        function loadFluxDataEmptyFile(testCase)
            f = tExperiments.writeFile(testCase, '');
            testCase.verifyError(@() loadFluxData(f), 'loadFluxData:emptyFile');
        end

        function loadFluxDataShortHeader(testCase)
            f = tExperiments.writeRows(testCase, {{'only','two'}});
            testCase.verifyError(@() loadFluxData(f), 'loadFluxData:noHeader');
        end

        function loadFluxDataNoExchangeColumns(testCase)
            f = tExperiments.writeRows(testCase, ...
                {{'condition','Ptot','grRate'}, {'c1','0.5','0.4'}});
            testCase.verifyError(@() loadFluxData(f), 'loadFluxData:noExchangeColumns');
        end

        function loadFluxDataRaggedRow(testCase)
            f = tExperiments.writeRows(testCase, {
                {'condition','Ptot','grRate','Glc (EX_a)','O2 (EX_b)'}, ...
                {'c1','0.5','0.4','-1','2'}, ...
                {'c2','0.5','0.4','-1'}});
            testCase.verifyError(@() loadFluxData(f), 'loadFluxData:raggedRow');
        end

        function loadFluxDataOldWeightColumnIsRejected(testCase)
            f = tExperiments.writeRows(testCase, {
                {'Condition','Ptot','grRate','glucose (r_1714)','evotuneRMSEweight'}, ...
                {'c1','0.5','0.4','-10','1'}});
            testCase.verifyError(@() loadFluxData(f), 'loadFluxData:renamedColumn');
        end

        function loadFluxDataMultipleConditions(testCase)
            f = tExperiments.writeRows(testCase, {
                {'condition','Ptot','grRate','glucose (EX_glc)','oxygen (EX_O2)','ethanol (EX_eth)'}, ...
                {'c1','0.5','0.4','-10.0','-2.0','5.0'}, ...
                {'c2','0.6','0.5','-8.5','-3.0','7.0'}});
            fd = loadFluxData(f);
            testCase.verifyEqual(fd.conds, {'c1';'c2'});
            testCase.verifyEqual(fd.Ptot, [0.5;0.6]);
            testCase.verifyEqual(fd.grRate, [0.4;0.5]);
            testCase.verifyEqual(fd.exchFluxes, [-10 -2 5; -8.5 -3 7]);
            testCase.verifyEqual(fd.exchMets, {'glucose','oxygen','ethanol'});
            testCase.verifyEqual(fd.exchRxnIDs, {'EX_glc','EX_O2','EX_eth'});
            testCase.verifyFalse(isfield(fd, 'weight'));
            testCase.verifyFalse(isfield(fd, 'source'));
        end

        function loadFluxDataHeaderForms(testCase)
            f = tExperiments.writeRows(testCase, {
                {'condition','Ptot','grRate','carbon dioxide (EX_co2)','X (EX_a(e))'}, ...
                {'c1','0.5','0.4','5.0','1'}});
            fd = loadFluxData(f);
            testCase.verifyEqual(fd.exchMets, {'carbon dioxide','X'});
            testCase.verifyEqual(fd.exchRxnIDs, {'EX_co2','EX_a(e)'});
        end

        function loadFluxDataHeaderWithoutParenthesesWarns(testCase)
            f = tExperiments.writeRows(testCase, {
                {'condition','Ptot','grRate','EX_glc'}, {'c1','0.5','0.4','-10.0'}});
            testCase.verifyWarning(@() loadFluxData(f), 'loadFluxData:headerForm');
            warning('off', 'loadFluxData:headerForm');
            testCase.addTeardown(@() warning('on', 'loadFluxData:headerForm'));
            fd = loadFluxData(f);
            testCase.verifyEqual(fd.exchMets, {'EX_glc'});
            testCase.verifyEqual(fd.exchRxnIDs, {'EX_glc'});
        end

        function loadFluxDataEmptyAndNonNumericCellsAreNaN(testCase)
            f = tExperiments.writeRows(testCase, {
                {'condition','Ptot','grRate','glc (EX_glc)','o2 (EX_O2)'}, ...
                {'c1','0.5','0.4','-10.0',''}, ...
                {'c2','0.6','0.5','NA','-3.0'}});
            fd = loadFluxData(f);
            testCase.verifyTrue(isnan(fd.exchFluxes(1, 2)));
            testCase.verifyTrue(isnan(fd.exchFluxes(2, 1)));
            testCase.verifyEqual(fd.exchFluxes(1, 1), -10);
            testCase.verifyEqual(fd.exchFluxes(2, 2), -3);
        end

        function loadFluxDataOptionalColumnsInAnyPosition(testCase)
            f = tExperiments.writeRows(testCase, {
                {'condition','Ptot','grRate','source','glc (EX_glc)','weight'}, ...
                {'c1','0.5','0.4','PMID:1','-10.0','1.5'}, ...
                {'c2','0.6','0.5','in-house','-8.5','0'}});
            fd = loadFluxData(f);
            testCase.verifyEqual(fd.weight, [1.5;0]);
            testCase.verifyEqual(fd.source, {'PMID:1';'in-house'});
            testCase.verifyEqual(fd.exchFluxes, [-10;-8.5]);
            testCase.verifyEqual(fd.exchRxnIDs, {'EX_glc'});
        end

        function loadFluxDataSkipsBlankLines(testCase)
            f = tExperiments.writeFile(testCase, ...
                sprintf('condition\tPtot\tgrRate\tglc (EX_glc)\n\nc1\t0.5\t0.4\t-10.0\n\n'));
            fd = loadFluxData(f);
            testCase.verifyEqual(fd.conds, {'c1'});
        end

        %% loadExperimentData

        function loadExperimentDataAllFiles(testCase)
            d = tempname; mkdir(d);
            testCase.addTeardown(@() rmdir(d, 's'));
            tExperiments.writeNamed(d, 'fluxData.tsv', sprintf(['Condition\tPtot\tgrRate\tglucose (r_1714)\tethanol (r_1761)\tweight\tsource\n' ...
                'glucose\tNaN\t0.36\t-13.33\t19.83\t1\tDLKcat\nethanol\tNaN\t0.12\tNaN\t-8.5\t1\tbrenda\n']));
            tExperiments.writeNamed(d, 'maxGrowth.tsv', sprintf(['Condition\tPtot\tgrRate\tglucose (r_1714)\tfructose (r_1709)\tweight\tsource\n' ...
                'glucose\tNaN\t0.41\t-1000\tNaN\t1\tDLKcat\nfructose\tNaN\t0.338\tNaN\t-1000\t1\tDLKcat\n']));
            tExperiments.writeNamed(d, 'zeroFlux.tsv', sprintf('Rxns\nr_1542\nr_1545\nr_1546\n'));
            ed = loadExperimentData(d);
            testCase.verifyEqual(ed.fluxData.conds, {'glucose';'ethanol'});
            testCase.verifyEqual(ed.fluxData.exchRxnIDs, {'r_1714','r_1761'});
            testCase.verifyEqual(ed.fluxData.source, {'DLKcat';'brenda'});
            testCase.verifyEqual(ed.fluxData.weight, [1;1]);
            testCase.verifyEqual(ed.maxGrowth.conds, {'glucose';'fructose'});
            testCase.verifyEqual(ed.zeroFlux, {'r_1542';'r_1545';'r_1546'});
        end

        function loadExperimentDataMissingFiles(testCase)
            d = tempname; mkdir(d);
            testCase.addTeardown(@() rmdir(d, 's'));
            ed = loadExperimentData(d);
            testCase.verifyEmpty(ed.fluxData);
            testCase.verifyEmpty(ed.maxGrowth);
            testCase.verifyEqual(ed.zeroFlux, {});
        end

        function loadExperimentDataPartial(testCase)
            d = tempname; mkdir(d);
            testCase.addTeardown(@() rmdir(d, 's'));
            tExperiments.writeNamed(d, 'fluxData.tsv', ...
                sprintf('Condition\tPtot\tgrRate\tglucose (r_1714)\nglucose\tNaN\t0.36\t-13.33\n'));
            ed = loadExperimentData(d);
            testCase.verifyNotEmpty(ed.fluxData);
            testCase.verifyEmpty(ed.maxGrowth);
            testCase.verifyEqual(ed.zeroFlux, {});
        end

        %% computeCarbonWeights

        function computeCarbonWeightsCountsCarbon(testCase)
            model = tExperiments.toyModel();
            w = computeCarbonWeights(model, {'EX_glc','EX_eth','EX_byp','biomass'}, 'biomass');
            testCase.verifyEqual(w('EX_glc'), 6);
            testCase.verifyEqual(w('EX_eth'), 2);
            testCase.verifyEqual(w('EX_byp'), 1);
            testCase.verifyEqual(w('biomass'), tExperiments.biomassCarbon);
        end

        function computeCarbonWeightsZeroCarbonIsOne(testCase)
            % O2 has no carbon; it still gets weight 1 so that its mismatches count.
            model = tExperiments.toyModel();
            w = computeCarbonWeights(model, {'EX_o2'}, 'biomass');
            testCase.verifyEqual(w('EX_o2'), 1);
        end

        function computeCarbonWeightsWithoutFormulas(testCase)
            model = rmfield(tExperiments.toyModel(), 'metFormulas');
            testCase.verifyWarning(@() computeCarbonWeights(model, {'EX_glc'}, 'biomass'), ...
                'computeCarbonWeights:noFormulas');
            warning('off', 'computeCarbonWeights:noFormulas');
            testCase.addTeardown(@() warning('on', 'computeCarbonWeights:noFormulas'));
            w = computeCarbonWeights(model, {'EX_glc'}, 'biomass');
            testCase.verifyEqual(w('EX_glc'), 1);
        end

        function computeCarbonWeightsMissingReaction(testCase)
            model = tExperiments.toyModel();
            testCase.verifyError(@() computeCarbonWeights(model, {'EX_absent'}, 'biomass'), ...
                'computeCarbonWeights:missingReaction');
        end

        %% simulateConditions

        function simulateConditionsEnzymeVersusCarbonLimited(testCase)
            model = tExperiments.toyModel();
            sims = simulateConditions(model, tExperiments.fluxDataToy(20), true, {'EX_byp'}, 'biomass');
            testCase.verifySize(sims, [2 1]);
            % Glucose: enzyme-limited growth of 36, the whole uptake is used.
            testCase.verifyTrue(sims(1).feasible);
            testCase.verifyEqual(sims(1).growth, tExperiments.enzymeCap, 'AbsTol', 1e-6);
            testCase.verifyEqual(sims(1).exchFluxes(strcmp(sims(1).exchRxnIDs, 'EX_glc')), ...
                -tExperiments.enzymeCap, 'AbsTol', 1e-6);
            testCase.verifyEqual(sims(1).blockFluxes, tExperiments.enzymeCap * tExperiments.byproductYield, 'AbsTol', 1e-6);
            % Ethanol: carbon-limited by the uptake of 20.
            testCase.verifyTrue(sims(2).feasible);
            testCase.verifyEqual(sims(2).growth, 20, 'AbsTol', 1e-6);
            testCase.verifyEqual(sims(2).exchFluxes(strcmp(sims(2).exchRxnIDs, 'EX_eth')), -20, 'AbsTol', 1e-6);
            testCase.verifyEqual(sims(2).blockFluxes, 20 * tExperiments.byproductYield, 'AbsTol', 1e-6);
            % The model is returned unchanged to the caller.
            testCase.verifyEqual(model.lb(strcmp(model.rxns, 'EX_glc')), -1000);
        end

        function simulateConditionsMaxGrowthOpensUptake(testCase)
            model = tExperiments.toyModel();
            sims = simulateConditions(model, tExperiments.maxGrowthToy(), false, {}, 'biomass');
            testCase.verifyEqual([sims.growth], [36 36], 'AbsTol', 1e-6);
            testCase.verifyEmpty(sims(1).exchRxnIDs);
            testCase.verifyEmpty(sims(1).blockRxnIDs);
        end

        function simulateConditionsUnmatchedCondition(testCase)
            model = tExperiments.toyModel();
            fd = tExperiments.fluxDataToy(20);
            fd.conds = {'glucose';'propionate'};
            testCase.verifyError(@() simulateConditions(model, fd, true, {}, 'biomass'), ...
                'simulateConditions:unmatchedCondition');
        end

        function simulateConditionsInfeasibleCondition(testCase)
            % A growth lower bound above what the enzyme pool allows.
            model = tExperiments.toyModel();
            model.lb(strcmp(model.rxns, 'biomass')) = 50;  % above the enzyme cap
            sims = simulateConditions(model, tExperiments.maxGrowthToy(), false, {}, 'biomass');
            testCase.verifyFalse(any([sims.feasible]));
        end

        function simulateConditionsUsesTheGivenOptimizeFunction(testCase)
            model = tExperiments.toyModel();
            calls = containers.Map({'n'}, {0});   % a handle, so the count survives the calls
            sims = simulateConditions(model, tExperiments.fluxDataToy(20), true, {}, 'biomass', ...
                'optimize', @(m) tExperiments.countedSolve(m, calls));
            testCase.verifyEqual(calls('n'), 2);
            testCase.verifyTrue(all([sims.feasible]));
        end

        function simulateConditionsAppliesTheHooks(testCase)
            % Oxygen measured at exactly 0 calls makeAnaerobic; a Ptot value
            % calls changeProteinBiomass with it.
            model = tExperiments.toyModel();
            fd = tExperiments.fluxDataToy(20);
            fd.exchMets = [fd.exchMets, {'oxygen'}];
            fd.exchRxnIDs = [fd.exchRxnIDs, {'EX_o2'}];
            fd.exchFluxes = [fd.exchFluxes, [0; NaN]];
            fd.Ptot = [NaN; 0.4];
            anaerobic = @(m) tExperiments.capReaction(m, 'R', 10);
            protein = @(m, ptot) tExperiments.capReaction(m, 'biomass', ptot*100);
            sims = simulateConditions(model, fd, true, {}, 'biomass', ...
                'makeAnaerobic', anaerobic, 'changeProteinBiomass', protein);
            testCase.verifyEqual(sims(1).growth, 10, 'AbsTol', 1e-6);   % capped by makeAnaerobic
            testCase.verifyEqual(sims(2).growth, 20, 'AbsTol', 1e-6);   % 40 cap above the carbon limit
            % Without hooks nothing is changed.
            sims = simulateConditions(model, fd, true, {}, 'biomass');
            testCase.verifyEqual(sims(1).growth, tExperiments.enzymeCap, 'AbsTol', 1e-6);
        end

        %% datasetRmse

        function datasetRmseOnlyTheByproductDeviates(testCase)
            model = tExperiments.toyModel();
            fd = tExperiments.fluxDataToy(20);
            sims = simulateConditions(model, fd, true, {'EX_byp'}, 'biomass');
            w = computeCarbonWeights(model, {'EX_glc','EX_eth','EX_byp','biomass'}, 'biomass');
            [rmse, perCond] = datasetRmse(fd, sims, true, w, 'biomass');
            % Growth and the own exchange match; the byproduct term is
            % measured 0 against the simulated value (carbon weight 1).
            expGlc = sqrt(mean([0 0 tExperiments.enzymeCap * tExperiments.byproductYield].^2));
            expEth = sqrt(mean([0 0 20 * tExperiments.byproductYield].^2));
            testCase.verifyEqual(perCond, [expGlc; expEth], 'AbsTol', 1e-6);
            testCase.verifyEqual(rmse, mean([expGlc expEth]), 'AbsTol', 1e-6);
        end

        function datasetRmseGrowthMismatch(testCase)
            model = tExperiments.toyModel();
            fd = tExperiments.fluxDataToy(18);
            sims = simulateConditions(model, fd, true, {'EX_byp'}, 'biomass');
            w = computeCarbonWeights(model, {'EX_glc','EX_eth','EX_byp','biomass'}, 'biomass');
            [~, perCond] = datasetRmse(fd, sims, true, w, 'biomass');
            growthDiff = (18 - 20) * tExperiments.biomassCarbon;
            bypDiff = 0 - 20 * tExperiments.byproductYield;
            testCase.verifyEqual(perCond(2), sqrt(mean([growthDiff 0 bypDiff].^2)), 'AbsTol', 1e-6);
        end

        function datasetRmseInfeasibleConditionGetsThePenalty(testCase)
            fd = tExperiments.fluxDataToy(20);
            sims = repmat(struct('feasible', false, 'growth', NaN, 'exchRxnIDs', {{}}, ...
                'exchFluxes', zeros(0,1), 'blockRxnIDs', {{}}, 'blockFluxes', zeros(0,1)), 2, 1);
            [rmse, perCond] = datasetRmse(fd, sims, true, containers.Map(), 'biomass', 'penalty', 99);
            testCase.verifyEqual(perCond, [99; 99]);
            testCase.verifyEqual(rmse, 99);
        end

        function datasetRmseMaxGrowthScoresGrowthOnly(testCase)
            fd = tExperiments.maxGrowthToy();
            sims = repmat(struct('feasible', true, 'growth', 0, 'exchRxnIDs', {{}}, ...
                'exchFluxes', zeros(0,1), 'blockRxnIDs', {{}}, 'blockFluxes', zeros(0,1)), 2, 1);
            sims(1).growth = tExperiments.enzymeCap;
            sims(2).growth = 30;
            [rmse, perCond] = datasetRmse(fd, sims, false, containers.Map(), 'biomass');
            second = abs((tExperiments.enzymeCap - 30) * tExperiments.biomassCarbon);
            testCase.verifyEqual(perCond, [0; second], 'AbsTol', 1e-9);
            testCase.verifyEqual(rmse, mean([0 second]), 'AbsTol', 1e-9);
        end

        function datasetRmseMultipliesByTheWeight(testCase)
            fd = tExperiments.maxGrowthToy();
            fd.weight = [0; 2];
            sims = repmat(struct('feasible', true, 'growth', 0, 'exchRxnIDs', {{}}, ...
                'exchFluxes', zeros(0,1), 'blockRxnIDs', {{}}, 'blockFluxes', zeros(0,1)), 2, 1);
            sims(1).growth = 30;
            sims(2).growth = 30;
            [rmse, perCond] = datasetRmse(fd, sims, false, containers.Map(), 'biomass');
            testCase.verifyEqual(rmse, mean(perCond .* [0; 2]), 'AbsTol', 1e-9);
        end

        function datasetRmseWrongSimCount(testCase)
            fd = tExperiments.maxGrowthToy();
            testCase.verifyError(@() datasetRmse(fd, struct('feasible', {}), false, ...
                containers.Map(), 'biomass'), 'datasetRmse:wrongSimCount');
        end

        %% experimentDistance

        function experimentDistanceCombinesBothDatasets(testCase)
            [ed, fluxSims, maxSims, w, expFlux, expMax] = tExperiments.distanceInputs();
            [rmse, detail] = experimentDistance(ed, fluxSims, maxSims, w, 'biomass');
            testCase.verifyEqual(rmse, mean([expFlux expMax]), 'AbsTol', 1e-9);
            testCase.verifyTrue(isfield(detail, 'fluxData') && isfield(detail, 'maxGrowth'));
        end

        function experimentDistanceMaxGrowthWeight(testCase)
            [ed, fluxSims, maxSims, w, expFlux, expMax] = tExperiments.distanceInputs();
            rmse = experimentDistance(ed, fluxSims, maxSims, w, 'biomass', 'maxGrowthWeight', 2);
            testCase.verifyEqual(rmse, (expFlux + 2*expMax) / 3, 'AbsTol', 1e-9);
            rmse = experimentDistance(ed, fluxSims, maxSims, w, 'biomass', 'maxGrowthWeight', 0.5);
            testCase.verifyEqual(rmse, (2*expFlux + expMax) / 3, 'AbsTol', 1e-9);
            % A lone dataset carries the whole score whatever the weight.
            solo = ed; solo.maxGrowth = [];
            rmse = experimentDistance(solo, fluxSims, [], w, 'biomass', 'maxGrowthWeight', 2);
            testCase.verifyEqual(rmse, expFlux, 'AbsTol', 1e-9);
        end

        function experimentDistanceMissingDatasetsContributeNothing(testCase)
            [ed, fluxSims, ~, w, expFlux] = tExperiments.distanceInputs();
            ed.maxGrowth = [];
            [rmse, detail] = experimentDistance(ed, fluxSims, [], w, 'biomass');
            testCase.verifyEqual(rmse, expFlux, 'AbsTol', 1e-9);
            testCase.verifyEqual(fieldnames(detail), {'fluxData'});
            none = struct('fluxData', [], 'maxGrowth', [], 'zeroFlux', {{}});
            testCase.verifyTrue(isnan(experimentDistance(none, [], [], w, 'biomass')));
        end

        function experimentDistanceNeedsTheSimulations(testCase)
            [ed, ~, maxSims, w] = tExperiments.distanceInputs();
            testCase.verifyError(@() experimentDistance(ed, [], maxSims, w, 'biomass'), ...
                'experimentDistance:missingSims');
        end

        %% applyFluxDataConstraints

        function applyFluxDataConstraintsRejectsInvalidSettings(testCase)
            m = tExperiments.simpleModel();
            fd = tExperiments.simpleFluxData({'c1'}, 0.4, -5, {'EX_glc'});
            testCase.verifyError(@() applyFluxDataConstraints(m, fd, 'biomass', 'maxMinGrowth', 'invalid'), ...
                'applyFluxDataConstraints:invalidMaxMinGrowth');
            testCase.verifyError(@() applyFluxDataConstraints(m, fd, 'biomass', 'looseStrictFlux', 200), ...
                'applyFluxDataConstraints:invalidLooseStrictFlux');
            testCase.verifyError(@() applyFluxDataConstraints(m, fd, 'biomass', 'condition', 'nonexistent'), ...
                'applyFluxDataConstraints:unknownCondition');
            testCase.verifyError(@() applyFluxDataConstraints(m, fd, 'biomass', 'condition', 5), ...
                'applyFluxDataConstraints:conditionOutOfRange');
            fd.exchRxnIDs = {'EX_nonexistent'};
            testCase.verifyError(@() applyFluxDataConstraints(m, fd, 'biomass'), ...
                'applyFluxDataConstraints:missingReaction');
        end

        function applyFluxDataConstraintsGrowthBounds(testCase)
            m = tExperiments.simpleModel();
            fd = tExperiments.simpleFluxData({'c1'}, 0.4, -5, {'EX_glc'});
            bio = strcmp(m.rxns, 'biomass');
            m2 = applyFluxDataConstraints(m, fd, 'biomass', 'maxMinGrowth', 'max');
            testCase.verifyEqual([m2.lb(bio) m2.ub(bio)], [0 0.4]);
            m2 = applyFluxDataConstraints(m, fd, 'biomass', 'maxMinGrowth', 'min');
            testCase.verifyEqual([m2.lb(bio) m2.ub(bio)], [0.4 1000]);
        end

        function applyFluxDataConstraintsZeroesTheCarbonSource(testCase)
            m = tExperiments.simpleModel();
            fd = tExperiments.simpleFluxData({'c1'}, 0.4, -5, {'EX_glc'});
            eth = strcmp(m.rxns, 'EX_eth');
            m2 = applyFluxDataConstraints(m, fd, 'biomass', 'cSource', 'EX_eth');
            testCase.verifyEqual([m2.lb(eth) m2.ub(eth)], [0 0]);
            % A carbon source that is not in the model is ignored, and
            % without cSource nothing is blocked.
            m2 = applyFluxDataConstraints(m, fd, 'biomass', 'cSource', 'EX_absent');
            testCase.verifyEqual(m2.lb(strcmp(m.rxns, 'EX_glc')), -5);
            m2 = applyFluxDataConstraints(m, fd, 'biomass');
            testCase.verifyEqual([m2.lb(eth) m2.ub(eth)], [-1000 1000]);
        end

        function applyFluxDataConstraintsLooseMode(testCase)
            m = tExperiments.simpleModel();
            fd = tExperiments.simpleFluxData({'c1'}, 0.4, [-5 3 NaN], {'EX_glc','EX_eth','EX_O2'});
            m2 = applyFluxDataConstraints(m, fd, 'biomass');
            testCase.verifyEqual(tExperiments.bounds(m2, 'EX_glc'), [-5 0]);
            testCase.verifyEqual(tExperiments.bounds(m2, 'EX_eth'), [0 3]);
            % NaN is no measurement: bounds are untouched.
            testCase.verifyEqual(tExperiments.bounds(m2, 'EX_O2'), tExperiments.bounds(m, 'EX_O2'));
        end

        function applyFluxDataConstraintsPercentageMode(testCase)
            m = tExperiments.simpleModel();
            fd = tExperiments.simpleFluxData({'c1'}, 0.4, [3 -5], {'EX_eth','EX_glc'});
            m2 = applyFluxDataConstraints(m, fd, 'biomass', 'looseStrictFlux', 10);
            testCase.verifyEqual(tExperiments.bounds(m2, 'EX_eth'), [2.85 3.15], 'AbsTol', 1e-12);
            % A negative value: the bounds are swapped to keep lb <= ub.
            testCase.verifyEqual(tExperiments.bounds(m2, 'EX_glc'), [-5.25 -4.75], 'AbsTol', 1e-12);
        end

        function applyFluxDataConstraintsExtremeValuesAreUnconstrained(testCase)
            m = tExperiments.simpleModel();
            fd = tExperiments.simpleFluxData({'c1'}, 0.4, [-1000 1000], {'EX_O2','EX_eth'});
            m2 = applyFluxDataConstraints(m, fd, 'biomass', 'looseStrictFlux', 10);
            testCase.verifyEqual(tExperiments.bounds(m2, 'EX_O2'), [-1000 0]);
            testCase.verifyEqual(tExperiments.bounds(m2, 'EX_eth'), [0 1000]);
        end

        function applyFluxDataConstraintsSelectsTheCondition(testCase)
            m = tExperiments.simpleModel();
            fd = tExperiments.simpleFluxData({'c1';'c2'}, [0.4;0.7], [-5;-8], {'EX_glc'});
            bio = strcmp(m.rxns, 'biomass');
            m2 = applyFluxDataConstraints(m, fd, 'biomass', 'condition', 2);
            testCase.verifyEqual(m2.ub(bio), 0.7);
            testCase.verifyEqual(m2.lb(strcmp(m.rxns, 'EX_glc')), -8);
            m2 = applyFluxDataConstraints(m, fd, 'biomass', 'condition', 'c2');
            testCase.verifyEqual(m2.lb(strcmp(m.rxns, 'EX_glc')), -8);
        end

    end

    methods (Static, Access = private)

        function file = writeFile(testCase, content)
            file = [tempname '.tsv'];
            fid = fopen(file, 'w');
            fwrite(fid, content);
            fclose(fid);
            testCase.addTeardown(@() delete(file));
        end

        function file = writeRows(testCase, rows)
            lines = cellfun(@(r) strjoin(r, sprintf('\t')), rows, 'UniformOutput', false);
            file = tExperiments.writeFile(testCase, [strjoin(lines, newline) newline]);
        end

        function writeNamed(dirName, fileName, content)
            fid = fopen(fullfile(dirName, fileName), 'w');
            fwrite(fid, content);
            fclose(fid);
        end

        function b = bounds(model, rxn)
            i = strcmp(model.rxns, rxn);
            b = [model.lb(i) model.ub(i)];
        end

        function sol = countedSolve(m, calls)
            calls('n') = calls('n') + 1;
            sol = solveLP(m);
        end

        function m = capReaction(m, rxn, value)
            m.ub(strcmp(m.rxns, rxn)) = value;
        end

        function fd = fluxDataToy(ethanolGrRate)
            % Glucose is enzyme-limited with matching measured uptake,
            % ethanol carbon-limited with an uptake of 20.
            fd.conds = {'glucose';'ethanol'};
            fd.Ptot = [NaN; NaN];
            fd.grRate = [tExperiments.enzymeCap; ethanolGrRate];
            fd.exchFluxes = [-tExperiments.enzymeCap NaN; NaN -20];
            fd.exchMets = {'glucose','ethanol'};
            fd.exchRxnIDs = {'EX_glc','EX_eth'};
            fd.weight = [1; 1];
            fd.source = {'test';'test'};
        end

        function fd = maxGrowthToy()
            fd.conds = {'glucose';'ethanol'};
            fd.Ptot = [NaN; NaN];
            fd.grRate = [tExperiments.enzymeCap; tExperiments.enzymeCap];
            fd.exchFluxes = [-1000 NaN; NaN -1000];
            fd.exchMets = {'glucose','ethanol'};
            fd.exchRxnIDs = {'EX_glc','EX_eth'};
        end

        function [ed, fluxSims, maxSims, w, expFlux, expMax] = distanceInputs()
            model = tExperiments.toyModel();
            fd = tExperiments.fluxDataToy(20);
            mg = tExperiments.maxGrowthToy();
            ed = struct('fluxData', fd, 'maxGrowth', mg, 'zeroFlux', {{'EX_byp'}});
            fluxSims = simulateConditions(model, fd, true, {'EX_byp'}, 'biomass');
            maxSims = simulateConditions(model, mg, false, {}, 'biomass');
            w = computeCarbonWeights(model, {'EX_glc','EX_eth','EX_byp','biomass'}, 'biomass');
            expFlux = datasetRmse(fd, fluxSims, true, w, 'biomass');
            expMax = datasetRmse(mg, maxSims, false, w, 'biomass');
        end

        function m = toyModel()
            rxns = {'EX_glc';'EX_eth';'EX_byp';'EX_o2';'TR_glc';'TR_eth';'R'; ...
                'biomass';'prot_pool_exchange';'usage_prot_E'};
            mets = {'glc_e';'eth_e';'byp_c';'o2_e';'A_c';'pool';'prot_E';'bio_met'};
            coeff = 100 / 3600;   % MW / (kcat * 3600), mg per mmol
            S = zeros(numel(mets), numel(rxns));
            S(1,1) = -1;                                 % EX_glc
            S(2,2) = -1;                                 % EX_eth
            S(3,3) = -1;                                 % EX_byp
            S(4,4) = -1;                                 % EX_o2
            S(1,5) = -1; S(5,5) = 1;                     % TR_glc
            S(2,6) = -1; S(5,6) = 1;                     % TR_eth
            S(5,7) = -1; S(7,7) = -coeff; S(8,7) = 1; S(3,7) = tExperiments.byproductYield;  % R
            S(8,8) = -1;                                 % biomass
            S(6,9) = 1;                                  % protein pool
            S(6,10) = -1; S(7,10) = 1;                   % usage of prot_E
            m = struct();
            m.id = 'toy';
            m.rxns = rxns;
            m.rxnNames = rxns;
            m.mets = mets;
            m.metNames = mets;
            m.metComps = ones(numel(mets), 1);
            m.comps = {'c'};
            m.compNames = {'cytosol'};
            m.S = sparse(S);
            m.lb = [-1000;-1000;0;-1000;0;0;0;0;0;0];
            m.ub = [0;0;1000;0;1000;1000;1000;1000;1;1000];
            m.rev = double(m.lb < 0);
            m.c = zeros(numel(rxns), 1); m.c(8) = 1;
            m.b = zeros(numel(mets), 1);
            m.genes = {}; m.grRules = repmat({''}, numel(rxns), 1);
            m.rxnGeneMat = sparse(numel(rxns), 0);
            m.metFormulas = {'C6H12O6';'C2H6O';'C1';'O2';'';'';'';''};
        end

        function m = simpleModel()
            m = struct();
            m.id = 'simple';
            m.rxns = {'EX_glc';'EX_O2';'EX_eth';'biomass'};
            m.rxnNames = m.rxns;
            m.mets = {'glc_e';'O2_e';'eth_e';'bio'};
            m.metNames = m.mets;
            m.metComps = ones(4, 1);
            m.comps = {'c'};
            m.compNames = {'cytosol'};
            m.S = sparse([-1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]);
            m.lb = [-1000;-1000;-1000;0];
            m.ub = [1000;1000;1000;1000];
            m.rev = [1;1;1;0];
            m.c = [0;0;0;1];
            m.b = zeros(4, 1);
            m.genes = {}; m.grRules = repmat({''}, 4, 1);
            m.rxnGeneMat = sparse(4, 0);
        end

        function fd = simpleFluxData(conds, grRate, fluxes, rxnIds)
            fd.conds = conds(:);
            fd.Ptot = NaN(numel(conds), 1);
            fd.grRate = grRate(:);
            fd.exchFluxes = reshape(fluxes, numel(conds), []);
            fd.exchMets = rxnIds;
            fd.exchRxnIDs = rxnIds;
        end

    end
end
