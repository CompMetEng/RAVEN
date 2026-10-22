classdef tBiomass < RavenTestCase
% tBiomass  Tests for the biomass-composition functions in biomass/.
%
%   The textbook E. coli model has a single lumped biomass reaction rather
%   than per-component pseudoreactions, so component-scaling functions are
%   tested on their missing-pseudoreaction error path.

    methods (Test)

        function getBiomassFractionsHandlesMissingComponents(testCase)
            cfg = testCase.biomassConfig();
            fractions = getBiomassFractions(testCase.model, cfg);
            testCase.verifyTrue(isfield(fractions, 'total'));
        end

        function setGAMReturnsModel(testCase)
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            evalc(['m2 = setGAM(testCase.model, 50, biomassRxn, ' ...
                '{''ATP'',''ADP'',''H2O'',''H+'',''Phosphate''});']);
            testCase.verifyClass(m2, 'struct');
        end

        function scaleBiomassFractionErrorsOnMissingComponent(testCase)
            cfg = testCase.biomassConfig();
            testCase.verifyError(@() scaleBiomassFraction(testCase.model, cfg, ...
                'protein', 0.5), ?MException);
        end

        function scaleBiomassFractionBalancesToOne(testCase)
            [m, cfg] = tBiomass.twoComponentToyModel();
            out = scaleBiomassFraction(m, cfg, 'protein', 0.7, 'balanceOut', 'DNA');
            f = getBiomassFractions(out, cfg);
            testCase.verifyEqual(f.protein, 0.7, 'AbsTol', 1e-9);
            testCase.verifyEqual(f.DNA, 0.3, 'AbsTol', 1e-9);
        end

        function scaleBiomassFractionBalancesToTargetTotal(testCase)
            [m, cfg] = tBiomass.twoComponentToyModel();
            out = scaleBiomassFraction(m, cfg, 'protein', 0.7, 'balanceOut', 'DNA', ...
                'targetTotal', 0.9);
            f = getBiomassFractions(out, cfg);
            testCase.verifyEqual(f.protein, 0.7, 'AbsTol', 1e-9);
            testCase.verifyEqual(f.DNA, 0.2, 'AbsTol', 1e-9);
            testCase.verifyEqual(f.total, 0.9, 'AbsTol', 1e-9);
        end

        function scaleBiomassFractionRefusesNegativeBalanceAtTargetTotal(testCase)
            % Protein 0.85 balances against DNA for a total of 1, but not
            % for a total of 0.8.
            [m, cfg] = tBiomass.twoComponentToyModel();
            testCase.verifyError(@() scaleBiomassFraction(m, cfg, 'protein', 0.85, ...
                'balanceOut', 'DNA', 'targetTotal', 0.8), 'scaleBiomassFraction:negativeBalance');
        end

        function scaleBiomassFractionRefusesNegativeBalance(testCase)
            % Balancing protein 1.2 against DNA would need DNA at -0.2,
            % which would turn the DNA pseudoreaction around.
            [m, cfg] = tBiomass.twoComponentToyModel();
            testCase.verifyError(@() scaleBiomassFraction(m, cfg, 'protein', 1.2, ...
                'balanceOut', 'DNA'), 'scaleBiomassFraction:negativeBalance');
        end

        function scaleBiomassPseudoreactionErrorsOnMissingComponent(testCase)
            cfg = testCase.biomassConfig();
            testCase.verifyError(@() scaleBiomassPseudoreaction(testCase.model, cfg, ...
                'protein', 0.9), ?MException);
        end

        function scaleBiomassPseudoreactionRebalancesProton(testCase)
            % Rescaling the substrates rebalances H+ so the pseudoreaction
            % stays charge neutral.
            [m, cfg] = tBiomass.protonToyModel();
            out = scaleBiomassPseudoreaction(m, cfg, 'protein', 0.5);
            % x_c: -0.5 * charge -2 = +1, so H+ must come in at -1
            testCase.verifyEqual(full(out.S(strcmp(out.mets,'x_c'), 1)), -0.5, 'AbsTol', 1e-9);
            testCase.verifyEqual(full(out.S(strcmp(out.mets,'h_c'), 1)), -1, 'AbsTol', 1e-9);
        end

        function scaleBiomassPseudoreactionRefusesMissingProton(testCase)
            [m, cfg] = tBiomass.protonToyModel();
            cfg.proton_met = 'nohydrogen_c';
            testCase.verifyError(@() scaleBiomassPseudoreaction(m, cfg, 'protein', 0.5), ...
                'scaleBiomassPseudoreaction:missingProton');
        end

        function setGAMRefusesMissingNgamRxn(testCase)
            % The model is not changed when the NGAM reaction is not found.
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            testCase.verifyError(@() setGAM(testCase.model, 50, biomassRxn, ...
                {'ATP','ADP','H2O','H+','Phosphate'}, 'ngamRxn', 'noSuchRxn', 'ngamValue', 1), ...
                'setGAM:missingNgamRxn');
        end

        function setGAMRefusesMissingNgamRxnWithoutValue(testCase)
            % The reaction is checked whenever it is given.
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            testCase.verifyError(@() setGAM(testCase.model, 50, biomassRxn, ...
                {'ATP','ADP','H2O','H+','Phosphate'}, 'ngamRxn', 'noSuchRxn'), ...
                'setGAM:missingNgamRxn');
        end

        function scaleBiomassPseudoreactionRefusesUnknownCharge(testCase)
            % An unset charge on a participating metabolite leaves the charge
            % balance unknown. Treating it as neutral would silently write the
            % wrong proton coefficient, so the rescale must refuse instead.
            [m, cfg] = tBiomass.protonToyModel();
            m.metCharges(strcmp(m.mets,'x_c')) = NaN;
            testCase.verifyError(@() scaleBiomassPseudoreaction(m, cfg, 'protein', 0.5), ...
                'scaleBiomassPseudoreaction:unknownCharge');
        end

        function fitParametersSetsUpperBoundPerXRxn(testCase)
            % model.ub for a second xRxn must come from its own column of
            % xValues, not the first column broadcast to every xRxn (which
            % also inverts that reaction's bounds whenever its value is
            % smaller than the first column's).
            m = struct();
            m.id='t'; m.rxns={'R1';'R2';'R3';'R4'}; m.rxnNames=m.rxns;
            m.mets={'m1';'m2'}; m.metNames=m.mets; m.metComps=[1;1];
            m.comps={'c'}; m.compNames={'c'};
            m.S=sparse([1 0 -1 0; 0 1 0 -1]); % R1:=>m1  R2:=>m2  R3:m1=>  R4:m2=>
            m.lb=[0;0;0;0]; m.ub=[1000;1000;1000;1000]; m.rev=[0;0;0;0];
            m.c=[0;0;0;0]; m.b=zeros(2,1);
            m.genes={}; m.grRules={'';'';'';''}; m.rxnGeneMat=sparse(4,0);
            pos.position={[]}; pos.isNegative={[]};
            evalc(['[~,fitnessScore] = fitParameters(m, {''R1'';''R2''}, [5 7], ' ...
                '{''R3'';''R4''}, [5 7], pos, ''fitToRatio'', false);']);
            % R3/R4 are only reachable via mass balance from R1/R2's fixed
            % (lb=ub) values, so a correct fit has exactly zero residual.
            testCase.verifyEqual(fitnessScore, 0, 'AbsTol', 1e-9);
        end

        function fitParametersRuns(testCase)
            % Minimal single-parameter fit of an exchange flux. fitParameters
            % uses fminsearch (base MATLAB), not quadprog, so no toolbox guard.
            biomassRxnIdx = find(testCase.model.c == 1, 1);
            biomassRxn = testCase.model.rxns{biomassRxnIdx};
            % parameterPositions.position/isNegative are cell arrays, one cell
            % per fitted parameter, of linear S-matrix indices (see
            % fitParameters' docstring) -- fit the coefficient of the first
            % reactant in the biomass reaction.
            metIdx = find(testCase.model.S(:,biomassRxnIdx) < 0, 1);
            pos.position = {sub2ind(size(testCase.model.S), metIdx, biomassRxnIdx)};
            pos.isNegative = {true};
            evalc(['p = fitParameters(testCase.model, {biomassRxn}, 0.5, ' ...
                '{biomassRxn}, 0.5, pos);']);
            testCase.verifyNotEmpty(p);
        end

    end

    methods (Static, Access = private)
        function [m, cfg] = twoComponentToyModel()
            % Protein and DNA pseudoreactions in grams, 0.5 g/gDW each.
            m = struct();
            m.id        = 'toy';
            m.rxns      = {'R1';'R2'};
            m.rxnNames  = {'protein pseudoreaction';'DNA pseudoreaction'};
            m.mets      = {'x_c';'y_c';'h_c';'prot_c';'dna_c'};
            m.metNames  = {'x';'y';'H+';'protein';'DNA'};
            m.metComps  = ones(5,1);
            m.metCharges= zeros(5,1);
            m.comps     = {'c'};
            m.compNames = {'cytosol'};
            m.S         = sparse([-0.5 0; 0 -0.5; 0 0; 1 0; 0 1]);
            m.lb = [0;0]; m.ub = [1000;1000]; m.rev = [0;0]; m.c = [0;0]; m.b = zeros(5,1);
            m.grRules = {'';''}; m.genes = {}; m.rxnGeneMat = sparse(2,0);

            cfg = struct();
            cfg.biomass_rxn = 'R1';
            cfg.proton_met  = 'h_c';
            cfg.components  = {struct('name','protein', ...
                'pseudoreaction_name','protein pseudoreaction','mass_strategy','grams'), ...
                struct('name','DNA', ...
                'pseudoreaction_name','DNA pseudoreaction','mass_strategy','grams')};
        end

        function [m, cfg] = protonToyModel()
            % One pseudoreaction: x_c -> protein, with h_c available to carry
            % the charge balance.
            m = struct();
            m.id        = 'toy';
            m.rxns      = {'R1'};
            m.rxnNames  = {'protein pseudoreaction'};
            m.mets      = {'x_c';'h_c';'prot_c'};
            m.metNames  = {'x';'h';'protein'};
            m.metComps  = [1;1;1];
            m.metCharges= [-2;1;0];
            m.comps     = {'c'};
            m.compNames = {'cytosol'};
            m.S         = sparse([-1;0;1]);
            m.lb = 0; m.ub = 1000; m.rev = 0; m.c = 1; m.b = zeros(3,1);
            m.grRules = {''}; m.genes = {}; m.rxnGeneMat = sparse(1,0);

            cfg = struct();
            cfg.biomass_rxn = 'R1';
            cfg.proton_met  = 'h_c';
            cfg.components  = {struct('name','protein', ...
                'pseudoreaction_name','protein pseudoreaction', ...
                'mass_strategy','mw')};
        end
    end

    methods (Access = private)
        function cfg = biomassConfig(testCase)
            biomassRxn = testCase.model.rxns{find(testCase.model.c == 1, 1)};
            cfg.biomass_rxn = biomassRxn;
            cfg.proton_met = 'h_c';
            cfg.components = {struct('name','protein', ...
                'pseudoreaction_name','no_such_pseudoreaction','mass_strategy','mw')};
        end
    end
end
