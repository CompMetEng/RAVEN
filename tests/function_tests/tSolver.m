classdef tSolver < RavenTestCase
% tSolver  Tests for the LP/QP solving and solver-configuration functions in solver/.

    methods (Test)

        function setRavenSolverAcceptsGlpk(testCase)
            testCase.verifyWarningFree(@() setRavenSolver('glpk'));
        end

        function setRavenSolverRejectsUnknown(testCase)
            testCase.verifyError(@() setRavenSolver('notarealsolver'), ?MException);
        end

        function solveLPFindsOptimum(testCase)
            sol = solveLP(testCase.model);
            testCase.verifyEqual(sol.stat, 1);
            testCase.verifyNumElements(sol.x, numel(testCase.model.rxns));
            testCase.verifyEqual(sol.f, 0.87392, 'AbsTol', 1e-3);
        end

        function solveLPMinimiseFluxes(testCase)
            sol = solveLP(testCase.model, 1);
            testCase.verifyEqual(sol.stat, 1);
        end

        function optimizeProbSolvesLP(testCase)
            res = optimizeProb(testCase.modelLP());
            testCase.verifyClass(res, 'struct');
            testCase.verifyTrue(isfield(res, 'obj'));
        end

        function optimizeProbReportsTimeLimit(testCase)
            res = optimizeProb(testCase.modelLP());
            testCase.verifyTrue(isfield(res, 'hitTimeLimit'));
            testCase.verifyTrue(islogical(res.hitTimeLimit));
            %The problem is solved long before the time limit is reached
            testCase.verifyFalse(res.hitTimeLimit);
        end

        function optimizeProbSolvesMILPWithGlpk(testCase)
            % A binary knapsack whose LP relaxation (9.33) differs from the
            % integer optimum (8, items 1 and 3). GLPK solves it, with a
            % warning that recommends Gurobi.
            evalc('setRavenSolver(''glpk'')');
            prob = struct();
            prob.c       = [-5; -4; -3];
            prob.A       = sparse([2 3 1]);
            prob.a       = prob.A;
            prob.b       = 4;
            prob.csense  = 'L';
            prob.lb      = zeros(3,1);
            prob.ub      = ones(3,1);
            prob.osense  = 1;
            prob.vartype = 'BBB';
            res = testCase.verifyWarning(@() optimizeProb(prob,[],false), 'RAVEN:glpkMILP');
            testCase.verifyEqual(res.stat, 1);
            testCase.verifyEqual(res.obj, -8, 'AbsTol', 1e-9);
            testCase.verifyEqual(res.full, [1; 0; 1], 'AbsTol', 1e-9);
            testCase.verifyEqual(res.mipgap, 0);
        end

        function optimizeProbSolvesMILPWithCobraGlpk(testCase)
            % The same knapsack through the COBRA Toolbox's GLPK interface.
            testCase.assumeDependency(~isempty(which('changeCobraSolver')), 'COBRA Toolbox');
            global CBT_MILP_SOLVER %#ok<GVMIS>
            origMILP = CBT_MILP_SOLVER;
            evalc('changeCobraSolver(''glpk'',''MILP'',0);');
            evalc('setRavenSolver(''cobra'')');
            testCase.addTeardown(@() evalc('setRavenSolver(''glpk'')'));
            if ~isempty(origMILP)
                testCase.addTeardown(@() changeCobraSolver(origMILP,'MILP',0));
            end
            prob = struct();
            prob.c       = [-5; -4; -3];
            prob.A       = sparse([2 3 1]);
            prob.a       = prob.A;
            prob.b       = 4;
            prob.csense  = 'L';
            prob.lb      = zeros(3,1);
            prob.ub      = ones(3,1);
            prob.osense  = 1;
            prob.vartype = 'BBB';
            res = testCase.verifyWarning(@() optimizeProb(prob,[],false), 'RAVEN:glpkMILP');
            testCase.verifyEqual(res.stat, 1);
            testCase.verifyEqual(res.obj, -8, 'AbsTol', 1e-6);
            testCase.verifyEqual(res.full, [1; 0; 1], 'AbsTol', 1e-6);
        end

        function getMinimalMediumRunsOnGlpk(testCase)
            evalc('setRavenSolver(''glpk'')');
            warnState = warning('off','RAVEN:glpkMILP');
            testCase.addTeardown(@() warning(warnState));
            [medium, idx] = getMinimalMedium(testCase.model, 'verbose', false);
            testCase.verifyNotEmpty(medium);
            m = testCase.model;
            inMedium = false(numel(m.rxns),1); inMedium(idx) = true;
            [~, exchIdx] = getExchangeRxns(m);
            m.lb(intersect(exchIdx, find(~inMedium & m.lb < 0))) = 0;
            sol = solveLP(m);
            testCase.verifyGreaterThan(sol.f, 0);
        end

        function checkSolutionOnFeasibleProblem(testCase)
            res = optimizeProb(testCase.modelLP());
            [isFeasible, isOptimal] = checkSolution(res);
            testCase.verifyTrue(logical(isFeasible));
            testCase.verifyTrue(logical(isOptimal));
        end

        function splitProbForConditioningReturnsProb(testCase)
            [p2, condInfo] = splitProbForConditioning(testCase.modelLP(), 1e6); %#ok<ASGLU>
            testCase.verifyClass(p2, 'struct');
        end

    end

    methods (Access = private)
        function prob = modelLP(testCase)
            % Build the COBRA-style LP problem for the test model exactly as
            % solveLP does, so the low-level solver functions get a valid input.
            m = testCase.model;
            nMets = size(m.S, 1);
            prob = [];
            prob.c = [m.c * -1; zeros(nMets, 1)];
            prob.b = zeros(nMets, 1);
            prob.lb = [m.lb; m.b(:, 1)];
            prob.ub = [m.ub; m.b(:, min(size(m.b, 2), 2))];
            prob.csense = repmat('E', 1, numel(prob.b));
            prob.osense = 1;
            prob.vartype = repmat('C', 1, numel(prob.c));
            prob.A = [m.S, -speye(nMets)];
            prob.a = m.S;
        end
    end
end
