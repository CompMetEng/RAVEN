classdef tSampling < RavenTestCase
% tSampling  Tests for the MCMC flux samplers (CHRR, ACHR, MVE) dispatched
% through randomSampling.
%
%   The maximum-volume ellipsoid solver is the numerical crux of CHRR, so it is
%   validated directly against analytic cases (box -> unit ball, scaled/sheared
%   box, triangle -> Steiner inellipse) before the end-to-end samplers are
%   exercised on the textbook model. All methods are LP-only, so the default
%   GLPK solver suffices (no MILP assumption needed).

    methods (Test)

        % ---- sampleMaxVolEllipse: analytic validation -------------------

        function sampleMaxVolEllipseUnitBox(testCase)
            % Box {-1 <= x_i <= 1} -> unit ball: centre 0, E*E' = I.
            A = [eye(3); -eye(3)];
            b = ones(6, 1);
            [x, E, converged] = sampleMaxVolEllipse(A, b, zeros(3, 1));
            testCase.verifyTrue(converged);
            testCase.verifyLessThan(norm(x), 1e-4);
            testCase.verifyLessThan(norm(E*E' - eye(3)), 1e-3);
        end

        function sampleMaxVolEllipseScaledBox(testCase)
            % Box {-r_i <= x_i <= r_i} -> E*E' = diag(r_i^2).
            r = [2; 0.5; 3];
            A = [eye(3); -eye(3)];
            b = [r; r];
            [x, E, converged] = sampleMaxVolEllipse(A, b, zeros(3, 1));
            testCase.verifyTrue(converged);
            testCase.verifyLessThan(norm(x), 1e-4);
            testCase.verifyLessThan(norm(E*E' - diag(r.^2)), 1e-3);
        end

        function sampleMaxVolEllipseShearedBox(testCase)
            % {-1 <= M*x <= 1} -> E*E' = inv(M)*inv(M)'.
            M = [1 0.5; 0 1];
            A = [M; -M];
            b = ones(4, 1);
            [x, E, converged] = sampleMaxVolEllipse(A, b, zeros(2, 1));
            Minv = inv(M);
            testCase.verifyTrue(converged);
            testCase.verifyLessThan(norm(x), 1e-4);
            testCase.verifyLessThan(norm(E*E' - Minv*Minv'), 1e-3);
        end

        function sampleMaxVolEllipseTriangle(testCase)
            % Standard simplex {x>=0, x1+x2<=1}: the max-area inscribed ellipse is
            % the Steiner inellipse, centred at the centroid (1/3,1/3) with
            % E2 = [[1/9,-1/18],[-1/18,1/9]] (det 1/108) -- NOT the incircle. This
            % confirms the solver maximises volume, not just the largest ball.
            A = [-1 0; 0 -1; 1 1];
            b = [0; 0; 1];
            [x, E, converged] = sampleMaxVolEllipse(A, b, [0.25; 0.25]);
            E2exp = [1/9 -1/18; -1/18 1/9];
            testCase.verifyTrue(converged);
            testCase.verifyLessThan(norm(x - [1/3; 1/3]), 1e-3);
            testCase.verifyLessThan(norm(E*E' - E2exp), 1e-3);
            testCase.verifyLessThan(abs(det(E*E') - 1/108), 1e-4);
        end

        function sampleMaxVolEllipseRejectsExterior(testCase)
            A = [eye(2); -eye(2)];
            b = ones(4, 1);
            testCase.verifyError(@() sampleMaxVolEllipse(A, b, [5; 0]), ...
                'RAVEN:sampling');
        end

        % ---- sampler implementations return diagnostics -----------------

        function sampleCHRRReturnsInfo(testCase)
            evalc('[sols, sInfo] = sampleCHRR(testCase.model, 10, 5, 20, 1);');
            nRxns = numel(testCase.model.rxns);
            testCase.verifySize(sols, [nRxns, 10]);
            testCase.verifyTrue(isfield(sInfo, 'nDimensions'));
            testCase.verifyTrue(isfield(sInfo, 'mveConverged'));
            % Samples must satisfy steady state and bounds.
            resid = testCase.model.S * sols;
            testCase.verifyLessThan(max(abs(resid(:))), 1e-6);
            testCase.verifyGreaterThanOrEqual(min(sols, [], 2), testCase.model.lb - 1e-6);
            testCase.verifyLessThanOrEqual(max(sols, [], 2), testCase.model.ub + 1e-6);
        end

        function sampleACHRReturnsSteadyState(testCase)
            lastwarn('');
            evalc('sols = sampleACHR(testCase.model, 10, 5, [], 1);');
            testCase.verifyNotEqual(tSampling.lastWarningId(), 'sampleACHR:steadyStateNotRestored');
            nRxns = numel(testCase.model.rxns);
            testCase.verifySize(sols, [nRxns, 10]);
            resid = testCase.model.S * sols;
            testCase.verifyLessThan(max(abs(resid(:))), 1e-6);
        end

        function sampleACHRKeepsSteadyStateOverManySteps(testCase)
            % The residual of S*v must stay at rounding level when the chain
            % runs for many steps: with 100 steps between samples, a step that
            % is not confined to the null space of S lets the residual grow
            % by orders of magnitude over 1000 samples.
            m = tSampling.parallelModel();
            evalc('warmup = sampleWarmupPoints(m, ''runParallel'', false);');
            lastwarn('');
            evalc('sols = sampleACHR(m, 1000, 100, warmup, 1);');
            testCase.verifyNotEqual(tSampling.lastWarningId(), 'sampleACHR:steadyStateNotRestored');
            resid = m.S * sols;
            testCase.verifyLessThan(max(abs(resid(:))), 1e-9);
            testCase.verifyGreaterThanOrEqual(min(sols, [], 2), m.lb - 1e-9);
            testCase.verifyLessThanOrEqual(max(sols, [], 2), m.ub + 1e-9);
            testCase.verifyGreaterThan(std(sols(2, :)), 0.1);
        end

        function sampleACHRMovesOnModelWithFixedReaction(testCase)
            % Fixed (and thereby forced) reactions keep an exact zero step; the chain
            % must still move at (nearly) every step, stay inside the bounds
            % and keep S*v at rounding level.
            m = tSampling.fixedReactionModel();
            evalc('warmup = sampleWarmupPoints(m, ''runParallel'', false);');
            lastwarn('');
            evalc('sols = sampleACHR(m, 3000, 1, warmup, 1);');
            testCase.verifyNotEqual(tSampling.lastWarningId(), 'sampleACHR:steadyStateNotRestored');
            moved = mean(any(diff(sols, 1, 2) ~= 0, 1));
            testCase.verifyGreaterThan(moved, 0.9);
            testCase.verifyLessThan(max(abs(m.S * sols), [], 'all'), 1e-9);
            testCase.verifyGreaterThanOrEqual(min(sols, [], 2), m.lb - 1e-9);
            testCase.verifyLessThanOrEqual(max(sols, [], 2), m.ub + 1e-9);
            testCase.verifyEqual(sols([3 5], :), 700*ones(2, 3000), 'AbsTol', 1e-9);
            evalc('sols = sampleACHR(m, 300, 100, warmup, 2);');
            testCase.verifyLessThan(max(abs(m.S * sols), [], 'all'), 1e-9);
            testCase.verifyGreaterThan(std(sols(2, :)), 0.1);
        end

        function sampleACHRKeepsNearlyConstantReactionConstant(testCase)
            % A reaction whose warmup values differ by less than 1e-9 does not move:
            % its recorded value is the warmup mean in every sample, and the
            % steady state is kept through the other reactions.
            m = tSampling.nearlyFixedReactionModel();
            evalc('warmup = sampleWarmupPoints(m, ''runParallel'', false);');
            s = rng; testCase.addTeardown(@() rng(s));
            rng(5);
            warmup(:, end) = 3 + 5e-10*rand(size(warmup, 1), 1);
            sols = sampleACHR(m, 200, 20, warmup, 1);
            testCase.verifyEqual(max(sols(end, :)) - min(sols(end, :)), 0);
            testCase.verifyLessThan(max(abs(m.S * sols), [], 'all'), 1e-9);
            testCase.verifyGreaterThan(std(sols(2, :)), 0.1);
        end

        function sampleACHRStaysAtSteadyStateOverMillionsOfSteps(testCase)
            % One sample after two million steps: the residual of S*v must not grow
            % between recorded samples, which needs the correction every few steps
            % (a correction only before the recorded sample leaves a residual of
            % about 0.3 here).
            m = tSampling.parallelModel();
            evalc('warmup = sampleWarmupPoints(m, ''runParallel'', false);');
            sols = sampleACHR(m, 1, 2e6, warmup, 1);
            testCase.verifyLessThan(max(abs(m.S * sols), [], 'all'), 1e-9);
            testCase.verifyGreaterThanOrEqual(sols, m.lb - 1e-9);
            testCase.verifyLessThanOrEqual(sols, m.ub + 1e-9);
        end

        function sampleACHRMovesFromWarmupWithRoundingNoise(testCase)
            % Warmup values that differ from a vertex by 1e-9 to 3e-9 must not stall
            % the chain: reactions with no warmup spread get an exact zero step,
            % directions below a tolerance relative to the largest component are
            % ignored, and the correction keeps the point inside the polytope. The
            % fractions of moving steps are 0.99 and 0.96; without one of these
            % three, they drop to 0.6 or below at one of the two noise levels.
            m = testCase.model;
            evalc('warmup = sampleWarmupPoints(m, ''runParallel'', false);');
            s = rng; testCase.addTeardown(@() rng(s));
            rng(3);
            noise = randn(size(warmup));
            sols = sampleACHR(m, 300, 20, warmup + 1e-9*noise, 1);
            testCase.verifyGreaterThan(mean(any(diff(sols, 1, 2) ~= 0, 1)), 0.97);
            sols = sampleACHR(m, 300, 20, warmup + 3e-9*noise, 1);
            testCase.verifyGreaterThan(mean(any(diff(sols, 1, 2) ~= 0, 1)), 0.93);
            testCase.verifyLessThan(max(abs(m.S * sols), [], 'all'), 1e-8);
        end

        function sampleACHRWarnsForWarmupThatIsInfeasibleForTheModel(testCase)
            % Warmup points of the model without the fixed ATP maintenance flux do
            % not satisfy the model with it; the samples keep a residual far above
            % rounding and sampleACHR says so.
            m = testCase.model;
            evalc('warmup = sampleWarmupPoints(m, ''runParallel'', false);');
            i = find(contains(lower(m.rxns), 'atpm'), 1);
            testCase.assertNotEmpty(i);
            m.lb(i) = 20;
            m.ub(i) = 20;
            testCase.verifyWarning(@() sampleACHR(m, 200, 20, warmup, 1), ...
                'sampleACHR:steadyStateNotRestored');
        end

        function sampleACHRWarnsForNoisyWarmup(testCase)
            % Warmup values that are off by 1e-5 leave a residual above 1e-9 after
            % the clip to the bounds.
            m = testCase.model;
            evalc('warmup = sampleWarmupPoints(m, ''runParallel'', false);');
            s = rng; testCase.addTeardown(@() rng(s));
            rng(3);
            noisy = warmup + 1e-5*randn(size(warmup));
            testCase.verifyWarning(@() sampleACHR(m, 300, 20, noisy, 1), ...
                'sampleACHR:steadyStateNotRestored');
        end

        function sampleACHRWarnsForNearlyDependentRows(testCase)
            % Two rows of S that differ by 1e-7 times a column: the correction
            % cannot separate them, the implied constraint is lost and the
            % residual stays near 3e-7.
            d = 1e-7;
            S = [1  -1      -1  0  0;
                 1  -(1+d)  -1  0  d;
                 0   1       1 -1  0;
                 0   0       0  1 -1];
            m = tSampling.modelFromMatrix(S, zeros(5, 1), [10; 6; 6; 1000; 1000]);
            evalc('warmup = sampleWarmupPoints(m, ''runParallel'', false);');
            testCase.verifyWarning(@() sampleACHR(m, 500, 100, warmup, 1), ...
                'sampleACHR:steadyStateNotRestored');
        end

        % ---- randomSampling dispatch ------------------------------------

        function randomSamplingDefaultsToACHR(testCase)
            % The default method must be ACHR: with a fixed seed the default call
            % must equal an explicit 'achr' call.
            evalc('sDefault = randomSampling(testCase.model, 8, ''thinning'', 5, ''seed'', 3);');
            evalc('sAchr = randomSampling(testCase.model, 8, ''method'', ''achr'', ''thinning'', 5, ''seed'', 3);');
            testCase.verifyEqual(sDefault, sAchr);
        end

        function randomSamplingCHRRDispatch(testCase)
            evalc(['sols = randomSampling(testCase.model, 10, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 1);']);
            nRxns = numel(testCase.model.rxns);
            testCase.verifySize(sols, [nRxns, 10]);
            resid = testCase.model.S * sols;
            testCase.verifyLessThan(max(abs(resid(:))), 1e-6);
        end

        function randomSamplingCHRRWithGurobi(testCase)
            % sampleChebyshevCenter must hand Gurobi a sparse constraint
            % matrix, which Gurobi requires.
            testCase.assumeSolver('gurobi');
            evalc('setRavenSolver(''gurobi'')');
            testCase.addTeardown(@() evalc('setRavenSolver(''glpk'')'));
            evalc(['sols = randomSampling(testCase.model, 10, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 1);']);
            testCase.verifySize(sols, [numel(testCase.model.rxns), 10]);
        end

        function randomSamplingCHRRExposesInfo(testCase)
            % The CHRR convergence diagnostics must be reachable through the
            % documented randomSampling entry point, not only sampleCHRR.
            evalc(['[sols, gr, sInfo] = randomSampling(testCase.model, 10, ' ...
                '''method'', ''chrr'', ''thinning'', 5, ''nBurnin'', 20, ''seed'', 1);']);
            testCase.verifyTrue(isfield(sInfo, 'mveConverged'));
        end

        function randomSamplingCHRRDeterministic(testCase)
            evalc(['s1 = randomSampling(testCase.model, 8, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 42);']);
            evalc(['s2 = randomSampling(testCase.model, 8, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 42);']);
            testCase.verifyEqual(s1, s2);
        end

        function randomSamplingCHRRVaries(testCase)
            % A correct CHRR run must mix: at least one reaction varies, and two
            % seeds give distinct sample matrices. Guards against a degenerate
            % sampler that collapses to a single point (which would still satisfy
            % S*v=0 and the bounds, and remain deterministic).
            evalc(['s1 = randomSampling(testCase.model, 50, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 1);']);
            evalc(['s2 = randomSampling(testCase.model, 50, ''method'', ''chrr'', ' ...
                '''thinning'', 5, ''nBurnin'', 20, ''seed'', 2);']);
            testCase.verifyGreaterThan(max(std(s1, 0, 2)), 1e-3);
            testCase.verifyNotEqual(s1, s2);
        end

        function randomSamplingUnknownMethodErrors(testCase)
            testCase.verifyError( ...
                @() randomSampling(testCase.model, 5, 'method', 'gibbs'), ...
                'RAVEN:badInput');
        end

        function randomSamplingLoopDetectionUsesModelBounds(testCase)
            % Loop detection must take its threshold from the model. R2/R3
            % form an a -> b -> a loop that runs up to this model's 100 cap,
            % while the linear path R1 -> R4 is held to 10. A threshold that
            % assumes ±1000 bounds matches neither, leaving the loop
            % reactions to be sampled as if they were loop-free.
            m = tSampling.loopModel();
            evalc(['[~, goodRxns] = randomSampling(m, 2, ''method'', ' ...
                '''randomObjective'', ''seed'', 1);']);
            testCase.verifyEqual(sort(goodRxns(:))', [1 4]);
        end

    end

    methods (Static, Access = private)

        function m = modelFromMatrix(S, lb, ub)
            % A model with the given stoichiometric matrix and bounds, without genes.
            [nMets, nRxns] = size(S);
            m = struct();
            m.id       = 'fromMatrix';
            m.rxns     = arrayfun(@(i) sprintf('R%d', i), (1:nRxns)', 'UniformOutput', false);
            m.rxnNames = m.rxns;
            m.mets     = arrayfun(@(i) sprintf('m%d', i), (1:nMets)', 'UniformOutput', false);
            m.metNames = m.mets;
            m.metComps = ones(nMets, 1);
            m.comps    = {'c'};
            m.compNames= {'cytosol'};
            m.S = sparse(S);
            m.lb  = lb;
            m.ub  = ub;
            m.rev = zeros(nRxns, 1);
            m.c   = zeros(nRxns, 1);
            m.b   = zeros(nMets, 1);
            m.genes = {}; m.grRules = repmat({''}, nRxns, 1); m.rxnGeneMat = sparse(nRxns, 0);
            m.metFormulas = repmat({'C'}, nMets, 1);
        end

        function id = lastWarningId()
            [~, id] = lastwarn;
        end

        function m = nearlyFixedReactionModel()
            % The parallel model plus a reaction without metabolites.
            m = tSampling.parallelModel();
            m.rxns{end+1, 1} = 'NF';
            m.rxnNames = m.rxns;
            m.S(:, end+1) = 0;
            m.lb(end+1, 1) = 0;
            m.ub(end+1, 1) = 10;
            m.rev(end+1, 1) = 0;
            m.c(end+1, 1) = 0;
            m.grRules{end+1, 1} = '';
            m.rxnGeneMat(end+1, :) = 0;
        end

        function m = parallelModel()
            m = struct();
            m.id       = 'parallel';
            m.rxns     = {'EX_in';'R1';'R2';'EX_out'};
            m.rxnNames = m.rxns;
            m.mets     = {'a';'b'};
            m.metNames = {'a';'b'};
            m.metComps = [1;1];
            m.comps    = {'c'};
            m.compNames= {'cytosol'};
            %             EX_in  R1  R2  EX_out
            m.S = sparse([ 1    -1  -1   0;    % a
                           0     1   1  -1]);  % b
            m.lb  = [0;0;0;0];
            m.ub  = [10;6;6;1000];
            m.rev = [0;0;0;0];
            m.c   = [0;0;0;0];
            m.b   = zeros(2,1);
            m.genes = {}; m.grRules = {'';'';'';''}; m.rxnGeneMat = sparse(4,0);
            m.metFormulas = {'C';'C'};
        end

        function m = fixedReactionModel()
            % EX_in -> a -> R1 -> b -> EX_out, a -> R2 -> c -> FIX, with FIX
            % fixed at 700. R2 is forced to 700 through c, so R2 and FIX do
            % not move, and the null-space basis carries rounding noise on them.
            m = struct();
            m.id       = 'fixedReaction';
            m.rxns     = {'EX_in';'R1';'R2';'EX_out';'FIX'};
            m.rxnNames = m.rxns;
            m.mets     = {'a';'b';'c'};
            m.metNames = m.mets;
            m.metComps = [1;1;1];
            m.comps    = {'c'};
            m.compNames= {'cytosol'};
            %             EX_in  R1  R2  EX_out  FIX
            m.S = sparse([ 1    -1  -1   0       0;    % a
                           0     1   0  -1       0;    % b
                           0     0   1   0      -1]);  % c
            m.lb  = [0;0;0;0;700];
            m.ub  = [1000;1000;1000;1000;700];
            m.rev = zeros(5, 1);
            m.c   = zeros(5, 1);
            m.b   = zeros(3, 1);
            m.genes = {}; m.grRules = repmat({''}, 5, 1); m.rxnGeneMat = sparse(5, 0);
            m.metFormulas = repmat({'C'}, 3, 1);
        end

        function m = loopModel()
            m = struct();
            m.id       = 'loop';
            m.rxns     = {'R1';'R2';'R3';'R4'};
            m.rxnNames = m.rxns;
            m.mets     = {'a';'b'};
            m.metNames = {'a';'b'};
            m.metComps = [1;1];
            m.comps    = {'c'};
            m.compNames= {'cytosol'};
            %             R1  R2  R3  R4
            m.S = sparse([ 1  -1   1   0;    % a
                           0   1  -1  -1]);  % b
            m.lb  = [0;0;0;0];
            m.ub  = [10;100;100;100];
            m.rev = [0;0;0;0];
            m.c   = [0;0;0;1];
            m.b   = zeros(2,1);
            m.genes = {}; m.grRules = {'';'';'';''}; m.rxnGeneMat = sparse(4,0);
            m.metFormulas = {'C';'C'};
        end

    end
end
