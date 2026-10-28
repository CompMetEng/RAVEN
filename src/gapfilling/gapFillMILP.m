function [addedRxns, newModel, exitFlag] = gapFillMILP(model, universalModel, varargin)
% gapFillMILP  Objective-based gap-filling via a global growth-floor MILP.
%
% Finds the smallest set of reactions from a universal database that makes
% the model's objective (biomass) >= minGrowth, following the database
% repair of Kumar et al. 2007 (BMC Bioinformatics 8:212,
% doi:10.1186/1471-2105-8-212). Reactions already in the draft model keep
% their bounds: a reaction with the wrong direction is a curation issue and
% is not changed here.
%
% The universal model should contain all candidate reactions, including
% exchange and transport reactions if those repair classes are desired.
%
% This is a single-level MILP. It is computationally more expensive than
% gapFillFastLP but directly optimises for objective feasibility.
%
% Parameters
% ----------
% model : struct
%     draft RAVEN model to be gap-filled.
% universalModel : struct
%     universal reaction database (RAVEN model).
%
% Optional name-value parameters
% --------------------------------
% 'minGrowth' (default [])
%     minimum objective value to achieve. If empty, runs FBA on the merged
%     model and uses 10% of the maximum as the floor.
% 'bigM' (default 1000)
%     Big-M constant for the coupling constraints of universal reactions
%     with an infinite bound. A finite bound is used as it is.
% 'params' (default [])
%     solver parameter struct passed to optimizeProb. Fields such as
%     params.TimeLimit and params.intTol can be used to tune the MILP.
% 'verbose' (default true)
%     print progress messages.
%
% Returns
% -------
% addedRxns : cell
%     reaction IDs from universalModel that are selected by the MILP.
% newModel : struct
%     model with addedRxns incorporated.
% exitFlag : double
%     1 = optimal, -1 = infeasible, -2 = timeout or other failure.
%
% Notes
% -----
% The MILP is solved with the solver set by setRavenSolver. RAVEN's bundled GLPK solves it for small
% models (with a RAVEN:glpkMILP warning); for larger models Gurobi is much
% faster and more reliable (setRavenSolver('gurobi')).
%
% See also
% --------
% gapFillFastLP, gapFillTopological, fillGaps

% ---- Parse optional arguments ----
p = inputParser();
addParameter(p, 'minGrowth', [],  @(x) isempty(x) || isnumeric(x));
addParameter(p, 'bigM',      1000,      @isnumeric);
addParameter(p, 'params',    [],        @(x) isempty(x) || isstruct(x));
addParameter(p, 'verbose',   true,      @islogical);
parse(p, varargin{:});
minGrowth = p.Results.minGrowth;
M         = p.Results.bigM;
params    = p.Results.params;
verbose   = p.Results.verbose;

% ---- Initialise outputs ----
addedRxns = {};
newModel  = model;
exitFlag  = -1;

% Check that the model has an objective
if all(model.c == 0)
    error('gapFillMILP: model has no objective function. Set model.c.');
end

% ---- Merge draft + universal database ----
mergedModel = mergeModels({model; universalModel});
nMergedMets = numel(mergedModel.mets);
nMergedRxns = numel(mergedModel.rxns);

% Identify draft reactions in the merged model
isDraft = ismember(mergedModel.rxns, model.rxns);
isUniv  = ~isDraft;

% Universal reactions that can carry flux (non-trivially blocked)
canCarry = (mergedModel.ub > 0 | mergedModel.lb < 0);
isUniv   = isUniv & canCarry;
univIdx  = find(isUniv);     % positions in merged model
nUniv    = numel(univIdx);

if verbose
    fprintf('gapFillMILP: merged model has %d mets, %d rxns (%d draft, %d universal).\n', ...
        nMergedMets, nMergedRxns, sum(isDraft), nUniv);
    fprintf('gapFillMILP: %d database candidates.\n', nUniv);
end

% ---- Determine minGrowth ----
if isempty(minGrowth)
    % Run FBA on merged model (unconstrained) to get max growth.
    probFBA.A      = mergedModel.S;
    probFBA.a      = probFBA.A;
    probFBA.b      = zeros(nMergedMets, 1);
    probFBA.csense = repmat('E', 1, nMergedMets);
    probFBA.c      = -mergedModel.c;   % negate: osense=1 minimises
    probFBA.osense = 1;
    probFBA.lb     = mergedModel.lb;
    probFBA.ub     = mergedModel.ub;
    probFBA.vartype = repmat('C', 1, nMergedRxns);
    solFBA = optimizeProb(probFBA, params, false);
    if ~checkSolution(solFBA)
        if verbose
            fprintf('gapFillMILP: FBA on merged model is infeasible — cannot gap-fill.\n');
        end
        return;
    end
    maxGrowth = mergedModel.c' * solFBA.full(1:nMergedRxns);
    if maxGrowth <= 0
        if verbose
            fprintf('gapFillMILP: merged model objective is non-positive (%.4g). Cannot gap-fill.\n', ...
                maxGrowth);
        end
        return;
    end
    minGrowth = 0.1 * maxGrowth;
    if verbose
        fprintf('gapFillMILP: setting minGrowth = %.4g (10%% of max %.4g).\n', minGrowth, maxGrowth);
    end
end

% ---- Build MILP ----
% Variable ordering:
%   cols 1 : nMergedRxns          → v (all fluxes, continuous)
%   cols nMergedRxns+1 : end      → y_db (binary, universal reactions)
nVar = nMergedRxns + nUniv;

% --- Variable bounds ---
lb = mergedModel.lb;
ub = mergedModel.ub;

% Universal reactions: allow them to be forced to 0 (extend lb toward 0)
lb(univIdx) = min(mergedModel.lb(univIdx), 0);

% Binary variable bounds
lb = [lb; zeros(nUniv, 1)];
ub = [ub; ones(nUniv, 1)];

% --- Objective: number of added universal reactions ---
c = zeros(nVar, 1);
c(nMergedRxns + 1 : end) = 1;

% --- Constraint blocks ---
% (1) Steady state: S * v = 0   [nMergedMets rows, 'E']
A_SS = [mergedModel.S, sparse(nMergedMets, nUniv)];

% (2) Growth floor: c' * v >= minGrowth   [1 row, 'G']
A_growth = [sparse(mergedModel.c'), sparse(1, nUniv)];

% (3) Universal upper coupling: v_k - ub_k * y_db_k <= 0   [nUniv rows, 'L']
%     (forces v_k <= 0 when y_db_k = 0)
% Coupling coefficients are the reaction bounds. An infinite bound is replaced
% by the Big-M constant: an unbounded universal reaction (ub = Inf, RAVEN's own
% default when no upper-bound annotation is set) would otherwise put an Inf
% coefficient into the constraint matrix, and Inf*0 (y_db_k = 0) is NaN, not 0.
yDbCols = nMergedRxns + (1:nUniv);
uniUb = mergedModel.ub(univIdx);
uniUb(isinf(uniUb)) = M;
uniLb = mergedModel.lb(univIdx);
uniLb(isinf(uniLb)) = -M;
A_uniUp = sparse(...
    [1:nUniv, 1:nUniv], ...
    [univIdx(:)', yDbCols], ...
    [ones(1, nUniv), -uniUb'], ...
    nUniv, nVar);

% (4) Universal lower coupling: v_k - lb_k * y_db_k >= 0   [nUniv rows, 'G']
%     (forces v_k >= 0 when y_db_k = 0; together with (3): v_k = 0 when y_db_k = 0)
A_uniLo = sparse(...
    [1:nUniv, 1:nUniv], ...
    [univIdx(:)', yDbCols], ...
    [ones(1, nUniv), -uniLb'], ...
    nUniv, nVar);

% --- Assemble ---
prob.A = [A_SS; A_growth; A_uniUp; A_uniLo];
prob.a = prob.A;
prob.b = [zeros(nMergedMets, 1); minGrowth; zeros(2 * nUniv, 1)];
prob.csense = [repmat('E', 1, nMergedMets), 'G', ...
               repmat('L', 1, nUniv), ...
               repmat('G', 1, nUniv)];
prob.c      = c;
prob.osense = 1;   % minimise the number of added reactions
prob.lb     = lb;
prob.ub     = ub;
prob.vartype = repmat('C', 1, nVar);
prob.vartype(nMergedRxns + 1 : end) = 'B';   % y_db is binary

if isempty(params)
    params.intTol   = 1e-9;
    params.TimeLimit = 600;
end

% ---- Solve ----
if verbose
    fprintf('gapFillMILP: solving MILP (%d vars, %d binary, %d constraints)...\n', ...
        nVar, nUniv, size(prob.A, 1));
end
sol = optimizeProb(prob, params, verbose);

if ~checkSolution(sol)
    if verbose
        fprintf('gapFillMILP: MILP infeasible or timed out.\n');
    end
    exitFlag = -1;
    return;
end

% ---- Extract solution ----
binThreshold = 1e-3;
yDb = sol.full(nMergedRxns + 1 : end);
dbSelected = univIdx(yDb > binThreshold);

% Map back to original IDs. dbSelected is already restricted to universal
% reactions via univIdx/isDraft above (an id-membership test against
% model.rxns, not a name it could have been renamed away from), so no
% further filtering against universalModel.rxns is needed here -- doing so
% would silently drop a universal reaction that mergeModels renamed
% because its id collided with one already in the draft model.
addedRxns = mergedModel.rxns(dbSelected);

if verbose
    fprintf('gapFillMILP: added %d universal reaction(s).\n', numel(addedRxns));
end

% ---- Build new model ----
newModel = model;

% Add universal reactions. Extracted from mergedModel by the dbSelected
% indices directly, not by looking addedRxns' names back up in
% universalModel: that lookup fails whenever mergeModels renamed one of
% them because its id collided with one already in the draft model, in
% which case setdiff would keep the (unrenamed) original in universalModel
% marked for removal, silently dropping the very reaction just selected.
if ~isempty(dbSelected)
    toRemove = true(nMergedRxns, 1);
    toRemove(dbSelected) = false;
    trimmedUniv = removeReactions(mergedModel, toRemove);
    newModel = mergeModels({newModel; trimmedUniv});
end

exitFlag = 1;
end
