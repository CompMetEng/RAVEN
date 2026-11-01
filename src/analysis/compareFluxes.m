function result = compareFluxes(model, fluxes1, fluxes2, varargin)
% compareFluxes  Find the reactions whose flux differs between two flux vectors.
%
% Compares two flux vectors from the same model and returns all reactions
% that changed, sorted by absolute flux difference. Classifies reactions as
% turned on, turned off, reversed direction, or changed magnitude.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% fluxes1 : double
%     reference flux vector (condition A, such as wild-type).
% fluxes2 : double
%     comparison flux vector (condition B, such as a mutant or knock-in).
%
% Name-Value Arguments
% --------------------
% cutoff : double
%     minimum |flux| for a reaction to be active, and minimum |flux2 - flux1|
%     for a reaction to be changed (default 1e-8).
% metaboliteList : cell
%     cell array of metabolite names. Only reactions involving at least one
%     of these metabolites are compared; all other reactions are excluded
%     from every field of the result. Matching is against model.metNames and
%     is case insensitive. A name that matches no metabolite is reported as
%     a warning (default all reactions).
% nMax : double
%     maximum rows printed in the summary table (default 20).
% verbose : logical
%     print the summary table (default true).
%
% Returns
% -------
% result : struct
%     structure with fields:
%
%     - turnedOn : cell array with the IDs of reactions that are active in
%       fluxes2 and not in fluxes1.
%     - turnedOff : cell array with the IDs of reactions that are active in
%       fluxes1 and not in fluxes2.
%     - flipped : cell array with the IDs of reactions that are active in
%       both conditions and change sign.
%     - changed : struct with fields for all reactions with |flux2 - flux1|
%       above cutoff, sorted by |flux2 - flux1| descending:
%
%       - rxn : cell array with the reaction IDs.
%       - flux1 : double array with the flux in condition A.
%       - flux2 : double array with the flux in condition B.
%       - absDelta : double array with |flux2 - flux1|.
%       - relChange : double array with (flux2 - flux1)/|flux1|. NaN when
%         |flux1| is below cutoff.
%       - type : cell array with 'on', 'off', 'flip' or '' (changed
%         magnitude) for each reaction.
%
% Examples
% --------
%     solA = solveLP(model);
%     modelAna = setParam(model, 'ub', 'O2exchange', 0);
%     solB = solveLP(modelAna);
%     result = compareFluxes(model, solA.x, solB.x)
%     result = compareFluxes(model, solA.x, solB.x, 'nMax', 50)
%     result = compareFluxes(model, solA.x, solB.x, 'metaboliteList', {'ATP','NADH'})
%
% See also
% --------
% walkFluxes, traceFluxPath, modelSummary

p = parseRAVENargs(varargin, {'cutoff',1e-8; 'metaboliteList',[]; 'nMax',20; 'verbose',true});
cutoff         = p.cutoff;
metaboliteList = p.metaboliteList;
nMax           = p.nMax;
verbose        = p.verbose;

f1    = fluxes1(:);
f2    = fluxes2(:);
nRxns = numel(model.rxns);

if numel(f1) ~= nRxns || numel(f2) ~= nRxns
    error('RAVEN:badInput', ...
        'Flux vectors must have %d elements (numel(model.rxns)).', nRxns);
end

active1 = abs(f1) > cutoff;
active2 = abs(f2) > cutoff;
delta   = f2 - f1;

% Restrict to reactions involving one of metaboliteList, if given. Reactions
% left out are treated as unchanged, so they appear in no output field.
inList = true(nRxns, 1);
if ~isempty(metaboliteList)
    metaboliteList = convertCharArray(metaboliteList);
    inList = false(nRxns, 1);
    notFound = {};
    for i = 1:numel(metaboliteList)
        metIdx = find(strcmpi(metaboliteList{i}, model.metNames));
        if isempty(metIdx)
            notFound{end+1} = metaboliteList{i}; %#ok<AGROW>
        else
            inList(any(model.S(metIdx,:) ~= 0, 1)) = true;
        end
    end
    if ~isempty(notFound)
        warning('RAVEN:warning', '%s', ...
            ravenList('No metabolite in the model is named:', notFound));
    end
    active1(~inList) = false;
    active2(~inList) = false;
    delta(~inList)   = 0;
end

onMask   = ~active1 &  active2;
offMask  =  active1 & ~active2;
flipMask =  active1 &  active2 & sign(f1) ~= sign(f2);

result.turnedOn  = model.rxns(onMask);
result.turnedOff = model.rxns(offMask);
result.flipped   = model.rxns(flipMask);

% All changed reactions sorted by |delta| descending
changedMask = abs(delta) > cutoff;
[~, si]     = sort(abs(delta(changedMask)), 'descend');
idx = find(changedMask);
idx = idx(si);

% Relative change: (f2-f1)/|f1|; NaN when reference flux is near-zero
denom = abs(f1(idx));
denom(denom < cutoff) = NaN;
relChange = delta(idx) ./ denom;

% Type label per changed reaction
changedTypes = repmat({''}, numel(idx), 1);
for k = 1:numel(idx)
    ri = idx(k);
    if onMask(ri)
        changedTypes{k} = 'on';
    elseif offMask(ri)
        changedTypes{k} = 'off';
    elseif flipMask(ri)
        changedTypes{k} = 'flip';
    end
end

result.changed.rxn       = model.rxns(idx);
result.changed.flux1     = f1(idx);
result.changed.flux2     = f2(idx);
result.changed.absDelta  = abs(delta(idx));
result.changed.relChange = relChange;
result.changed.type      = changedTypes;

if verbose
    printTable_(model, result, nMax, nnz(inList));
end
end

%--------------------------------------------------------------------------
function printTable_(model, r, nMax, nConsidered)
    W   = 72;
    bar = repmat('=', 1, W);
    sep = repmat('-', 1, W);
    nC  = numel(r.changed.rxn);

    fprintf('\n%s\n', bar);
    fprintf('  Flux comparison: %d of %d reactions changed\n', nC, nConsidered);
    nOn  = numel(r.turnedOn);
    nOff = numel(r.turnedOff);
    nFlp = numel(r.flipped);
    if nOn + nOff + nFlp > 0
        fprintf('  (on: %d  |  off: %d  |  reversed: %d)\n', nOn, nOff, nFlp);
    end

    if nC == 0
        fprintf('%s\n\n', bar);
        return
    end

    nShow    = min(nC, nMax);
    showRxns = r.changed.rxn(1:nShow);
    eqns     = constructEquations(model, showRxns);

    fprintf('%s\n', sep);
    fprintf('  %-20s  %10s  %10s  %10s  %7s  %-5s  %-25s  %s\n', ...
        'Reaction', 'Cond.1', 'Cond.2', '|dflux|', 'drel%', 'Type', 'Name', 'Equation');
    fprintf('  %s\n', sep);

    for i = 1:nShow
        rel = r.changed.relChange(i);
        if isnan(rel)
            relStr = '    new';
        else
            relStr = sprintf('%+6.0f%%', rel * 100);
        end
        ri = find(strcmp(model.rxns, showRxns{i}), 1);
        nm = '';
        if isfield(model,'rxnNames') && ri <= numel(model.rxnNames)
            nm = char(model.rxnNames{ri});
        end
        if numel(nm) > 25, nm = [nm(1:22) '...']; end
        fprintf('  %-20s  %+10.4g  %+10.4g  %10.4g  %7s  %-5s  %-25s  %s\n', ...
            showRxns{i}, r.changed.flux1(i), r.changed.flux2(i), ...
            r.changed.absDelta(i), relStr, r.changed.type{i}, nm, eqns{i});
    end
    if nC > nMax
        fprintf('  ... (%d more, use nMax to show all)\n', nC - nMax);
    end
    fprintf('%s\n\n', bar);
end
