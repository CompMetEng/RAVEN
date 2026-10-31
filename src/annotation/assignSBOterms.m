function model = assignSBOterms(model, varargin)
% assignSBOterms  Assign SBO terms to metabolites and reactions.
%
% Assigns SBO terms to metabolites and reactions following a generic rule
% set that opts parameterises. The function is the MATLAB counterpart of
% yeast-GEM's addSBOterms.
%
% SBO is written via editMiriam(..., 'fill'), which preserves pre-existing
% SBO annotations. The terms are stored in model.metMiriams and
% model.rxnMiriams under the namespace 'sbo' (there is no separate SBO
% field); they can be read with extractMiriam(model.rxnMiriams, 'sbo').
% exportModel writes them as the SBML sboTerm attribute.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% opts : struct
%     structure with any of the following fields. Missing or empty fields
%     take the defaults shown:
%
%     - biomassMetNames : {'biomass','DNA','RNA','protein','carbohydrate',
%       'lipid','cofactor','ion'}
%     - biomassMetSuffixes : {' backbone',' chain'}
%     - biomassRxnName : 'biomass pseudoreaction'
%     - ngamRxnName : 'non-growth associated maintenance reaction'
%     - pseudoreactionSubstrings : {'pseudoreaction','SLIME rxn'}
%     - onlyLastReactionForPseudo : false. If true, the biomass, NGAM and
%       pseudoreaction terms are applied only to the last reaction of the
%       model. This reproduces yeast-GEM's addSBOterms.m, which applies the
%       pseudoreaction rule to the last reaction only, because of a typo in
%       the loop header (for i=numel(...)).
%     - transportRxnIdxs : indexes of the transport reactions (default:
%       the reactions identified by getTransportRxns).
%
% Returns
% -------
% model : struct
%     model structure with the SBO terms added to model.metMiriams and
%     model.rxnMiriams.
%
% Examples
% --------
%     model = assignSBOterms(model);
%     model = assignSBOterms(model, struct('onlyLastReactionForPseudo', true));
%
% Notes
% -----
% Metabolites:
%
%     SBO:0000649 (Biomass) when the metabolite name is in
%     opts.biomassMetNames, or ends with any of opts.biomassMetSuffixes.
%     Otherwise SBO:0000247 (Simple chemical).
%
% Reactions (default → override → pseudoreaction override):
%
%     SBO:0000176 (Metabolic reaction) default.
%     Single-reactant reactions become:
%         SBO:0000627 (exchange) if the lone metabolite is extracellular
%             (compartment 'e' or compartment name 'extracellular'),
%         SBO:0000632 (sink) if coef < 0,
%         SBO:0000628 (demand) otherwise.
%     Transport reactions (opts.transportRxnIdxs or, by default, the
%         reactions identified by getTransportRxns: the same metabolite
%         name more than once in a single reaction) → SBO:0000655.
%     Reactions whose name matches opts.biomassRxnName → SBO:0000629.
%     Reactions whose name matches opts.ngamRxnName → SBO:0000630.
%     Reactions whose name contains any of opts.pseudoreactionSubstrings
%         → SBO:0000395.

p=parseRAVENargs(varargin, {'opts',[]});
opts=p.opts;
if isempty(opts)
    opts = struct();
end
opts = applyDefaults(opts);

% Metabolite SBO ------------------------------------------------------
metsSBO = cell(size(model.mets));
for i = 1:length(model.mets)
    metName = model.metNames{i};
    if any(strcmp(opts.biomassMetNames, metName)) || endsWithAny(metName, opts.biomassMetSuffixes)
        metsSBO{i} = 'SBO:0000649';
    else
        metsSBO{i} = 'SBO:0000247';
    end
end

% Reaction SBO --------------------------------------------------------
rxnSBO = cell(size(model.rxns));
rxnSBO(:) = {'SBO:0000176'};

% Single-reactant reactions
reactantNumber = sum(model.S ~= 0, 1);
singleRxns = find(reactantNumber == 1);
for k = 1:numel(singleRxns)
    idx = singleRxns(k);
    metRow = find(model.S(:, idx));
    compName = model.compNames{model.metComps(metRow)};
    compShort = model.comps{model.metComps(metRow)};
    if strcmp(compShort, 'e') || strcmp(compName, 'extracellular')
        rxnSBO{idx} = 'SBO:0000627';
    elseif sum(model.S(:, idx)) < 0
        rxnSBO{idx} = 'SBO:0000632';
    else
        rxnSBO{idx} = 'SBO:0000628';
    end
end

% Transport reactions
if isfield(opts, 'transportRxnIdxs') && ~isempty(opts.transportRxnIdxs)
    transportIdxs = opts.transportRxnIdxs;
else
    transportIdxs = getTransportRxns(model);
end
rxnSBO(transportIdxs) = {'SBO:0000655'};

% Pseudoreaction overrides
if opts.onlyLastReactionForPseudo
    pseudoTargets = numel(model.rxns);
else
    pseudoTargets = 1:numel(model.rxns);
end
for ii = pseudoTargets
    name = model.rxnNames{ii};
    if strcmp(name, opts.biomassRxnName)
        rxnSBO{ii} = 'SBO:0000629';
    elseif strcmp(name, opts.ngamRxnName)
        rxnSBO{ii} = 'SBO:0000630';
    else
        for k = 1:numel(opts.pseudoreactionSubstrings)
            if contains(name, opts.pseudoreactionSubstrings{k})
                rxnSBO{ii} = 'SBO:0000395';
                break;
            end
        end
    end
end

model = editMiriam(model, 'met', 'all', 'sbo', metsSBO, 'fill');
model = editMiriam(model, 'rxn', 'all', 'sbo', rxnSBO, 'fill');
end


function opts = applyDefaults(opts)
defaults = struct( ...
    'biomassMetNames', {{'biomass','DNA','RNA','protein','carbohydrate','lipid','cofactor','ion'}}, ...
    'biomassMetSuffixes', {{' backbone',' chain'}}, ...
    'biomassRxnName', 'biomass pseudoreaction', ...
    'ngamRxnName', 'non-growth associated maintenance reaction', ...
    'pseudoreactionSubstrings', {{'pseudoreaction','SLIME rxn'}}, ...
    'onlyLastReactionForPseudo', false);
fields = fieldnames(defaults);
for k = 1:numel(fields)
    f = fields{k};
    if ~isfield(opts, f) || isempty(opts.(f))
        opts.(f) = defaults.(f);
    end
end
end


function tf = endsWithAny(s, suffixes)
tf = false;
for i = 1:numel(suffixes)
    if endsWith(s, suffixes{i})
        tf = true;
        return;
    end
end
end
