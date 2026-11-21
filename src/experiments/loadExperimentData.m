function experimentData = loadExperimentData(dataDir)
% loadExperimentData  Read the measured conditions of a model from a folder.
%
% Suppose you measured the growth rate and some exchange fluxes of an organism in one or more
% growth conditions, and want to know how well a model reproduces them: to compare model versions,
% or as the objective when parameters of a model are adjusted to the data. The functions of the
% experiments folder simulate each condition in the model and condense the mismatch into a single
% number, the distance.
%
% Two kinds of conditions are scored:
%
% - A flux condition limits the uptake of the carbon source to its measured rate. The growth rate
%   and the measured exchange fluxes of the model are compared with the measurements. The
%   conditions are stored in fluxData.tsv.
% - A maximum-growth condition leaves the uptake of the carbon source unlimited. Only the maximal
%   growth rate of the model is compared with the measured growth rate. The conditions are stored
%   in maxGrowth.tsv.
%
% Both files have the layout described in loadFluxData. A third file, zeroFlux.tsv, lists reactions
% that carry no flux in any condition: the header Rxns and one reaction id per line. A model that
% carries flux in them is penalised in the flux conditions. Each file is optional; a missing file
% gives [] (or an empty cell array for zeroFlux), not an error.
%
% IMPORTANT: the model objective must be the maximisation of the growth. The functions do not set
% the objective: simulateConditions solves each condition with the objective that the model has
% (model.c, which is maximised), and the growth rate is the flux of the biomass reaction in that
% solution. In a maximum-growth condition the number that is compared with the measurement is
% therefore the maximum of the model objective, which is the maximal growth rate only if the
% objective is the biomass reaction. Set the objective with setParam(model, 'obj', biomassRxn, 1)
% first; simulateConditions warns once per call if the objective is not the biomass reaction.
%
% Workflow:
%
% 1. loadExperimentData reads the three files of a folder.
% 2. computeCarbonWeights gives each reaction the weight that puts all fluxes on the same carbon
%    scale.
% 3. simulateConditions simulates the conditions of one file in the model.
% 4. experimentDistance compares the simulations of both kinds of conditions with the
%    measurements and returns the distance. datasetRmse does the same for one file.
%
% applyFluxDataConstraints is independent of the distance: it constrains a model to one measured
% condition, for example before a flux variability analysis.
%
% Fluxes are in mmol/gDCW/h and growth rates in 1/h. A negative exchange flux is an uptake and a
% positive one an excretion.
%
% Parameters
% ----------
% dataDir : char
%     folder that holds the files.
%
% Returns
% -------
% experimentData : struct
%     structure with the fields:
%
%     - fluxData : the structure of loadFluxData for fluxData.tsv, or [].
%     - maxGrowth : the structure of loadFluxData for maxGrowth.tsv, or [].
%     - zeroFlux : cell column with the reaction ids that carry no flux in any condition, or {}.
%
% Examples
% --------
%     experimentData = loadExperimentData('data');
%
% See also
% --------
% loadFluxData, computeCarbonWeights, simulateConditions, experimentDistance

experimentData.fluxData   = loadOptionalFluxData(fullfile(dataDir, 'fluxData.tsv'));
experimentData.maxGrowth  = loadOptionalFluxData(fullfile(dataDir, 'maxGrowth.tsv'));
experimentData.zeroFlux   = loadZeroFlux(fullfile(dataDir, 'zeroFlux.tsv'));
end

function fluxData = loadOptionalFluxData(file)
fluxData = [];
if isfile(file)
    fluxData = loadFluxData(file);
end
end

function zeroFlux = loadZeroFlux(file)
zeroFlux = {};
if ~isfile(file)
    return
end
lines = splitlines(string(fileread(file)));
lines = strtrim(erase(lines, char(13)));
lines = lines(2:end);
lines = lines(strlength(lines) > 0);
zeroFlux = cellstr(lines);
end
