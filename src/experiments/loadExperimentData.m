function experimentData = loadExperimentData(dataDir)
% loadExperimentData  Read the measured conditions of a model from a folder.
%
% The folder holds up to three tab-delimited files, each optional:
%
% - fluxData.tsv : flux conditions, in which the carbon-source uptake is limited to the measured
%   rate and the model is compared with the measured growth rate and exchange fluxes.
% - maxGrowth.tsv : maximum-growth conditions, in which the uptake is unlimited and only the
%   maximal growth rate is compared. Same layout as fluxData.tsv (see loadFluxData).
% - zeroFlux.tsv : the header Rxns and one reaction id per line, for reactions that carry no flux
%   in any condition. Flux through them is penalised in the flux conditions.
%
% A missing file gives [] ({} for zeroFlux). The files feed simulateConditions and
% experimentDistance, which condense the mismatch between model and measurements into one number:
% computeCarbonWeights, simulateConditions (once per file), experimentDistance.
%
% The model objective must be the maximisation of growth: the functions do not set it, and in a
% maximum-growth condition the compared value is the maximum of the model objective.
% simulateConditions warns if the objective is not the growth reaction.
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
%     - zeroFlux : cell column with the reaction ids of zeroFlux.tsv, or {}.
%
% Examples
% --------
%     experimentData = loadExperimentData('data');
%
% See also
% --------
% loadFluxData, computeCarbonWeights, simulateConditions, experimentDistance, applyFluxDataConstraints

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
