function experimentData = loadExperimentData(dataDir)
% loadExperimentData  Load the measured conditions of a model from a folder.
%
% Reads fluxData.tsv, maxGrowth.tsv and zeroFlux.tsv from a folder. The
% first two have the layout described in loadFluxData; fluxData.tsv holds
% conditions with measured uptake and excretion fluxes, and maxGrowth.tsv
% holds one carbon source per condition, which is opened fully when the
% maximum growth rate is simulated. zeroFlux.tsv has the header Rxns and one
% reaction id per line. Each file is optional: a missing file gives [] (or an
% empty cell array for zeroFlux), not an error.
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
%     - zeroFlux : cell column with the reaction ids that carry no flux in
%       any condition, or {}.
%
% Examples
% --------
%     experimentData = loadExperimentData('data');
%
% See also
% --------
% loadFluxData, simulateConditions, experimentDistance

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
