function map = defaultCompartmentMap()
% defaultCompartmentMap  Default map from compartment labels to model compartment ids.
%
% Returns a containers.Map from lower-case compartment labels of predictors and databases (as
% used by DeepLoc, MULocDeep, COMPARTMENTS and UniProt) to compartment ids, for yeast and fungal
% models (for example yeast-GEM codes). The map is passed to parseScores or getUniProtScores as
% 'compartmentMap'.
%
% The map holds three decisions:
%
% - Compartments present. A label absent from the map is dropped along with its score column, so
%   the map lists only the compartment ids of the model. 'plastid' is absent because fungi have
%   no plastid; a plant model adds it.
% - Labels that name the same compartment. Predictors differ in wording, so 'cytoplasm' and
%   'cytosol' both map to 'c', and 'mitochondrion', 'mitochondria' and 'mitochondrial' all map
%   to 'm'. Merged columns are combined by their maximum.
% - Distinct compartments that are collapsed. 'lysosome' maps to 'v' because a yeast model has
%   no lysosome and the vacuole is its functional equivalent. A model that separates the two
%   needs a different mapping.
%
% Changes to the mapping go into a copy of the map, which leaves this function unedited, for
% example:
%
%     map = defaultCompartmentMap;
%     map('cell wall') = 'ce';
%     GSS = parseScores(file, 'compartmentMap', map);
%
% Returns
% -------
% map : containers.Map
%     label (char, lower case) -> compartment id (char).
%
% See also
% --------
% parseScores, getUniProtScores, predictLocalization

labels = {'cytoplasm','cytosol','nucleus','nucleoplasm','mitochondrion','mitochondria', ...
          'mitochondrial','peroxisome','endoplasmic reticulum','golgi apparatus','golgi', ...
          'vacuole','lysosome/vacuole','lysosome','extracellular','extracellular space', ...
          'extracellular region','secreted','cell membrane','plasma membrane','cell envelope', ...
          'lipid particle','lipid droplet'};
ids    = {'c','c','n','n','m','m', ...
          'm','p','er','g','g', ...
          'v','v','v','e','e', ...
          'e','e','ce','ce','ce', ...
          'lp','lp'};
map = containers.Map(labels, ids);
end
