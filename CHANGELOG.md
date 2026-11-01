# Changelog

## 3.0.0b2 — 2026-11-01

Second beta toward 3.0. Changes below are relative to 3.0.0b1.

### Repository layout and downloads

* The MATLAB code is in `src/` (was: folders in the repository root), the tests in `tests/` (was: `testing/`), the generated help pages in `docs/html/` (was: `doc/`), and the tutorials in `tutorials/` with data under `tutorials/data/{models,omics,sequences}` (was: `tutorial/`). The external software is in `src/external` (was: `src/software`). `checkRaven`, `downloadRavenBinaries`, `getBlast`, `getDiamond`, `getKEGGModelForOrganism`, `optimizeProb`, `setRavenSolver` and `updateDocumentation` use the new paths. This matters for code that adds RAVEN folders to the path or addresses binaries by path.
* `parseHPA`, `parseHPArna` and `scoreModel` are in `INIT/` (was: `omics/`, which no longer exists). `tOmics.m` is merged into `tINIT.m`.
* `fetchToolboxAsset` replaces `fetchRavenDataAsset`, and assets are downloaded from `CompMetEng/toolbox-assets` (was: raven-data). `downloadRavenBinaries`, the KEGG functions, `getPhylDist` and `checkRaven` use it.
* `fetchToolboxAsset` downloads to a `.part` file and checks the SHA-256 sum listed in the `manifest-v3` of toolbox-assets. A mismatch deletes the file and raises an error; an asset that the manifest does not list, or a manifest that cannot be read, gives a warning.
* The repository, documentation and version-check links point to `CompMetEng`.

### I/O

* `writeYAMLmodel` writes an ec reaction without enzymes (a row of zeros in `model.ec.rxnEnzMat`) with no `enzymes` entry (was: `- enzymes: !!omap` with nothing under it, which is not valid YAML). `readYAMLmodel` reads a missing `enzymes` entry as a reaction without enzymes, and still reads `- enzymes: !!omap` with no entries and `- enzymes: !!omap []` from files written earlier.

* `readYAMLmodel` reads notes maps (as written by cobrapy and raven-toolbox) into `metNotes` and `rxnNotes` (was: the error "Unknown entry"), reads folded scalars whose continuation starts with `-`, and keeps one `model.subSystems` entry per reaction (was: empty entries were dropped, so later subsystems sat on the wrong reaction).
* `writeYAMLmodel` writes the enzymes of each `ec-rxns` entry indented under `- enzymes` (was: at the indentation of the key, which ruamel.yaml could not read), quotes the YAML 1.1 booleans (`y`, `n`, `yes`, `no`, `on`, `off` in any case), and writes non-integers with the shortest digits that read back exactly. `printOrange` prints `\n` as a line break.
* `exportForGit`: the txt table is headed `Rxn id`, writes bounds and objective with the digits needed to read them back and pads no field. The output folder is created with its parents also when `subDirs` is false, and the logical flags are validated.

### Manipulation

* `mergeCompartments`: a merged metabolite is the copy already in compartment `s`, if any,
  otherwise the copy whose id sorts first, with every field of that copy (was: the first copy
  in model order), so `ACCOA_c` and `ACCOA_m` become `ACCOA_c` whatever their order.
* `simplifyModel` with `constrainReversible`: on a model without a reversible reaction this step
  changes nothing (was: an error "logical indices ... outside of the array bounds", as
  `getAllowedBounds` took the empty reaction list for all reactions). `prepINITModel` therefore
  runs on a model whose reactions are all irreversible after the first simplification, with the
  result of raven-toolbox's `prep_init_model`.
* `removeGenes`: the help states that, when at least one of the requested genes is in the model,
  genes no remaining reaction uses are removed as well; when none is, no gene, reaction or bound
  changes.

### Gap-filling

* `gapFillTopological`: without `targets` and without an objective reaction it raises
  `RAVEN:badInput` (was: an error without identifier). Help states the default seeds and
  targets.
* `checkProduction` with `checkNeededForProduction`: the returned `model` holds excretion
  reactions only (was: the excretion reaction of the last metabolite in `notProduced` was
  returned as an uptake reaction). This matters for code that uses the fifth output.

### Curation

* `renameModelGenes` with a file path as `geneTable`: the table is read and the genes are
  renamed (was: an error `MATLAB:invalidConversion` when the message that names the file
  was printed).

### Queries

* `getIndexes` with type `'metnames'` and `returnLogical` true: raises `RAVEN:badInput` (was:
  the error "Unable to use a value of type cell as an index"). A name can match several
  metabolites, so no single logical vector is defined.

### Analysis

* `findGeneDeletions` with `'dgd'`: `grRatioMuts` has one element per gene pair, as the help
  states (was: twice as many, the second half zero). This matters for code that indexes
  `grRatioMuts` by the number of its elements.

### Context-specific models (ftINIT / tINIT)

* `ftINIT` with `resolveTies`: each main-extraction step takes its fluxes from the
  least-total-flux solution for the tie-broken on/off state, so the direction and magnitude
  at which the next step forces a reaction no longer depend on the solver's vertex.
  That phase stops at a relative or absolute gap of 1e-6.
* `ftINIT` with `resolveTies`: among the least-total-flux solutions, each main-extraction
  step takes its fluxes from the one with the least id-rank-weighted backward flux, so the
  direction of a flux that carries the same total either way no longer depends on the
  solver. The id-rank phases rank a merged reaction by the smallest original reaction id in
  its group, so they no longer depend on which member names it.
* `ftINITInternalAlg` with `resolveTies`: when the primary solve stops at the time limit, the
  tie-break phases are skipped and the primary solution is kept, as in raven-toolbox (was: the
  phases ran on the unproven incumbent). The count and id-rank phases get half of
  `params.TimeLimit` each, rounded down and at least 1 s (was: the full limit).
* `ftINIT` with `resolveTies`: the last step skips the minimal-flux and direction phases, as
  no step reads its fluxes (new `ftINITInternalAlg` argument `pinFluxes`). The reaction set
  does not change; the fluxes and `metProduction` of the last step are those of the id-rank
  phase.
* `ftINITInternalAlg`: new optional positional arguments `rankKeys` (14th, the key each
  reaction is sorted by in the id-rank phases, default `model.rxns`) and `pinFluxes` (15th,
  default `true`). The "tie resolution did not converge" warning of `resolveTies` also names
  `phase2c-minflux` and `phase2d-direction` when the minimal-flux or direction phase does not
  prove optimal.
* `mergeLinear`: groups are numbered 1, 2, ... in the order of the first reaction of each
  group in `origRxnIds` (was: the order in which the reactions were merged).
* `mergeLinear`: a reaction flipped after the first merge pass is the one marked in
  `reversedRxns` (was: the reaction at the same index of the original model).

### Utilities

* `grRuleToDNF`: help states the grammar it reads, that of `parseGrRule` (`and`/`or` in any
  case, `&`, `&&`, `|`, `||`, operators as whole words, so `G1&G2` is one gene id), and that a
  malformed rule raises `RAVEN:badGrRule`.

* `checkRaven`: when the binaries cannot be downloaded, the message names the
  `RAVEN-<tag>-full.zip` bundle attached to each release for an offline install (was:
  "the `*-binaries` RAVEN release").

### Homology and curation

* `getModelFromHomology` has the options `complexPolicy` (`flag`, `keep` or `drop`: how a reaction is treated when an `and`-linked subunit has no ortholog), `keepGeneFree` (reactions of the template without a grRule are transferred) and `preserveNotes` (`rxnNotes`, `rxnReferences` and `rxnConfidenceScores` are kept). A gene that several template genes reach appears once per grRule. The draft no longer carries the version, date and description of the templates. It no longer fails when a preferred order removes no gene.
* `getModelFromHomology` raises an error that names the gene ids of the new organism that split in a grRule (whitespace, parentheses, `|` or `&`, for example the UniProt id `sp|P33284|HXK_KLULA`) (was: the ids were accepted).
* New: `hitsFromSuppOrth` reads a `suppOrth.tsv` into a blast structure for `getModelFromHomology`, in both directions, with filters on support, predictors, identity and E-value. `runSuppOrth` runs the `suppOrth` command and reads its output.
* `curateModelFromTables` reads the tables by header name in any order, with optional columns, and writes every supplied field to existing entities (`metNotes` and `rxnReferences` were never written). The MIRIAM columns follow the core columns, and grRules go to the right reaction when existing and new reactions are mixed. The option `overwriteWithEmpty` is new, an empty cell counts as not given, and the input is validated: a file, column, gene-rule or repeated-id error raises `RAVEN:badInput`.
* `scaleBiomassFraction` raises an error when the required fraction is 0 or less, or when the balancing component has no current fraction (was: the fractions were silently inverted). The new argument `targetTotal` sets the balanced total of the biomass.
* `setGAM` and `scaleBiomassPseudoreaction` raise an error before the model changes when the NGAM reaction or the proton metabolite is missing. `setGAM` refuses `ngamValue` without `ngamRxn`, and an empty `ngamRxn` counts as not given.
* `sortModel` keeps the directions of the reactions unless `sortReversible` is requested. `generateNewIds` takes the numerically largest id and counts only finite integer suffixes.
* `replaceMets`: the final `contractModel` call is restricted to the surviving metabolite (was: unrelated duplicate reactions were merged as well).
* `addIdentifierPrefix` prefixes whole gene ids in grRules (was: an id with a dash, dot or colon, for example `5-HT1A`, got a prefix at every word start, so `exportModel` failed).
* `getIndexes` with an ambiguous `name[comp]` raises an error that names the matching metabolites.
* `gapFillMILP` keeps the lower bound of unreversed candidate reactions.
* `gapReport` always reports the metabolites that are produced or consumed with the boundary closed. It closes the reactions of `getExchangeRxns` when the model has no `unconstrained` field, and stops when the closed model is infeasible.

### Solvers

* `optimizeProb` returns the logical `res.hitTimeLimit` for every solver. `getMinNrFluxes`, `ftINITFillGapsMILP`, `runINIT` and `ftINITInternalAlg` report a time limit as such, with a warning (was: every solver except Gurobi reported "The problem is infeasible"). The exit code -1 means that no usable solution came back, and -2 means a solution that is not proven optimal.
* The MILPs of `gapFillMILP`, `fillGaps`, `getMinimalMedium` and `getMinNrFluxes` run on the bundled GLPK and on the GLPK of the COBRA Toolbox, with a `RAVEN:glpkMILP` warning (was: refused).

### Sampling and reporting

* `reporterMetabolites` uses the population standard deviation for the background, sorts the results by Z-score and reports one-sided P-values, which gives the numbers of raven-toolbox (was: scores differed). The score is 0 when the spread of the background is 0 or at rounding level, P-values outside [0, 1] are refused, and metabolites without a score are dropped from all result fields.
* `sampleChebyshevCenter` and `randomSampling` pass a sparse matrix, so the method `'chrr'` works with Gurobi.
* `sampleACHR` projects the directions onto the null space of S and keeps the steps within the bounds. It warns when the samples keep a residual (was: the residual of `S*v` grew exponentially, up to 9 after 100 steps).
* `haveFlux` returns all reactions as blocked when the simplification removes every reaction (was: a GLPK error on an empty matrix).

### Tasks

* `parseTaskList` reads `false` and `0` as false, trims `;` lists, ignores the case of headers, keeps empty CHANGED bounds as `NaN`, and has separate spare columns, so that bounds default to 0 and 1000. Comment rows and orphan rows are skipped with a warning, a non-numeric bound raises `RAVEN:badInput`, and a binary `.xls` file is refused with `RAVEN:badInput` (save it as `.xlsx` or `.txt`).
* `checkTasks` warns when a task metabolite matches several model metabolites. `removeLowScoreGenes` keeps the first gene of the rule on a tie.
* `ftINIT` removes the genes of removed reactions from the output.

### Documentation

* Help blocks of 176 files are rewritten in impersonal, present-tense wording with British spelling; the code is unchanged. Statements that differed from the code are corrected, for example `simplifyModel` `deleteDuplicates` (first occurrence kept), `removeBadRxns` `refModel`, `replaceMets` outputs, `getMinNrFluxes` `params`, `checkRaven` `developMode` and the `downloadRavenBinaries` offline bundle. `reverseRxns` has the RAVEN 3 help structure.

## 3.0.0b1 — 2026-09-07

First beta toward 3.0. Changes below are relative to RAVEN 2.9.3. The last 2.x release is 2.11.4:
a reorganized codebase, fewer toolbox dependencies, hardened I/O and MILP algorithms, plus
capabilities backported from the `raven-toolbox` (Python) port.

### Repository & tooling

* Repo reorganized from a flat layout into functional folders (`analysis/`, `annotation/`,
  `biomass/`, `comparison/`, `conversion/`, `gapfilling/`, `io/`, `manipulation/`, `queries/`,
  `solver/`, `utils/`, etc.).
* Function help blocks reformatted to NumPy-style docstrings; design-rationale, version-history
  and cross-language comments trimmed to current-state prose.
* Allow for positional-or-named function arguments via `parseRAVENargs`. Essential arguments
  are positional, non-compulsory arguments are named.
* Dropped the Statistics, Bioinformatics, Text Analytics and Optimization Toolbox dependencies.
* All 317 `dispEM` call sites across 67 files replaced with native `warning()`/`error()`; new
  `ravenList` for formatted item lists (fixed: was silently dropping its own indent).
* `checkRaven` replaces `checkInstallation`, moved to the repo root, faster checks.
* Universal progress reporting (`progressReport`).
* Poor performing visualization functions removed.

### Reconstruction (KEGG & homology)

* KEGG HMM libraries and external binaries moved out of git to on-demand, SHA256-verified fetch
  from the `raven-data` release repo.
* KEGG artifact set updated 116 → 118.
* Legacy KEGG artifact-builder functions removed, moved to raven-toolbox (Python).
* Recalibration of default HMM cutoffs.
* MetaCyc reconstruction removed, as it persistently yielded high false positive annotations.
* `getModelFromHomology`: recalibration of `minLen` default, token-based GPR rewrite (replaces
  a buggy cleanup regex).
* New `assignCompartments`: functionality-constrained compartment-assignment MILP (WIP).

### Curation

* Genome data utilities for homology-based reconstruction: `getGeneData` (gained a `UniProt`
  output column from the GFF `Dbxref` field; fixed crash when a CDS's `Parent` resolves to no
  gene), `downloadGenomeData`, `processProteinFastaFile`.
* `curateModelFromTables`: fixed missing metabolite-support check.

### I/O

* Legacy pre-L3V1/FBCv2 SBML import removed, only L3V1/FBCv2 remains supported.
* `importExcelModel` removed, to discourage the practice of curating models in Excel.
* New `model.ec` export to Excel (ENZYMES/ENZRXNS sheets).
* Apache POI Excel I/O replaced with a dependency-free OOXML writer.
* Changed YAML file format, to better align with cobrapy and have full compatibility with
  raven-toolbox (Python): fixed quoting; multi-line entries; position-dependent defaults;
  `objective_coefficient` import; empty/missing gene sections; compartment annotations;
  EC-code annotation round-trip.
* Support for `rxnDeltaG`/`metDeltaG` in SBML.

### Manipulation

* `addRxns`: new string keyword for `metsBy`; fixed gene ids shredded by `and`/`or`.
* `addTransport`: fixed row-oriented `metNames`; crash on a compartment with no counterpart
  metabolite; `rxnConfidenceScores`/`rxnDeltaG` defaults.
* `mergeCompartments`: fixed inverted unconstrained warning; stale `compMiriams`.
* `mergeModels`: fixed `metFrom` handling; same-batch id-rename collisions.
* `convertToIrrev`: fixed wrong reverse-copy index.
* `setExchangeBounds`: fixed missing multi-exchange metabolites.
* `closeModel`: fixed `b` growth and `metNotes` padding.
* `replaceMets`: fixed three bugs.
* `deleteUnusedGenes` is now a thin wrapper over `removeReactions`.
* `removeLowScoreGenes`: fixed gene field alignment.
* `permuteModel`: fixed compartment remap.
* `simplifyModel`: fixed dropping task constraints during gap-filling.

### Queries

* `checkModelStruct`: struct-array issues output; fixed two bugs; grRules start/end `and`/`or`
  handling.
* `findPotentialErrors` made public, with parse-tree detection; grRules now parsed into a tree,
  also fixing `removeGenes` and `expandModel`.
* `getExchangeRxns`: fixed index return in both branches, and exchange/transport
  misclassification for boundary-metabolite models.
* `getIndexes`: hash-map lookup, `islogical` fix.
* `getElementalBalance`: empty-reaction fix, new charge-balance check.
* `constructS`/`parseFormulas`: fixed output corruption.
* `findDuplicateRxns` rewritten with `unique(…,'rows')`: O(n log n) instead of O(n²).

### Gap-filling

* New LP/MILP/topological gap-filling algorithms: `gapFillFastCore`, `gapFillSwiftCore`,
 `gapFillMILP` and `gapFillTopological`, in addition to the original `fillGaps` function.
* `makeSomething`/`consumeSomething` merged into `findLeakMetabolite(model, direction, …)`
  (originals kept as thin wrappers).

### Context-specific models (ftINIT / tINIT)

* `getINITModel`/`runINIT` deprecated, in favor of `ftINIT` (a separate, faster
  reimplementation that does not call `runINIT`). Code is still available for now, as
  published Human-GEM-derived models were built with tINIT.
* `ftINIT` gained per-reaction big-M, essential-reaction clamping, tree-based gene removal;
  `runINIT` gained per-reaction big-M.
* ftINIT gained `resolveTies`/`proveAbsGap` for deterministic MILP tie-breaking.
* `ftINIT`: fixed metabolomics `spdiags` crash; MILP handling under the free solver; gap-fill
  failures reported as success; inverted `allowExcretion` sense in `ftINITInternalAlg`.
* `ftINIT` task gap-filling now edits the reference model in place instead of copying it.
* `getMinNrFluxes`: fixed crash on all-non-negative scores.
* `sortReactionOrder`: fixed scoring the wrong columns.
* `checkTasks`: task loop parallelized via `parfor`.

### Analysis & flux sampling

* FSEOF regression-based target selection with knockdown/knockout classification.
* `reporterMetabolites` exact closed-form background correction.
* `getFluxZ`: fixed inverted sign in the zero-variance branch.
* New analysis utilities: `modelSummary`, `walkFluxes`, `traceFluxPath`, `compareFluxes`.
* New `getMinimalMedium`: MILP-based minimal medium finder.
* `randomSampling`: FVA-based `goodRxns` selection, reproducible seed parameter; fixed
  undefined variable; per-reaction LPs parallelized via `parfor`; CHRR and ACHR flux samplers
  added.
* `compareRxnsGenesMetsComps` folded into `compareMultipleModels`.
* New `diffModels`: structured model comparison.

### Biomass

* `fitParameters`: fixed every `xRxn`'s upper bound set from the first.
* `guessComposition`: fixed unknown metabolite's own coefficient.

### Localization

* WoLF PSORT localization loader removed, replaced by DeepLoc2, MULocDeep, COMPARTMENTS,
  UniProt (WIP).

### Solver

* `optimizeProb`: fixed dead no-op block; SCIP translation; gurobi results; MILP guard.
