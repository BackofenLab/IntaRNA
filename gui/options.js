import { CSV_CATALOG } from './columns.js';

// Browser controls mirror CommandLineParsing.cpp; process/file-management options
// are replaced by local text inputs, downloads and a single cancellable worker.
const number = (key, label, value, min, max, extra = {}) => ({ key, label, value, min, max, type: 'number', ...extra });
const select = (key, label, value, choices, extra = {}) => ({ key, label, value, choices, type: 'select', ...extra });
const flag = (key, label, extra = {}) => ({ key, label, value: false, type: 'checkbox', ...extra });
const text = (key, label, extra = {}) => ({ key, label, value: '', type: 'text', ...extra });
const prediction = { prediction: true };
const strandHelp = (strand) => ({
  Set: `Selects which ${strand} FASTA records participate in the computation; an empty value uses every record.`,
  Acc: `Chooses whether ${strand} accessibility is computed, ignored, or imported as unpaired probabilities or opening energies.`,
  AccW: `Sets the sliding window used to compute local ${strand} accessibility; zero folds the complete sequence.`,
  AccL: `Limits the base-pair span considered while computing ${strand} accessibility; zero uses the accessibility-window size.`,
  AccConstr: `Constrains nucleotides to be paired, unpaired or blocked during ${strand} accessibility folding.`,
  AccFile: `Provides precomputed ${strand} unpaired probabilities or opening energies instead of calculating accessibility.`,
  IntLenMax: `Limits how many ${strand} nucleotides a reported interaction may span; zero uses the complete available region.`,
  IntLoopMax: `Limits the number of unpaired ${strand} nucleotides allowed in one interaction interior loop.`,
  Region: `Restricts ${strand} interaction prediction to the listed inclusive sequence-coordinate ranges.`,
  RegionLenMax: `Splits the ${strand} into automatically detected accessible regions no longer than this value; zero disables the split.`,
  PfScale: `Sets the ${strand} partition-function scaling factor used for accessibility computations.`,
  Shape: `Supplies nucleotide-position and reactivity pairs that constrain ${strand} accessibility folding.`,
  ShapeMethod: `Chooses how ${strand} SHAPE pseudo-energies are added to the folding model.`,
  ShapeConversion: `Chooses how ${strand} SHAPE reactivities are converted before pseudo-energy calculation.`,
});
const qHelp = strandHelp('query'), tHelp = strandHelp('target');
const FIELD_HELP = Object.freeze({
  ...Object.fromEntries(Object.entries(qHelp).map(([suffix, help]) => [`q${suffix}`, help])),
  ...Object.fromEntries(Object.entries(tHelp).map(([suffix, help]) => [`t${suffix}`, help])),
  accNoLP: 'Excludes lonely base pairs from the intramolecular structures used to calculate accessibility on both strands.',
  accNoGUend: 'Excludes GU pairs at intramolecular helix ends when calculating accessibility on both strands.',
  noSeed: 'Disables the requirement that every predicted interaction contain a seed region.',
  seedTQ: 'Defines exact target-and-query seed structures instead of identifying seeds automatically.',
  seedBP: 'Sets the minimum number of base pairs required in an automatically identified seed.',
  seedMaxUP: 'Limits the total unpaired nucleotides across both strands within a seed.',
  seedQMaxUP: 'Limits unpaired query nucleotides within a seed; minus one uses the combined seed limit.',
  seedTMaxUP: 'Limits unpaired target nucleotides within a seed.',
  seedMaxE: 'Rejects seeds whose complete energy exceeds this threshold in kcal/mol.',
  seedMaxEhybrid: 'Rejects seeds whose intermolecular hybridization energy exceeds this threshold in kcal/mol.',
  seedMinPu: 'Requires at least this joint probability that the seed nucleotides are unpaired before interaction.',
  seedQRange: 'Restricts query seed positions to the listed inclusive sequence-coordinate ranges.',
  seedTRange: 'Restricts target seed positions to the listed inclusive sequence-coordinate ranges.',
  seedNoGU: 'Disallows GU wobble pairs anywhere inside an automatically identified seed.',
  seedNoGUend: 'Disallows GU wobble pairs at either end of an automatically identified seed.',
  helixMinBP: 'Sets the minimum number of base pairs in each helix block of the helix-block interaction model.',
  helixMaxBP: 'Sets the maximum number of base pairs in each helix block of the helix-block interaction model.',
  helixMaxIL: 'Limits the number of unpaired nucleotides allowed in an interior loop inside a helix block.',
  helixMinPu: 'Requires at least this probability that the nucleotides of a helix block are initially unpaired.',
  helixMaxE: 'Rejects helix blocks whose energy exceeds this threshold in kcal/mol.',
  helixFullE: 'Applies the helix energy limit to the full helix energy instead of its hybridization contribution.',
  model: 'Chooses the thermodynamic interaction model used for prediction.',
  mode: 'Chooses the algorithm used to search the selected interaction model.',
  kineticScore: 'Chooses how candidate extensions are scored by the greedy kinetic seed-extension algorithm.',
  energyVRNA: 'Selects the ViennaRNA nearest-neighbor energy parameter set.',
  temperature: 'Sets the folding and interaction temperature in degrees Celsius.',
  energyAdd: 'Adds this constant energy in kcal/mol to every reported interaction.',
  energyNoDangles: 'Omits dangling-end contributions from interaction-energy calculations.',
  windowWidth: 'Divides long sequences into prediction windows of this width; zero disables prediction windows.',
  windowOverlap: 'Sets how many nucleotides adjacent prediction windows share.',
  outNumber: 'Limits the number of reported interactions per query-target pair; zero calculates accessibility without reporting interactions.',
  outNumberPerSeed: 'Keeps this many best candidate sites per seed before final filtering; zero disables per-seed selection.',
  outOverlap: 'Controls on which strand separately reported interactions may overlap.',
  outMaxE: 'Reports only interactions with energy below this exclusive threshold in kcal/mol.',
  outDeltaE: 'Reports interactions no more than this many kcal/mol above the best interaction for the same sequence pair.',
  outMinPu: 'Requires at least this joint probability that the interaction sites are initially unpaired.',
  outBestSeedOnly: 'Retains only the energetically best seed when an interaction contains multiple valid seeds.',
  outNoLP: 'Excludes interactions containing isolated intermolecular base pairs.',
  outNoGUend: 'Excludes interactions with a GU wobble pair at either interaction-helix end.',
  outPerRegion: 'Reports and filters interaction candidates independently for every selected region combination.',
  outPairwise: 'Pairs query and target FASTA records by position instead of predicting every record combination.',
  outCsvSort: 'Sorts the downloaded CSV rows by the selected native output column.',
  outSpotProbSpots: 'Restricts interaction-spot probability output to the listed target-and-query coordinate pairs.',
});
const strand = (prefix, label) => [
  text(`${prefix}Set`, 'FASTA record subset', { ...prediction, hint: 'Record numbers, e.g. 1-3,5-5. Empty selects all.' }),
  select(`${prefix}Acc`, 'Accessibility', 'C', { C: 'Compute', N: 'Ignore', P: 'Import unpaired probabilities', E: 'Import opening energies' }),
  number(`${prefix}AccW`, 'Accessibility window', 150, 0, 99999, { exclude: [1, 2], hint: '0 uses the full sequence; otherwise at least 3 nt.' }),
  number(`${prefix}AccL`, 'Maximum base-pair span', 100, 0, 99999, { exclude: [1, 2], hint: '0 uses the accessibility window.' }),
  text(`${prefix}AccConstr`, 'Accessibility constraints', { hint: 'One sequence only: ./x/p/b per base, or ranges such as b:3-4,p:6-7.' }),
  text(`${prefix}AccFile`, 'Accessibility data', { type: 'textarea', hint: 'Paste RNAplfold text for one sequence when importing P or E.' }),
  number(`${prefix}IntLenMax`, 'Maximum interaction length', 0, 0, 99999, { ...prediction, hint: '0 uses the full sequence, bounded by the accessibility window.' }),
  number(`${prefix}IntLoopMax`, 'Maximum unpaired bases per interior loop', prefix === 'q' ? 16 : 10, 0, 30, prediction),
  text(`${prefix}Region`, 'Interaction regions', { ...prediction, hint: 'Inclusive sequence coordinates, e.g. 1-20,30-50.' }),
  number(`${prefix}RegionLenMax`, 'Automatic region length limit', 0, 0, 99999, { ...prediction, hint: '0 disables automatic region detection.' }),
  number(`${prefix}PfScale`, `${label} partition-function scale`, 1.07, 1, 99999, { step: 'any' }),
];
const fields = new Map([
  ...strand('q', 'Query'), ...strand('t', 'Target'),
  flag('accNoLP', 'Exclude lonely pairs in accessibility'),
  flag('accNoGUend', 'Exclude GU helix ends in accessibility'),
  flag('noSeed', 'Do not require a seed', prediction),
  text('seedTQ', 'Explicit seeds', { ...prediction, hint: 'Target&query dot-bar encodings, e.g. 1|||&1|||; separate seeds with commas.' }),
  number('seedBP', 'Seed base pairs', 7, 2, 20, prediction),
  number('seedMaxUP', 'Maximum unpaired seed bases (both strands)', 0, 0, 20, prediction),
  number('seedQMaxUP', 'Maximum unpaired query seed bases', -1, -1, 20, { ...prediction, hint: '-1 uses the combined limit.' }),
  number('seedTMaxUP', 'Maximum unpaired target seed bases', -1, -1, 20, prediction),
  number('seedMaxE', 'Maximum seed energy (kcal/mol)', 0, -999, 999, { ...prediction, step: 'any' }),
  number('seedMaxEhybrid', 'Maximum seed hybridization energy (kcal/mol)', 999, -999, 999, { ...prediction, step: 'any' }),
  number('seedMinPu', 'Minimum seed unpaired probability', 0, 0, 1, { ...prediction, step: 'any' }),
  text('seedQRange', 'Query seed regions', prediction), text('seedTRange', 'Target seed regions', prediction),
  flag('seedNoGU', 'Exclude GU pairs in seeds', prediction), flag('seedNoGUend', 'Exclude GU pairs at seed ends', prediction),
  number('helixMinBP', 'Minimum helix base pairs', 2, 2, 4, prediction),
  number('helixMaxBP', 'Maximum helix base pairs', 10, 2, 20, prediction),
  number('helixMaxIL', 'Maximum helix interior-loop size', 0, 0, 2, prediction),
  number('helixMinPu', 'Minimum helix unpaired probability', 0, 0, 1, { ...prediction, step: 'any' }),
  number('helixMaxE', 'Maximum helix energy (kcal/mol)', 0, -999, 999, { ...prediction, step: 'any' }),
  flag('helixFullE', 'Use full helix energy for the energy limit', prediction),
  select('model', 'Interaction model', 'X', { X: 'Seed extension (MFE)', S: 'Single site (MFE)', P: 'Ensemble energy', B: 'Helix blocks (MFE)' }, prediction),
  select('mode', 'Prediction algorithm', 'H', { H: 'Heuristic', M: 'Exact', S: 'Seed only', K: 'Greedy seed extension' }, prediction),
  select('kineticScore', 'Greedy extension score', 'A', { A: 'Energy change', B: 'Change / (1 + s1 + s2)', C: 'Change / (1 + 2 max(s1,s2))' }, prediction),
  select('energyVRNA', 'ViennaRNA parameter set', 'Turner04', { Turner04: 'Turner 2004', Turner99: 'Turner 1999', Andronescu07: 'Andronescu 2007' }),
  number('temperature', 'Temperature (°C)', 37, 0, 100, { step: 'any' }),
  number('energyAdd', 'Additional energy (kcal/mol)', 0, -999, 999, { step: 'any' }),
  flag('energyNoDangles', 'Exclude dangling-end energy'),
  number('windowWidth', 'Prediction window width', 0, 0, 99999, { ...prediction, hint: '0 disables prediction windows; otherwise at least 10 nt.' }),
  number('windowOverlap', 'Prediction window overlap', 150, 0, 99999, prediction),
  number('outNumber', 'Maximum interactions per sequence pair', 1, 0, 100, { ...prediction, hint: 'Browser limit: 100; 0 computes accessibility only.' }),
  number('outNumberPerSeed', 'Best sites retained per seed', 0, 0, 1000, { ...prediction, hint: '0 disables per-seed candidate selection; positive values can increase runtime.' }),
  select('outOverlap', 'Allow overlapping interactions', 'B', { B: 'Both strands', N: 'Neither strand', T: 'Target only', Q: 'Query only' }, prediction),
  number('outMaxE', 'Maximum interaction energy (exclusive, kcal/mol)', 0, -999, 999, { ...prediction, step: 'any' }),
  number('outDeltaE', 'Maximum difference from MFE (kcal/mol)', 100, 0, 100, { ...prediction, step: 'any' }),
  number('outMinPu', 'Minimum unpaired probability', 0, 0, 1, { ...prediction, step: 'any' }),
  flag('outBestSeedOnly', 'Report only the best seed', prediction),
  flag('outNoLP', 'Exclude lonely interaction pairs', prediction), flag('outNoGUend', 'Exclude GU interaction helix ends', prediction),
  flag('outPerRegion', 'Report each region combination separately', prediction), flag('outPairwise', 'Match FASTA records pairwise', prediction),
  select('outCsvSort', 'CSV download sort column', '', { '': 'Native computation order', ...Object.fromEntries(CSV_CATALOG.map((column) => [column.key, column.label])) }),
  text('outSpotProbSpots', 'Interaction probability spots', { ...prediction, hint: 'Optional target&query coordinate pairs, e.g. 3&7,5&9. Empty selects all spots.' }),
  ...['q', 't'].flatMap((p) => [
    text(`${p}Shape`, `${p === 'q' ? 'Query' : 'Target'} SHAPE data`, { type: 'textarea', hint: 'One sequence only. Paste position and reactivity, one pair per line; positions start at 1.' }),
    { ...text(`${p}ShapeMethod`, 'SHAPE method'), value: 'D', hint: 'D, Dm1.8b-0.6, Z, Zb0.89 or W.' },
    { ...text(`${p}ShapeConversion`, 'SHAPE conversion'), value: 'Os1.6i-2.29', hint: 'M, C0.25, S, Ls0.68i0.2 or Os1.6i-2.29.' },
  ]),
].map((field) => [field.key, { ...field, help: FIELD_HELP[field.key] }]));

// Explicit rows separate visual order from native parameter defaults.
const rows = (keys) => keys.map((row) => row.map((key) => fields.get(key)));
const strandRows = (p) => rows([
  [`${p}Set`, `${p}Region`, `${p}IntLenMax`, `${p}IntLoopMax`, `${p}RegionLenMax`],
  [`${p}Acc`, `${p}AccW`, `${p}AccL`, `${p}PfScale`],
]);
const strandGroup = (p, label) => ({
  key: `${label.toLowerCase()}-options`, label, strand: p,
  rows: strandRows(p),
  fields: [...strandRows(p).flat(), ...rows([
    ...(p === 'q' ? [['accNoLP', 'accNoGUend']] : []),
    [`${p}AccConstr`, `${p}AccFile`, `${p}Shape`, `${p}ShapeMethod`, `${p}ShapeConversion`],
  ]).flat()],
});
export const GROUPS = [
  { key: 'seed-options', label: 'Seed', prediction: true, rows: rows([
    ['noSeed'], ['seedTQ'], ['seedBP', 'seedMaxUP', 'seedQMaxUP', 'seedTMaxUP'],
    ['seedMaxE', 'seedMinPu', 'seedMaxEhybrid'], ['seedQRange', 'seedTRange'], ['seedNoGU', 'seedNoGUend'],
  ]) },
  { key: 'output-options', label: 'Output', rows: rows([
    ['outNumber', 'outNumberPerSeed', 'outOverlap'],
    ['outMaxE', 'outDeltaE', 'outMinPu'],
    ['outBestSeedOnly', 'outNoLP', 'outNoGUend', 'outPerRegion', 'outPairwise'],
    ['outSpotProbSpots'],
  ]) },
  { key: 'csv-options', label: 'CSV', fields: [fields.get('outCsvSort')] },
  { key: 'interaction-options', label: 'Model & Energy', rows: rows([
    ['model', 'mode'], ['energyVRNA', 'temperature', 'energyAdd'], ['energyNoDangles'],
    ['windowWidth', 'windowOverlap'], ['kineticScore'],
    ['helixMinBP', 'helixMaxBP', 'helixMaxIL', 'helixMinPu', 'helixMaxE'], ['helixFullE'],
  ]) },
  strandGroup('q', 'Query'), strandGroup('t', 'Target'),
].map((group) => ({ ...group, fields: group.fields ?? group.rows.flat() }));
export const FIELDS = GROUPS.flatMap((group) => group.fields);
export const PARAMETER_DEFAULTS = Object.freeze(Object.fromEntries(FIELDS.map((field) => [field.key, field.value])));
export const UI_DEFAULTS = Object.freeze({ ...PARAMETER_DEFAULTS, outNumber: 10 });
export const EXTRA_OUTPUTS = [
  ['qMinE', 'Query minimum energies', true, 'CSV position profile of the lowest interaction energy covering each query nucleotide. Useful for locating favorable query sites.'],
  ['tMinE', 'Target minimum energies', true, 'CSV position profile of the lowest interaction energy covering each target nucleotide. Useful for locating favorable target sites.'],
  ['qAcc', 'Query opening energies', false, 'RNAplfold-style opening-energy matrix for the query. Use it to inspect or reproduce query accessibility.'],
  ['tAcc', 'Target opening energies', false, 'RNAplfold-style opening-energy matrix for the target. Use it to inspect or reproduce target accessibility.'],
  ['qPu', 'Query unpaired probabilities', false, 'RNAplfold-style query unpaired-probability matrix, converted from rounded opening energies.'],
  ['tPu', 'Target unpaired probabilities', false, 'RNAplfold-style target unpaired-probability matrix, converted from rounded opening energies.'],
  ['qSpotProb', 'Query interaction probabilities', true, 'CSV position profile estimating the probability that each query nucleotide is covered by an interaction.'],
  ['tSpotProb', 'Target interaction probabilities', true, 'CSV position profile estimating the probability that each target nucleotide is covered by an interaction.'],
  ['pMinE', 'Pair minimum energies', true, 'CSV target-by-query matrix of the minimum interaction energy covering each nucleotide pair. Useful as an energy heatmap.'],
  ['spotProb', 'Interaction spot probabilities', true, 'CSV target-by-query matrix, or selected spots, estimating interaction coverage probabilities. Useful as a probability heatmap.'],
  ['bpProb', 'Actual base-pair probabilities (model P)', true, 'CSV target-by-query matrix of actual base-pair probabilities for the model P interaction ensemble.'],
  ['bpsvg', 'Base-pair probability SVG (model P)', true, 'SVG dot plot of model P base-pair probabilities with accessibility frames and seed/MFE outlines.'],
  ['normal', 'ASCII interaction structures', false, 'Compact text report with ASCII interaction structures and their overall energies.'],
  ['detailed', 'Detailed interaction report', false, 'Expanded text report with interaction structures, coordinates and energy contributions.'],
  ['ensemble', 'Ensemble report', false, 'Text report of interaction and intramolecular ensemble energies and partition functions.'],
].map(([key, label, predictionOnly, description]) => ({ key, label, prediction: predictionOnly, description }));
export const ENSEMBLE_COLUMNS = ['Eall', 'Zall', 'EallTotal', 'P_E'];
