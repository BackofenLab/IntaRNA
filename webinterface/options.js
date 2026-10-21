import { CSV_CATALOG } from './columns.js';

// Browser controls mirror CommandLineParsing.cpp; process/file-management options
// are replaced by local text inputs, downloads and a single cancellable worker.
const number = (key, label, value, min, max, extra = {}) => ({ key, label, value, min, max, type: 'number', ...extra });
const select = (key, label, value, choices, extra = {}) => ({ key, label, value, choices, type: 'select', ...extra });
const flag = (key, label, extra = {}) => ({ key, label, value: false, type: 'checkbox', ...extra });
const text = (key, label, extra = {}) => ({ key, label, value: '', type: 'text', ...extra });
const prediction = { prediction: true };
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
export const GROUPS = [
  { key: 'query-options', label: 'Query', fields: strand('q', 'Query') },
  { key: 'target-options', label: 'Target', fields: strand('t', 'Target') },
  { key: 'seed-options', label: 'Seed', prediction: true, fields: [
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
  ] },
  { key: 'helix-options', label: 'Helix (model B)', prediction: true, fields: [
    number('helixMinBP', 'Minimum helix base pairs', 2, 2, 4, prediction),
    number('helixMaxBP', 'Maximum helix base pairs', 10, 2, 20, prediction),
    number('helixMaxIL', 'Maximum helix interior-loop size', 0, 0, 2, prediction),
    number('helixMinPu', 'Minimum helix unpaired probability', 0, 0, 1, { ...prediction, step: 'any' }),
    number('helixMaxE', 'Maximum helix energy (kcal/mol)', 0, -999, 999, { ...prediction, step: 'any' }),
    flag('helixFullE', 'Use full helix energy for the energy limit', prediction),
  ] },
  { key: 'shape-options', label: 'SHAPE', fields: ['q', 't'].flatMap((p) => [
    text(`${p}Shape`, `${p === 'q' ? 'Query' : 'Target'} SHAPE data`, { type: 'textarea', hint: 'One sequence only. Paste position and reactivity, one pair per line; positions start at 1.' }),
    { ...text(`${p}ShapeMethod`, 'SHAPE method'), value: 'D', hint: 'D, Dm1.8b-0.6, Z, Zb0.89 or W.' },
    { ...text(`${p}ShapeConversion`, 'SHAPE conversion'), value: 'Os1.6i-2.29', hint: 'M, C0.25, S, Ls0.68i0.2 or Os1.6i-2.29.' },
  ]) },
  { key: 'interaction-options', label: 'Interaction and energy', fields: [
    select('model', 'Interaction model', 'X', { X: 'Seed extension (MFE)', S: 'Single site (MFE)', P: 'Ensemble energy', B: 'Helix blocks (MFE)' }, prediction),
    select('mode', 'Prediction algorithm', 'H', { H: 'Heuristic', M: 'Exact', S: 'Seed only', K: 'Greedy seed extension' }, prediction),
    select('kineticScore', 'Greedy extension score', 'A', { A: 'Energy change', B: 'Change / (1 + s1 + s2)', C: 'Change / (1 + 2 max(s1,s2))' }, prediction),
    select('energyVRNA', 'ViennaRNA parameter set', 'Turner04', { Turner04: 'Turner 2004', Turner99: 'Turner 1999', Andronescu07: 'Andronescu 2007' }),
    number('temperature', 'Temperature (°C)', 37, 0, 100, { step: 'any' }),
    number('energyAdd', 'Additional energy (kcal/mol)', 0, -999, 999, { step: 'any' }),
    flag('energyNoDangles', 'Exclude dangling-end energy'), flag('accNoLP', 'Exclude lonely pairs in accessibility'),
    flag('accNoGUend', 'Exclude GU helix ends in accessibility'),
    number('windowWidth', 'Prediction window width', 0, 0, 99999, { ...prediction, hint: '0 disables prediction windows; otherwise at least 10 nt.' }),
    number('windowOverlap', 'Prediction window overlap', 150, 0, 99999, prediction),
  ] },
  { key: 'output-options', label: 'Output', fields: [
    number('outNumber', 'Maximum interactions per sequence pair', 1, 0, 100, { ...prediction, hint: 'Browser limit: 100; 0 computes accessibility only.' }),
    select('outOverlap', 'Allow overlapping interactions', 'B', { B: 'Both strands', N: 'Neither strand', T: 'Target only', Q: 'Query only' }, prediction),
    number('outMaxE', 'Maximum interaction energy (exclusive, kcal/mol)', 0, -999, 999, { ...prediction, step: 'any' }),
    number('outDeltaE', 'Maximum difference from MFE (kcal/mol)', 100, 0, 100, { ...prediction, step: 'any' }),
    number('outMinPu', 'Minimum unpaired probability', 0, 0, 1, { ...prediction, step: 'any' }),
    flag('outBestSeedOnly', 'Report only the best seed', prediction),
    flag('outNoLP', 'Exclude lonely interaction pairs', prediction), flag('outNoGUend', 'Exclude GU interaction helix ends', prediction),
    flag('outPerRegion', 'Report each region combination separately', prediction), flag('outPairwise', 'Match FASTA records pairwise', prediction),
    select('outCsvSort', 'CSV download sort column', '', { '': 'Native computation order', ...Object.fromEntries(CSV_CATALOG.map((column) => [column.key, column.label])) }),
    text('outSpotProbSpots', 'Interaction probability spots', { ...prediction, hint: 'Optional target&query coordinate pairs, e.g. 3&7,5&9. Empty selects all spots.' }),
  ] },
];
export const FIELDS = GROUPS.flatMap((group) => group.fields);
export const PARAMETER_DEFAULTS = Object.freeze(Object.fromEntries(FIELDS.map((field) => [field.key, field.value])));
export const EXTRA_OUTPUTS = [
  ['qMinE', 'Query minimum energies', true], ['tMinE', 'Target minimum energies', true],
  ['qAcc', 'Query opening energies'], ['tAcc', 'Target opening energies'],
  ['qPu', 'Query unpaired probabilities'], ['tPu', 'Target unpaired probabilities'],
  ['qSpotProb', 'Query interaction probabilities', true], ['tSpotProb', 'Target interaction probabilities', true],
  ['pMinE', 'Pair minimum energies', true], ['spotProb', 'Interaction spot probabilities', true],
  ['bpProb', 'Actual base-pair probabilities (model P)', true], ['bpsvg', 'Base-pair probability SVG (model P)', true],
  ['normal', 'ASCII interaction structures'], ['detailed', 'Detailed interaction report'], ['ensemble', 'Ensemble report'],
].map(([key, label, predictionOnly = false]) => ({ key, label, prediction: predictionOnly }));
export const ENSEMBLE_COLUMNS = ['Eall', 'Zall', 'EallTotal', 'P_E'];
