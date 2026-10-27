export const COPOMUS_MEASURES = Object.freeze([
  ['mfeCover', 'MFE coverage', 'Prefers mutations whose designed positions lie inside both the wild-type and compensatory minimum-free-energy interactions.'],
  ['E', 'Interaction energy', 'Checks that both matching pairs interact favorably while each single mutant weakens the interaction by the configured separation threshold.'],
  ['minDeltaE', 'Minimum energy difference', 'Ranks by the smallest energy separation between a matching pair and either single-mutant combination.'],
  ['EDqi', 'Query accessibility', 'Prefers mutations whose compensatory query nucleotide is less accessible than its wild-type counterpart.'],
  ['cEDqi', 'Constant query accessibility', 'Prefers mutations that preserve the query nucleotide accessibility within one kcal/mol.'],
  ['Eqi', 'Query energy profile', 'Checks the interaction-energy profile at the mutated query position across all four sequence combinations.'],
  ['Eti', 'Target energy profile', 'Checks the interaction-energy profile at the mutated target position across all four sequence combinations.'],
]);
export const COPOMUS_FILTERS = Object.freeze([
  ['GU', 'Exclude GU/UG pairs', 'Removes candidate positions whose wild-type base pair is GU or UG.'],
  ['AU', 'Exclude AU/UA pairs', 'Removes candidate positions whose wild-type base pair is AU or UA.'],
  ['GC', 'Exclude GC/CG pairs', 'Removes candidate positions whose wild-type base pair is GC or CG.'],
  ['lp', 'Exclude unstackable lonely pairs', 'Removes candidates that cannot stack with a neighboring canonical base pair in the input sequences.'],
  ['lpMfe', 'Exclude lonely pairs in the MFE', 'Removes candidates that form a lonely base pair in the wild-type minimum-free-energy interaction.'],
  ['he', 'Exclude potential helix ends', 'Removes candidates that could occur at either end of an interaction helix.'],
  ['heMfe', 'Exclude helix ends in the MFE', 'Removes candidates located at a helix end in the wild-type minimum-free-energy interaction.'],
]);
export const COPOMUS_DEFAULTS = Object.freeze({
  measures: ['mfeCover', 'E', 'minDeltaE'],
  candidateSelection: 'mfe', candidateFilters: [], generator: 'flip',
  mutationEncoding: '', alpha: 1, beta: 1,
});

const choices = (catalog) => new Set(catalog.map(([key]) => key));
const measureChoices = choices(COPOMUS_MEASURES), filterChoices = choices(COPOMUS_FILTERS);
const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
const sequenceOffset = (index, start) => index - start - (start < 0 && index > 0 ? 1 : 0);

export function normalizeCopomus(options = {}, context) {
  if (!isObject(options)) throw new Error('Mutation Design settings must be an object.');
  const normalized = { ...COPOMUS_DEFAULTS, ...options };
  for (const [key, catalog, label, nonempty] of [
    ['measures', measureChoices, 'measures', true],
    ['candidateFilters', filterChoices, 'candidate filters', false],
  ]) {
    const values = normalized[key];
    if (!Array.isArray(values) || nonempty && !values.length || new Set(values).size !== values.length || values.some((value) => !catalog.has(value))) throw new Error(`Choose available CopomuS ${label}, without duplicates.`);
    normalized[key] = [...values];
  }
  if (!['mfe', 'mfeSO'].includes(normalized.candidateSelection)) throw new Error('Choose an available CopomuS candidate-selection method.');
  if (!['flip', 'any'].includes(normalized.generator)) throw new Error('Choose an available CopomuS mutation generator.');
  if (typeof normalized.mutationEncoding !== 'string' || normalized.mutationEncoding.length > 10000 || normalized.mutationEncoding.includes('\0')) throw new Error('Mutation encoding must be at most 10,000 text characters.');
  normalized.mutationEncoding = normalized.mutationEncoding.trim().toUpperCase().replaceAll(' ', '');
  for (const key of ['alpha', 'beta']) {
    const value = typeof normalized[key] === 'number' ? normalized[key] : typeof normalized[key] === 'string' && normalized[key].trim() ? Number(normalized[key]) : NaN;
    if (!Number.isFinite(value)) throw new Error(`${key}: enter a finite number.`);
    normalized[key] = value;
  }
  if (!context) return normalized;
  const { query, target, qIdxPos0, tIdxPos0 } = context;
  if (query.sequence.includes('N') || target.sequence.includes('N')) throw new Error('Mutation Design requires unambiguous A, C, G or U sequences.');
  if (!normalized.mutationEncoding) return normalized;
  for (const entry of normalized.mutationEncoding.split(',')) {
    const match = /^([ACGU])(-?\d+)([ACGU])&([ACGU])(-?\d+)([ACGU])$/.exec(entry);
    if (!match) throw new Error(`Invalid mutation encoding: “${entry}”. Use query-first substitutions such as G1C&C10G.`);
    const [, queryWild, queryPosition, queryMutated, targetWild, targetPosition, targetMutated] = match;
    const queryIndex = sequenceOffset(Number(queryPosition), qIdxPos0), targetIndex = sequenceOffset(Number(targetPosition), tIdxPos0);
    if (Number(queryPosition) === 0 || Number(targetPosition) === 0 || queryIndex < 0 || queryIndex >= query.sequence.length || targetIndex < 0 || targetIndex >= target.sequence.length) throw new Error('A mutation position is outside its sequence or uses skipped coordinate zero.');
    if (query.sequence[queryIndex] !== queryWild || target.sequence[targetIndex] !== targetWild) throw new Error('Mutation encoding does not match the wildtype sequences.');
    if (queryWild === queryMutated || targetWild === targetMutated) throw new Error('Both encoded bases must change.');
  }
  return normalized;
}
