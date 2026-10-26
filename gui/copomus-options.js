export const COPOMUS_MEASURES = Object.freeze([
  ['mfeCover', 'MFE coverage'], ['E', 'Interaction energy'],
  ['minDeltaE', 'Minimum energy difference'], ['EDqi', 'Query accessibility'],
  ['cEDqi', 'Constant query accessibility'], ['Eqi', 'Query energy profile'],
  ['Eti', 'Target energy profile'],
]);
export const COPOMUS_FILTERS = Object.freeze([
  ['GU', 'Exclude GU/UG pairs'], ['AU', 'Exclude AU/UA pairs'],
  ['GC', 'Exclude GC/CG pairs'], ['lp', 'Exclude unstackable lonely pairs'],
  ['lpMfe', 'Exclude lonely pairs in the MFE'], ['he', 'Exclude potential helix ends'],
  ['heMfe', 'Exclude helix ends in the MFE'],
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
