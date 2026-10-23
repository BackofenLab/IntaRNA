import { FIELDS, PARAMETER_DEFAULTS, EXTRA_OUTPUTS, ENSEMBLE_COLUMNS } from './options.js';
import { CSV_CATALOG } from './columns.js';

// Bound total work as well as per-strand length, including all-vs-all FASTA jobs.
export const LIMITS = Object.freeze({ minLength: 2, maxLength: 2000, maxRecords: 20, maxPairs: 1000000 });
export const CSV_COLUMNS = 'id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E';
export const CONTEXT_COLUMNS = ['id1', 'start1', 'end1', 'id2', 'start2', 'end2', 'E'];
export const DEFAULTS = Object.freeze({ query: '>query\nGGGAGGAGGG', target: '>target\nCCCUCCUCCC', task: 'prediction' });
export const EVALUATION_EXAMPLE = '1||||||||||&1||||||||||';

export function parseFasta(text, label, multiple = true) {
  if (typeof text !== 'string' || text.length > 60000 || !text.trim()) throw new Error(`${label}: enter a sequence or FASTA records.`);
  const lines = text.trim().split(/\r?\n/).map((line) => line.trim()).filter(Boolean);
  const records = [];
  let current;
  for (const line of lines) {
    if (line.startsWith('>')) {
      current = { id: line.slice(1).trim(), sequence: '' };
      if (!current.id || current.id.length > 200 || /[\x00-\x1f\x7f]/.test(current.id)) throw new Error(`${label}: the FASTA header needs a name (at most 200 characters).`);
      // The native FASTA parser takes the identifier before the first space.
      current.id = current.id.split(/\s+/)[0];
      records.push(current);
    } else {
      if (!current) { current = { id: label.toLowerCase(), sequence: '' }; records.push(current); }
      current.sequence += line.replace(/\s/g, '').toUpperCase().replaceAll('T', 'U');
    }
  }
  if ((!multiple && records.length !== 1) || records.length > LIMITS.maxRecords) throw new Error(`${label}: ${multiple ? 'use at most 20 FASTA records' : 'evaluation requires exactly one FASTA record'}.`);
  const ids = new Set();
  for (const record of records) {
    if (!record.sequence || /[^ACGUN]/.test(record.sequence)) throw new Error(`${label}: use only A, C, G, U, T or N in every sequence.`);
    if (record.sequence.length < LIMITS.minLength || record.sequence.length > LIMITS.maxLength) throw new Error(`${label}: use 2–2,000 nucleotides per record.`);
    // Native CSV replaces semicolons with underscores. Keep pair grouping unique.
    const id = record.id.replaceAll(';', '_');
    if (ids.has(id)) throw new Error(`${label}: FASTA identifiers must be unique (also after replacing semicolons with underscores).`);
    ids.add(id);
  }
  return { records, fasta: records.map(({ id, sequence }) => `>${id}\n${sequence}\n`).join('') };
}
export function parseSequence(text, label) {
  const parsed = parseFasta(text, label, false);
  return { ...parsed.records[0], fasta: parsed.fasta };
}

// IntaRNA omits coordinate zero when an index range crosses from negative to positive.
export function sequenceOffset(index, start) {
  return index - start - (start < 0 && index > 0 ? 1 : 0);
}
function ranges(value, key, length, start = 1) {
  if (!value) return [];
  return value.split(',').map((part) => {
    const match = /^(-?[1-9]\d*)-(-?[1-9]\d*)$/.exec(part);
    if (!match) throw new Error(`${key}: use inclusive ranges such as 1-5,8-10 (coordinate zero is skipped).`);
    const from = sequenceOffset(Number(match[1]), start), to = sequenceOffset(Number(match[2]), start);
    if (from < 0 || to < from || to >= length) throw new Error(`${key}: ranges must be ascending and within the selected sequence.`);
    return [from, to];
  });
}
function selectedRecords(parsed, subset, key) {
  const intervals = ranges(subset, key, parsed.records.length);
  return intervals.length ? parsed.records.filter((_, index) => intervals.some(([a, b]) => a <= index && index <= b)) : parsed.records;
}
function oneRecord(records, key) {
  if (records.length !== 1) throw new Error(`${key}: this setting requires exactly one selected sequence.`);
  return records[0];
}
function validateConstraints(value, records, start, key) {
  const record = oneRecord(records, key);
  if (/^[.xpb]+$/.test(value)) {
    if (value.length !== record.sequence.length) throw new Error(`${key}: provide one constraint character per nucleotide.`);
    return;
  }
  const sections = value.match(/[xpb]:[^xpb]+/g);
  if (!sections || sections.join('') !== value) throw new Error(`${key}: use ./x/p/b per base or ranges such as b:3-4,p:6-7.`);
  for (const section of sections) ranges(section.slice(2).replace(/,$/, ''), key, record.sequence.length, start);
}
function validateAccessibility(value, record, mode, key) {
  const lines = value.split(/\r?\n/);
  if (!/^#[\w\s]+$/.test(lines[0]) || !/^\s*#i.\s+l=1(\s+\d+)*\s*$/.test(lines[1] ?? '')) throw new Error(`${key}: paste RNAplfold text with its comment and length header lines.`);
  const lengths = lines[1].match(/\d+/g).map(Number);
  if (lengths.length < 2 || lengths.some((n, i) => n !== i + 1)) throw new Error(`${key}: accessibility lengths must be consecutive, starting at 1, with at least two columns.`);
  const data = lines.slice(2).filter((line) => line.trim());
  if (data.length < record.sequence.length) throw new Error(`${key}: provide one data row per nucleotide.`);
  for (let index = 0; index < record.sequence.length; index++) {
    const cells = data[index].trim().split(/\s+/);
    const values = cells.slice(1, 1 + Math.min(index + 1, lengths.length));
    if (Number(cells[0]) !== index + 1 || values.length !== Math.min(index + 1, lengths.length) || values.some((v) => !Number.isFinite(Number(v)) || Number(v) < 0 || mode === 'P' && Number(v) > 1)) throw new Error(`${key}: rows must be consecutive; required values must be finite ${mode === 'P' ? 'probabilities in [0,1]' : 'opening energies ≥ 0'}.`);
  }
}
function complementary(a, b) { return ['AU', 'UA', 'GC', 'CG', 'GU', 'UG'].includes(a + b); }
export function validateInteractions(encoding, query, target, qStart, tStart, separator = ':') {
  if (typeof encoding !== 'string' || !encoding || encoding.length > 10000 || encoding.split(separator).length > 100) throw new Error('Provide 1–100 target&query dot-bar interactions.');
  for (const entry of encoding.split(separator)) {
    const strands = entry.split('&');
    if (strands.length !== 2) throw new Error('Each interaction needs target&query dot-bar encodings.');
    const positions = strands.map((strand, i) => {
      const match = /^(-?[1-9]\d*)([.|]+)$/.exec(strand);
      if (!match) throw new Error('Use a nonzero start index followed by dots and pairing bars, e.g. 1|||&1|||.');
      if (separator === ',' && (!match[2].startsWith('|') || !match[2].endsWith('|'))) throw new Error('Explicit seed strands must start and end with a pairing bar.');
      const record = i ? query : target, start = i ? qStart : tStart;
      const offset = sequenceOffset(Number(match[1]), start);
      if (offset < 0 || offset + match[2].length > record.sequence.length) throw new Error('A dot-bar interaction extends outside its sequence.');
      return [...match[2]].flatMap((c, p) => c === '|' ? [offset + p] : []);
    });
    if (!positions[0].length || positions[0].length !== positions[1].length) throw new Error('Interaction strands must contain the same nonzero number of pairing bars.');
    for (let p = 0; p < positions[0].length; p++) {
      if (!complementary(target.sequence[positions[0][p]], query.sequence[positions[1].at(-1 - p)])) throw new Error('The supplied interaction contains a non-complementary base pair.');
    }
  }
}

// Validate scalar controls independently so unfinished forms can be saved/loaded.
export function normalizeParameters(request) {
  if (!['prediction', 'evaluation'].includes(request.task ?? 'prediction')) throw new Error('Choose prediction or evaluation.');
  const evaluating = request.task === 'evaluation';
  const parameters = {};
  const supplied = request.parameters ?? {};
  for (const field of FIELDS) {
    if (evaluating && field.prediction) continue;
    const value = supplied[field.key] ?? field.value;
    if (field.type === 'checkbox') {
      if (typeof value !== 'boolean') throw new Error(`${field.label}: choose on or off.`);
      parameters[field.key] = value;
    } else if (field.type === 'number') {
      const n = typeof value === 'number' ? value : (typeof value === 'string' && value.trim() ? Number(value) : NaN);
      if (!Number.isFinite(n) || n < field.min || n > field.max || (field.step !== 'any' && !Number.isInteger(n)) || field.exclude?.includes(n)) throw new Error(`${field.label}: use ${field.min}–${field.max}${field.exclude ? ', excluding 1 and 2' : ''}${field.step !== 'any' ? ' (whole numbers)' : ''}.`);
      parameters[field.key] = n;
    } else {
      if (typeof value !== 'string' || value.length > (field.type === 'textarea' ? 1000000 : 10000) || /\x00/.test(value)) throw new Error(`${field.label}: invalid or oversized text input.`);
      if (field.choices && !Object.hasOwn(field.choices, value)) throw new Error(`${field.label}: choose an available value.`);
      parameters[field.key] = value.trim();
    }
  }
  for (const key of ['qIdxPos0', 'tIdxPos0']) {
    const raw = request[key] ?? 1;
    const n = typeof raw === 'number' || typeof raw === 'string' && raw.trim() ? Number(raw) : NaN;
    if (!Number.isInteger(n) || n < -2000000000 || n > 2000000000) throw new Error(`${key}: use an integer between -2,000,000,000 and 2,000,000,000.`);
    parameters[key] = n;
  }
  return parameters;
}

export function prepareRequest(request) {
  const task = request?.task ?? 'prediction';
  if (!['prediction', 'evaluation'].includes(task)) throw new Error('Choose prediction or evaluation.');
  const evaluating = task === 'evaluation';
  const query = parseFasta(request.query, 'Query', !evaluating), target = parseFasta(request.target, 'Target', !evaluating);
  const parameters = normalizeParameters(request);
  const qRecords = selectedRecords(query, parameters.qSet, 'qSet'), tRecords = selectedRecords(target, parameters.tSet, 'tSet');
  if (parameters.outPairwise && qRecords.length !== tRecords.length) throw new Error('Pairwise output requires the same number of selected query and target records.');
  const work = parameters.outPairwise ? qRecords.reduce((sum, q, i) => sum + q.sequence.length * tRecords[i].sequence.length, 0) : qRecords.reduce((sum, q) => sum + q.sequence.length, 0) * tRecords.reduce((sum, t) => sum + t.sequence.length, 0);
  if (work > LIMITS.maxPairs) throw new Error('The total query × target length across selected sequence pairs must be at most 1,000,000.');

  // Validate sequence-dependent settings before loading the engine.
  for (const [p, records] of [['q', qRecords], ['t', tRecords]]) {
    const start = parameters[`${p}IdxPos0`];
    if (parameters[`${p}AccW`] && parameters[`${p}AccL`] > parameters[`${p}AccW`]) throw new Error(`${p}AccL must not exceed ${p}AccW.`);
    if (parameters[`${p}RegionLenMax`] && (parameters[`${p}Acc`] === 'N' || parameters[`${p}Region`])) throw new Error(`${p}RegionLenMax: automatic regions require accessibility and cannot be combined with manual regions.`);
    for (const record of records) ranges(parameters[`${p}Region`], `${p}Region`, record.sequence.length, start);
    const seedRange = parameters[`seed${p.toUpperCase()}Range`];
    if (seedRange) ranges(seedRange, `seed${p.toUpperCase()}Range`, oneRecord(records, 'Seed regions').sequence.length, start);
    if (parameters[`${p}AccConstr`]) validateConstraints(parameters[`${p}AccConstr`], records, start, `${p}AccConstr`);
    if (['P', 'E'].includes(parameters[`${p}Acc`])) {
      const record = oneRecord(records, 'Accessibility import');
      if (!parameters[`${p}AccFile`]) throw new Error(`${p}Acc: paste accessibility data for the selected import mode.`);
      validateAccessibility(parameters[`${p}AccFile`], record, parameters[`${p}Acc`], `${p}AccFile`);
    }
    if (parameters[`${p}Shape`]) {
      const record = oneRecord(records, 'SHAPE');
      for (const line of parameters[`${p}Shape`].split(/\r?\n/).filter((line) => line.trim() && !line.trim().startsWith('#'))) {
        const cells = line.trim().split(/\s+/);
        if (cells.length !== 2 || !/^\d+$/.test(cells[0]) || Number(cells[0]) < 1 || Number(cells[0]) > record.sequence.length || !Number.isFinite(Number(cells[1]))) throw new Error(`${p}Shape: use position (1-based) and finite reactivity on each line.`);
      }
    }
    if (!/^(D(m-?\d+(\.\d+)?)?(b-?\d+(\.\d+)?)?|Z(b\d+(\.\d+)?)?|W)$/.test(parameters[`${p}ShapeMethod`])) throw new Error(`${p}ShapeMethod: invalid SHAPE method encoding.`);
    if (!/^(M|C(\d+(\.\d+)?)?|S|[LO](s-?\d+(\.\d+)?i-?\d+(\.\d+)?)?)$/.test(parameters[`${p}ShapeConversion`])) throw new Error(`${p}ShapeConversion: invalid SHAPE conversion encoding.`);
  }
  const csvColumns = request.csvColumns ?? CSV_COLUMNS.split(',');
  if (!Array.isArray(csvColumns) || !csvColumns.length || new Set(csvColumns).size !== csvColumns.length || csvColumns.some((key) => !CSV_CATALOG.some((c) => c.key === key))) throw new Error('Select at least one available CSV column, without duplicates.');
  if (parameters.outCsvSort && !csvColumns.includes(parameters.outCsvSort)) throw new Error('The CSV sort column must be selected for output.');
  const outputs = request.outputs ?? [];
  if (!Array.isArray(outputs) || new Set(outputs).size !== outputs.length || outputs.some((key) => !EXTRA_OUTPUTS.some((entry) => entry.key === key && (!evaluating || !entry.prediction)))) throw new Error('Choose available additional output formats for this task.');
  if (outputs.includes('ensemble') && (qRecords.length !== 1 || tRecords.length !== 1)) throw new Error('Ensemble reports require one query and one target.');
  const needsEnsemble = csvColumns.some((c) => ENSEMBLE_COLUMNS.includes(c)) || outputs.some((c) => ['ensemble', 'spotProb', 'qSpotProb', 'tSpotProb'].includes(c));
  if (evaluating) {
    validateInteractions(request.rri, query.records[0], target.records[0], parameters.qIdxPos0, parameters.tIdxPos0);
  } else {
    if (parameters.model === 'B' && parameters.mode !== 'H') throw new Error('Helix-block model B supports heuristic mode only.');
    if (parameters.noSeed && ['S', 'K'].includes(parameters.mode)) throw new Error('Seed-only and greedy prediction require a seed.');
    if (parameters.mode === 'K' && (parameters.model !== 'X' || needsEnsemble)) throw new Error('Greedy prediction requires model X and cannot produce ensemble/probability output.');
    if (parameters.seedTQ) validateInteractions(parameters.seedTQ, oneRecord(qRecords, 'Explicit seeds'), oneRecord(tRecords, 'Explicit seeds'), parameters.qIdxPos0, parameters.tIdxPos0, ',');
    const minimumSeedLength = parameters.seedTQ ? Math.min(...parameters.seedTQ.split(',').map((seed) => seed.split('&')[0].split('|').length - 1)) : parameters.seedBP;
    if (!parameters.noSeed && ['q', 't'].some((p) => parameters[`${p}IntLenMax`] && parameters[`${p}IntLenMax`] < minimumSeedLength)) throw new Error('Maximum interaction lengths must cover the required seed base pairs.');
    if (!parameters.noSeed && !parameters.seedTQ && [...qRecords, ...tRecords].some((r) => r.sequence.length < parameters.seedBP)) throw new Error('Every selected sequence must be at least as long as the seed; reduce Seed base pairs or disable the seed.');
    if (parameters.model === 'B' && (parameters.helixMinBP > parameters.helixMaxBP || !parameters.noSeed && minimumSeedLength > parameters.helixMaxBP)) throw new Error('Helix maximum must cover the helix minimum and required seed base pairs.');
    if (outputs.some((key) => ['bpProb', 'bpsvg'].includes(key)) && (parameters.model !== 'P' || !['H', 'M'].includes(parameters.mode) || parameters.windowWidth || parameters.mode === 'M' && !parameters.noSeed && (parameters.seedTQ || parameters.seedMaxUP || parameters.seedQMaxUP > 0 || parameters.seedTMaxUP > 0))) throw new Error('Base-pair probabilities require model P, heuristic/exact mode and no windows; exact mode supports no seed or stack-only seeds.');
    if (parameters.windowWidth) {
      const maxLength = Math.max(...['q', 't'].map((p) => parameters[`${p}Acc`] === 'N' ? parameters[`${p}IntLenMax`] : Math.max(parameters[`${p}AccW`], parameters[`${p}IntLenMax`])));
      if (parameters.windowWidth < 10 || parameters.windowOverlap >= parameters.windowWidth || parameters.windowOverlap < maxLength || ['q', 't'].some((p) => !parameters[`${p}IntLenMax`] && (parameters[`${p}Acc`] === 'N' || !parameters[`${p}AccW`])) || parameters.outNumber > 1 && parameters.outOverlap !== 'B' || needsEnsemble || parameters.model === 'P') throw new Error('Prediction windows require width ≥ 10, overlap below width and at least the bounded interaction lengths, overlapping output, and no ensemble statistics/model P.');
    }
    // Non-overlap constraints in native merged regions require a single region on
    // the relevant strand; automatic region detection is only safe per region.
    if (!parameters.outPerRegion && parameters.outOverlap !== 'B') {
      for (const p of ['q', 't']) {
        if ((parameters.outOverlap === 'N' || parameters.outOverlap === (p === 'q' ? 'T' : 'Q')) && (parameters[`${p}Region`].includes(',') || parameters[`${p}RegionLenMax`])) throw new Error('Non-overlapping output with multiple/automatic regions requires output per region.');
      }
    }
    if (parameters.outSpotProbSpots && !/^(-?[1-9]\d*&-?[1-9]\d*)(,-?[1-9]\d*&-?[1-9]\d*)*$/.test(parameters.outSpotProbSpots)) throw new Error('Probability spots must be target&query coordinate pairs, e.g. 3&7,5&9.');
    for (const spot of parameters.outSpotProbSpots ? parameters.outSpotProbSpots.split(',') : []) {
      const [t, q] = spot.split('&').map(Number);
      for (const [records, coordinate, start] of [[tRecords, t, parameters.tIdxPos0], [qRecords, q, parameters.qIdxPos0]]) {
        const offset = sequenceOffset(coordinate, start);
        if (records.some((record) => offset < 0 || offset >= record.sequence.length)) throw new Error('Probability spots must fall within every selected sequence.');
      }
    }
  }
  return { task, query, target, parameters, csvColumns, outputs, rri: evaluating ? request.rri : '' };
}

export function makeArguments(request, format = 'C') {
  const args = ['--query=/query.fa', '--target=/target.fa', `--outMode=${format}`, `--out=${format === 'C' ? '/results.csv' : '/report.txt'}`];
  if (request.task === 'evaluation') args.push('--personality=IntaRNAeval');
  for (const [key, value] of Object.entries(request.parameters)) {
    if (key === 'outCsvSort' && format !== 'C' || key === 'outSpotProbSpots' || key === 'kineticScore' && request.parameters.mode !== 'K') continue;
    // Omit defaults: the native parser distinguishes explicit values in several
    // cross-option checks, especially accessibility/region constraints.
    if (value === '' || value === PARAMETER_DEFAULTS[key]) continue;
    const fileInput = /^(q|t)(Shape|AccFile)$/.test(key);
    args.push(`--${key}=${fileInput ? `/${key}.txt` : value}`);
  }
  if (request.rri) args.push(`--rri=${request.rri}`);
  if (format === 'C') {
    args.push(`--outCsvCols=${[...new Set([...request.csvColumns, ...CONTEXT_COLUMNS, 'hybridDP'])].join(',')}`);
    for (const key of request.outputs.filter((key) => !['normal', 'detailed', 'ensemble'].includes(key))) {
      const spots = key === 'spotProb' && request.parameters.outSpotProbSpots ? `${request.parameters.outSpotProbSpots}:` : '';
      args.push(`--out=${key}:${spots}/output/${key}${key === 'bpsvg' ? '.svg' : '.csv'}`);
    }
  }
  return args;
}

export const JSON_SCHEMA = 'IntaRNA-parameters-v1';

// One parameter document crosses the JS/native boundary. File contents use the
// native JSON embedding contract, so no browser-specific input paths are needed.
export function makeParameterDocument(request, format = 'C') {
  const parameters = {};
  for (const argument of makeArguments(request, format)) {
    const separator = argument.indexOf('='), key = argument.slice(2, separator);
    let value = argument.slice(separator + 1);
    if (key === 'query' || key === 'target') value = { text: request[key].fasta };
    else if (/^[qt](Shape|AccFile)$/.test(key)) value = { text: request.parameters[key] + '\n' };
    else if (Object.hasOwn(request.parameters, key)) value = request.parameters[key];
    if (key === 'out') (parameters.out ??= []).push(value);
    else if (!(['qIdxPos0', 'tIdxPos0'].includes(key) && value === 1)) parameters[key] = value;
  }
  return { schema: JSON_SCHEMA, parameters };
}
