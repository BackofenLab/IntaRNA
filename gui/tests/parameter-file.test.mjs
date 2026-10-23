import { test } from 'node:test';
import assert from 'node:assert/strict';
import { DEFAULTS, EVALUATION_EXAMPLE, prepareRequest, CSV_COLUMNS } from '../protocol.js';
import { FIELDS } from '../options.js';
import { exportParameters, serializeParameters, importParameters } from '../parameter-file.js';
const provenance = 'IntaRNA test\nViennaRNA test\n';

test('JSON round trip retains embedded inputs, typed options, selected columns and reports', () => {
  for (const request of [
    { ...DEFAULTS, query: '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG', qIdxPos0: -3, tIdxPos0: -5, parameters: { temperature: 22, noSeed: true, mode: 'M', qRegion: '-3-7', outSpotProbSpots: '-5&-3' }, csvColumns: ['E', 'id2'], outputs: ['normal', 'qPu', 'spotProb'] },
    { ...DEFAULTS, task: 'evaluation', rri: EVALUATION_EXAMPLE, parameters: { qShape: '1 0.2\n2 0.3' }, outputs: ['detailed', 'ensemble'] },
  ]) {
    const document = exportParameters(request, provenance);
    const imported = importParameters(JSON.stringify(document));
    const actual = prepareRequest(imported), expected = prepareRequest(request);
    actual.outputs.sort(); expected.outputs.sort();
    assert.deepEqual(actual, expected);
    assert.deepEqual(exportParameters(imported, provenance), document);
    assert.ok(!JSON.stringify(document).includes('/query.fa'));
  }
});

test('empty inputs and unfinished forms remain editable after Load/Store', () => {
  const request = { query: '', target: '', parameters: { seedBP: 12 }, csvColumns: [], outputs: [] };
  const document = serializeParameters(request, provenance);
  const loaded = importParameters(document);
  assert.equal(loaded.query, ''); assert.equal(loaded.target, '');
  assert.equal(loaded.parameters.seedBP, 12); assert.deepEqual(loaded.csvColumns, []);
  assert.throws(() => prepareRequest(loaded), /enter a sequence/);
  const defaults = exportParameters(DEFAULTS, provenance).parameters;
  for (const field of FIELDS) assert.equal(defaults[field.key], undefined, field.key);
  assert.equal(defaults.qIdxPos0, undefined);
  assert.equal(defaults.outCsvCols, undefined);
});

test('invalid/native-only settings are rejected before replacing the form', () => {
  const base = exportParameters(DEFAULTS, provenance);
  for (const values of [{ threads: 8 }, { parameterFile: 'other.json' }, { query: '/input.fa' }, { qShape: '/data.txt' }, { noSeed: 'false' }, { temperature: null }, { out: ['../file'] }, { out: ['unknown:data'] }, { personality: 'IntaRNAens' }, { energyVRNA: { text: 'x' } }, { outMode: 'D' }, { query: { text: 'a', path: 'x' } }]) {
    assert.throws(() => importParameters({ ...base, parameters: { ...base.parameters, ...values } }), JSON.stringify(values));
  }
  for (const document of [[], null, { ...base, schema: 'new' }, { ...base, other: true }, { ...base, gui: [] }, { ...base, gui: { reports: ['unknown'] } }]) assert.throws(() => importParameters(document));
  assert.throws(() => importParameters('{invalid'));
  assert.throws(() => importParameters(' '.repeat(5000001)), /5 MB/);
  assert.deepEqual(importParameters('\ufeff'+JSON.stringify(base)).csvColumns, CSV_COLUMNS.split(','));
});
