import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { DEFAULTS, EVALUATION_EXAMPLE, parseSequence, parseFasta, prepareRequest, makeArguments, CSV_COLUMNS, sequenceOffset } from '../protocol.js';
import { FIELDS, GROUPS } from '../options.js';
import { CSV_CATALOG, CSV_GROUPS, defaultColumns } from '../columns.js';
const request = (parameters = {}, extra = {}) => ({ ...DEFAULTS, parameters, ...extra });

test('FASTA normalization, identifiers, bounded multi-record work and subsets', () => {
  assert.equal(parseSequence('>test\r\nacgt n\r\nACG', 'Query').fasta, '>test\nACGUNACG\n');
  for (const input of ['', 'A', '>\nACGUACGU', '>one\n', 'ACGU-ACGU', '>a\nACGU\n>a\nACGU', '>a;b\nACGU\n>a_b\nACGU', '>a first\nACGU\n>a second\nACGU']) assert.throws(() => parseFasta(input, 'Query'));
  const fasta = '>a\nGGGAGGAGGG\n>b\nGGGAGGAGGG';
  assert.equal(prepareRequest(request({}, { query: fasta })).query.records.length, 2);
  assert.throws(() => prepareRequest(request({}, { task: 'evaluation', query: fasta, rri: EVALUATION_EXAMPLE })), /exactly one/);
  assert.throws(() => prepareRequest(request({}, { query: '>a\n'+'A'.repeat(5000)+'\n>b\n'+'A'.repeat(2000), target: 'C'.repeat(300) })), /1,000,000/);
  assert.doesNotThrow(() => prepareRequest(request({ qSet: '2-2' }, { query: '>a\n'+'A'.repeat(5000)+'\n>b\n'+'A'.repeat(2000), target: 'C'.repeat(300) })));
  assert.throws(() => prepareRequest(request({ outPairwise: true }, { query: fasta })), /same number/);
  assert.throws(() => prepareRequest(request({ qSet: '3-3' })), /within/);
  assert.throws(() => prepareRequest(request({}, { query: 'A'.repeat(5001) })));
});

test('CLI defaults, ranges and strict typed controls reject invalid values', () => {
  const prepared = prepareRequest(DEFAULTS);
  assert.equal(prepared.parameters.outNumber, 1);
  assert.equal(prepared.parameters.outNumberPerSeed, 0);
  assert.equal(prepared.parameters.qIntLoopMax, 16);
  assert.equal(prepared.parameters.tIntLoopMax, 10);
  assert.equal(prepared.csvColumns.join(','), CSV_COLUMNS);
  assert.ok(!makeArguments(prepared).some((arg) => arg.startsWith('--personality=')));
  assert.ok(!makeArguments(prepared).some((arg) => arg.startsWith('--kineticScore=')));
  for (const parameters of [{ outNumber: 101 }, { outNumber: 1.5 }, { outNumberPerSeed: -1 }, { outNumberPerSeed: 1001 }, { outNumberPerSeed: 1.5 }, { temperature: '' }, { temperature: Infinity }, { qAccW: 2 }, { qAccW: 50 }, { model: '--help' }, { noSeed: 'false' }, { seedBP: 11 }]) assert.throws(() => prepareRequest(request(parameters)));
  for (const qIdxPos0 of ['--help', NaN, 2e9+1, '', false]) assert.throws(() => prepareRequest(request({}, { qIdxPos0 })));
  assert.doesNotThrow(() => prepareRequest(request({}, { qIdxPos0: -2e9 })));
  assert.doesNotThrow(() => prepareRequest(request({}, { qIdxPos0: 0 })));
  assert.equal(sequenceOffset(1, -3), 3);
  assert.equal(sequenceOffset(-1, -3), 2);
  assert.throws(() => prepareRequest(request({}, { csvColumns: [] })), /at least one/);
  assert.throws(() => prepareRequest(request({}, { outputs: ['../other'] })), /available/);
});

test('evaluation validates antiparallel supplied structures in signed coordinates', () => {
  const evalRequest = request({ seedBP: 'ignored', model: 'ignored' }, { task: 'evaluation', rri: EVALUATION_EXAMPLE });
  const prepared = prepareRequest(evalRequest);
  assert.equal(prepared.parameters.model, undefined);
  assert.ok(makeArguments(prepared).includes('--personality=IntaRNAeval'));
  assert.ok(makeArguments(prepared).includes(`--rri=${EVALUATION_EXAMPLE}`));
  for (const rri of ['', '1|||&1||', '1|||||||||||&1|||||||||||', '0|||&1|||', '1|||&1|||:', '1|||&1...']) assert.throws(() => prepareRequest({ ...evalRequest, rri }));
  assert.doesNotThrow(() => prepareRequest({ ...evalRequest, qIdxPos0: -3, tIdxPos0: -5, rri: '-5||||||||||&-3||||||||||' }));
  assert.throws(() => prepareRequest({ ...evalRequest, outputs: ['bpsvg'] }), /available/);
});

test('cross-option constraints, regions, probabilities, SHAPE and local file arguments', () => {
  for (const parameters of [{ noSeed: true, mode: 'S' }, { noSeed: true, outNumberPerSeed: 1 }, { mode: 'K', model: 'P' }, { model: 'B', mode: 'M' }, { model: 'B', helixMaxBP: 4 }, { qRegion: '8-2' }, { qAccConstr: 'b:1-20' }, { qAccConstr: 'bbb' }, { qShapeMethod: 'X' }, { tShapeConversion: 'bad' }, { qShape: '1 nope' }, { qAcc: 'P' }, { qAcc: 'P', qAccFile: 'bad' }, { outSpotProbSpots: '20&1' }, { windowWidth: 100 }, { qIntLenMax: 3 }, { qAcc: 'N', qRegionLenMax: 20 }, { qRegion: '1-10', qRegionLenMax: 20 }, { seedTQ: '1.|||||||||&1|||||||||.' }]) assert.throws(() => prepareRequest(request(parameters)));
  assert.doesNotThrow(() => prepareRequest(request({ noSeed: true, qAccConstr: 'b:1-2,p:8-9', qRegion: '1-10' })));
  assert.doesNotThrow(() => prepareRequest(request({ outNumberPerSeed: 1, outNumber: 3, outOverlap: 'N', windowWidth: 200 })));
  assert.doesNotThrow(() => prepareRequest(request({ outNumberPerSeed: 1, outOverlap: 'N', qRegion: '1-4,6-10' })));
  assert.throws(() => prepareRequest(request({ mode: 'K' }, { csvColumns: ['Eall'] })), /ensemble/);
  assert.throws(() => prepareRequest(request({}, { outputs: ['bpsvg'] })), /model P/);
  assert.throws(() => prepareRequest(request({ outCsvSort: 'ED1' })), /sort column/);
  const prepared = prepareRequest(request({ model: 'P', mode: 'M', qShape: '1 0.2\n2 -999', outSpotProbSpots: '1&2' }, { outputs: ['bpProb', 'bpsvg', 'spotProb'] }));
  const args = makeArguments(prepared);
  assert.ok(args.includes('--qShape=/qShape.txt'));
  assert.ok(args.includes('--out=spotProb:1&2:/output/spotProb.csv'));
  assert.ok(args.includes('--out=bpsvg:/output/bpsvg.svg'));
  assert.ok(!makeArguments(prepared, 'D').some((arg) => arg.startsWith('--outCsvCols')));
});

test('output controls preserve the reviewed three-row filter order', () => {
  const rows = GROUPS.find((group) => group.key === 'output-options').rows.map((row) => row.map((field) => field.key));
  assert.deepEqual(rows.slice(0, 3), [
    ['outNumber', 'outNumberPerSeed', 'outOverlap'],
    ['outMaxE', 'outDeltaE', 'outMinPu'],
    ['outBestSeedOnly', 'outNoLP', 'outNoGUend', 'outPerRegion', 'outPairwise'],
  ]);
});

test('browser catalog covers every native CSV column and every control names a CLI option', () => {
  const header = readFileSync(new URL('../../src/IntaRNA/OutputHandlerCsv.h', import.meta.url), 'utf8');
  const ids = [...header.matchAll(/colType2string\[\w+\] = "(\w+)"/g)].map((match) => match[1]);
  assert.deepEqual(CSV_CATALOG.map((c) => c.key).toSorted(), ids.toSorted());
  const source = readFileSync(new URL('../../src/bin/CommandLineParsing.cpp', import.meta.url), 'utf8');
  for (const field of FIELDS.filter((field) => field.key !== 'outSpotProbSpots')) assert.ok(source.includes(`"${field.key}"`), field.key);
});

test('JSON worker boundary embeds all inputs with typed non-default settings', async () => {
  const { makeParameterDocument, JSON_SCHEMA } = await import('../protocol.js');
  const prepared = prepareRequest(request({ temperature: 20, noSeed: true, qShape: '1 0.2' }, { outputs: ['qPu'] }));
  const document = JSON.parse(JSON.stringify(makeParameterDocument(prepared)));
  assert.equal(document.schema, JSON_SCHEMA);
  assert.equal(document.parameters.query.text, prepared.query.fasta);
  assert.equal(document.parameters.target.text, prepared.target.fasta);
  assert.equal(document.parameters.qShape.text, '1 0.2\n');
  assert.equal(document.parameters.temperature, 20);
  assert.equal(document.parameters.noSeed, true);
  assert.equal(document.parameters.qIdxPos0, undefined);
  assert.deepEqual(document.parameters.out, ['/results.csv', 'qPu:/output/qPu.csv']);
  const worker = readFileSync(new URL('../intarna-worker.js', import.meta.url), 'utf8');
  assert.ok(worker.includes("callMain(['--parameterFile=/parameter.json'])"));
});


test('browser batches use longest selected strands and round-trip full FASTA inputs', async () => {
  const batch = Array.from({ length: 50 }, (_, i) => `>q${i}\n${'A'.repeat(5000)}\n`).join('');
  assert.equal(parseFasta(batch, 'Query').records.length, 50);
  assert.throws(() => parseFasta(batch + '>extra\nAA', 'Query'), /50 FASTA/);
  assert.doesNotThrow(() => prepareRequest(request({}, { query: batch, target: 'C'.repeat(200) })));
  assert.throws(() => prepareRequest(request({}, { query: batch, target: 'C'.repeat(201) })), /1,000,000/);
  const crossed = { query: '>q1\n' + 'A'.repeat(5000) + '\n>q2\nAA', target: '>t1\nCC\n>t2\n' + 'C'.repeat(5000) };
  assert.throws(() => prepareRequest(request({ outPairwise: true, noSeed: true }, crossed)), /1,000,000/);
  const { serializeParameters, importParameters } = await import('../parameter-file.js');
  assert.equal(importParameters(serializeParameters(request({}, { query: batch }), 'test')).query, batch);
  const parsed = parseFasta('>one\tdescription here\nACGU\n>two   \nACGU', 'Query');
  assert.deepEqual(parsed.records.map((record) => record.id), ['one', 'two']);
  assert.equal(parsed.warnings.length, 1);
  assert.match(parsed.warnings[0], /one/);
  assert.equal(parsed.fasta, '>one\nACGU\n>two\nACGU\n');
});

test('CSV groups cover native columns once and task defaults serialize explicitly', async () => {
  assert.deepEqual(CSV_GROUPS.flatMap((group) => group.keys).toSorted(), CSV_CATALOG.map((column) => column.key).toSorted());
  const energy = [...CSV_GROUPS.find((group) => group.label === 'Energy').keys, 'Etotal', 'Eall', 'Eall1', 'Eall2', 'EallTotal'];
  for (const key of energy) assert.ok(defaultColumns('evaluation').includes(key), key);
  for (const key of ['seedStart1', 'seedStart2', 'seedEnd1', 'seedEnd2']) assert.ok(!defaultColumns('evaluation').includes(key), key);
  for (const key of ['ED1', 'ED2', 'seedE', 'seedStart1', 'seedEnd1', 'seedStart2', 'seedEnd2']) assert.ok(defaultColumns('prediction').includes(key), key);
  const { exportParameters, importParameters } = await import('../parameter-file.js');
  const input = request({ outNumber: 10, outNumberPerSeed: 2 }, { csvColumns: defaultColumns('prediction') });
  const saved = exportParameters(input, 'test');
  assert.equal(saved.parameters.outNumber, 10);
  assert.equal(saved.parameters.outNumberPerSeed, 2);
  assert.equal(saved.parameters.outCsvCols, input.csvColumns.join(','));
  assert.deepEqual(prepareRequest(importParameters(saved)), prepareRequest(input));
});
