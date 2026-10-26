import { test } from 'node:test';
import assert from 'node:assert/strict';
import { DEFAULTS, prepareRequest } from '../protocol.js';
import { archiveFiles } from '../archive.js';
import { visualizationUrl } from '../visualization.js';

test('visualization retains full strand orientation, bulges, names and indices crossing zero', () => {
  const request = prepareRequest({ ...DEFAULTS, query: '>other\nGGGAGGAGGG\n>q;β\nGGGAGGAGGG', qIdxPos0: -3, tIdxPos0: -5 });
  const url = new URL(visualizationUrl({ id1: 'target', id2: 'q_β', start1: -3, end1: 3, start2: -2, end2: 4, E: -1, hybridDP: '((..((&)).)).' }, request));
  assert.equal(url.origin + url.pathname, 'https://backofenlab.github.io/vaRRI/');
  assert.equal(url.searchParams.get('sequence'), 'CCCUCCUCCC&GGGAGGAGGG');
  assert.equal(url.searchParams.get('structure'), '..((..((..&.)).))....');
  assert.equal(url.searchParams.get('seqName2'), 'q_β');
  assert.equal(url.searchParams.get('startIndex1'), '-5');
  assert.equal(url.searchParams.get('startIndex2'), '-3');
  assert.equal(url.searchParams.get('showRenderingOnly'), 'true');
  assert.equal(url.searchParams.has('regionHighlights'), false);
  assert.throws(() => visualizationUrl({ id1: 'missing', hybridDP: '(&)' }, request), /not found/);
  assert.throws(() => visualizationUrl({ id1: 'target', id2: 'q_β', start1: -5, end1: -4, hybridDP: '(((&)))' }, request), /coordinates/);
  assert.throws(() => visualizationUrl({ id1: 'target', hybridDP: '(&)' }, { ...request, parameters: { ...request.parameters, tIdxPos0: 0 } }), /nonzero/);
});

test('visualization highlights every available seed region in transparent purple', () => {
  const request = prepareRequest({ ...DEFAULTS, query: 'G'.repeat(20), target: 'C'.repeat(20), qIdxPos0: -3, tIdxPos0: -5 });
  const row = {
    id1: 'target', id2: 'query', start1: -5, end1: 15, start2: -3, end2: 17,
    E: -1, hybridDP: `${'('.repeat(10)}${'.'.repeat(10)}&${'.'.repeat(10)}${')'.repeat(10)}`,
    seedStart1: '-3:5', seedEnd1: '2:8', seedStart2: '-1:7', seedEnd2: '4:10',
  };
  const url = new URL(visualizationUrl(row, request));
  assert.equal(url.searchParams.get('regionHighlights'), '-3-2&-1-4:800080:0.3,5-8&7-10:800080:0.3');
  row.seedEnd2 = 'NAN';
  assert.equal(new URL(visualizationUrl(row, request)).searchParams.has('regionHighlights'), false);
});

test('archive preserves effective non-default settings, selected CSV, auxiliary inputs and reports', () => {
  const request = prepareRequest({ ...DEFAULTS, task: 'evaluation', rri: '1||||||||||&3||||||||||', qIdxPos0: 3, parameters: { temperature: 21, qShape: '1 0.2\n2 0.3' }, csvColumns: ['E'] });
  request.outputs = ['detailed'];
  const files = archiveFiles(request, { csv: 'id1;start1;end1;id2;start2;end2;E;hybridDP\ntarget;1;10;query;3;12;-21.4;((((((((((&))))))))))\n', files: [{ name: 'intarna-detailed.txt', text: 'report' }] }, 'IntaRNA 3.4.1\nViennaRNA 2.7.0\n');
  const entries = Object.fromEntries(files.map((file) => [file.name, file.text]));
  assert.deepEqual(Object.keys(entries), ['parameter.json', 'intarna-interactions.csv', 'intarna-detailed.txt']);
  const document = JSON.parse(entries['parameter.json']);
  assert.equal(document.parameters.query.text, '>query\nGGGAGGAGGG\n');
  assert.equal(entries['intarna-interactions.csv'], 'E\n-21.4\n');
  assert.equal(document.parameters.qShape.text, '1 0.2\n2 0.3');
  assert.ok(document._comment.includes('ViennaRNA 2.7.0'));
  for (const [key, value] of Object.entries({ personality: 'IntaRNAeval', temperature: 21, qIdxPos0: 3, outCsvCols: 'E' })) assert.equal(document.parameters[key], value);
  assert.equal(document.parameters.tIdxPos0, undefined);
  assert.equal(document.parameters.seedBP, undefined);
  assert.deepEqual(document.gui.reports, ['detailed']);
  assert.ok(document._comment.some((line) => line.includes('IntaRNA --parameterFile=parameter.json')));
  assert.equal(entries['intarna-detailed.txt'], 'report');
  assert.throws(() => archiveFiles(request, { csv: '', files: [] }, ''), /version/);
});
