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
  assert.throws(() => visualizationUrl({ id1: 'missing', hybridDP: '(&)' }, request), /not found/);
  assert.throws(() => visualizationUrl({ id1: 'target', id2: 'q_β', start1: -5, end1: -4, hybridDP: '(((&)))' }, request), /coordinates/);
  assert.throws(() => visualizationUrl({ id1: 'target', hybridDP: '(&)' }, { ...request, parameters: { ...request.parameters, tIdxPos0: 0 } }), /nonzero/);
});

test('archive preserves effective non-default settings, selected CSV, auxiliary inputs and reports', () => {
  const request = prepareRequest({ ...DEFAULTS, task: 'evaluation', rri: '1||||||||||&3||||||||||', qIdxPos0: 3, parameters: { temperature: 21, qShape: '1 0.2\n2 0.3' }, csvColumns: ['E'] });
  request.outputs = ['detailed'];
  const files = archiveFiles(request, { csv: 'id1;start1;end1;id2;start2;end2;E;hybridDP\ntarget;1;10;query;3;12;-21.4;((((((((((&))))))))))\n', files: [{ name: 'intarna-detailed.txt', text: 'report' }] }, 'IntaRNA 3.4.1\nViennaRNA 2.7.0\n');
  const entries = Object.fromEntries(files.map((file) => [file.name, file.text]));
  assert.equal(entries['query.fa'], '>query\nGGGAGGAGGG\n');
  assert.equal(entries['intarna-interactions.csv'], 'E\n-21.4\n');
  assert.equal(entries['qShape.txt'], '1 0.2\n2 0.3\n');
  assert.match(entries['parameters.cfg'], /# ViennaRNA 2.7.0/);
  for (const line of ['personality = IntaRNAeval', 'temperature = 21', 'qIdxPos0 = 3', 'qShape = qShape.txt', 'outCsvCols = E']) assert.ok(entries['parameters.cfg'].includes(line), line);
  assert.ok(!entries['parameters.cfg'].includes('tIdxPos0 ='));
  assert.ok(!entries['parameters.cfg'].includes('seedBP ='));
  assert.match(entries['detailed.cfg'], /out = intarna-detailed.txt/);
  assert.match(entries['parameters.cfg'], /IntaRNA --parameterFile=parameters.cfg/);
  assert.equal(entries['intarna-detailed.txt'], 'report');
  assert.throws(() => archiveFiles(request, { csv: '', files: [] }, ''), /version/);
});
