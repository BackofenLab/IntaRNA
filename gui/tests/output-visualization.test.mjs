import { test } from 'node:test';
import assert from 'node:assert/strict';
import { outputFileContext, outputVisualization } from '../output-visualization.js';

const request = {
  target: { records: [{ id: 'target-1' }, { id: 'target-2' }] },
  query: { records: [{ id: 'query-1' }, { id: 'query-2' }] },
};

test('profile plots retain signed categorical indices and sequence names', () => {
  const spec = outputVisualization({
    name: 'qMinE-t1q2.csv',
    text: 'idx;query-2;minE\n-2;G;-3.2\n-1;C;NA\n1;A;-2.1\n',
  }, request);
  assert.equal(spec.kind, 'profile');
  assert.deepEqual(spec.traces[0].x, ['-2', '-1', '1']);
  assert.deepEqual(spec.traces[0].y, [-3.2, null, -2.1]);
  assert.equal(spec.layout.xaxis.type, 'category');
  assert.deepEqual(spec.layout.xaxis.categoryarray, ['-2', '-1', '1']);
  assert.match(spec.layout.title.text, /query-2/);
  assert.equal(spec.context.label, 'target-1 : query-2');
});

test('pair energy and probability matrices use literal category axes', () => {
  const energy = outputVisualization({
    name: 'pMinE-t2q1.csv',
    text: 'minE;A_-1;C_1\nG_5;-4.2;NA\nU_6;-2.0;-1.5\n',
  }, request);
  assert.equal(energy.kind, 'matrix');
  assert.deepEqual(energy.traces[0].z, [[-4.2, null], [-2, -1.5]]);
  assert.equal(energy.layout.xaxis.type, 'category');
  assert.equal(energy.layout.yaxis.type, 'category');
  assert.match(energy.layout.title.text, /target-2 : query-1/);

  const spots = outputVisualization({
    name: 'spotProb-t1q2.csv',
    text: 'spotProb;A_-1;C_1\nG_5;0.1;0.2\n',
  }, request);
  assert.equal(spots.kind, 'matrix');
  assert.equal(spots.traces[0].type, 'heatmap');
  assert.deepEqual(spots.layout.xaxis.categoryarray, ['A_-1', 'C_1']);
});

test('selected spot probabilities use a categorical bar plot', () => {
  const onePair = { target: { records: [{ id: 'target' }] }, query: { records: [{ id: 'query' }] } };
  const spec = outputVisualization({ name: 'spotProb.csv', text: 'spot;probability\n0&0;0.4\n5&-2;0.6\n' }, onePair);
  assert.equal(spec.kind, 'spots');
  assert.equal(spec.traces[0].type, 'bar');
  assert.deepEqual(spec.traces[0].x, ['0&0', '5&-2']);
  assert.equal(spec.context.label, 'target : query');
});

test('output headers identify the sequence or sequence pair they cover', () => {
  assert.equal(outputFileContext('tMinE-t2q1.csv', request).label, 'target-2 : query-1');
  assert.equal(outputFileContext('qAcc-s2.csv', request).label, 'query-2');
  assert.equal(outputFileContext('intarna-detailed.txt', request).label, '');
});

test('non-plot files are ignored and malformed plot data fails clearly', () => {
  assert.equal(outputVisualization({ name: 'qAcc.csv', text: '# accessibility' }, request), null);
  assert.throws(() => outputVisualization({ name: 'qMinE-t1q1.csv', text: 'idx;query-1;wrong\n1;A;2\n' }, request), /Unexpected qMinE/);
  assert.throws(() => outputVisualization({ name: 'pMinE-t1q1.csv', text: 'minE;A_1\nG_1;not-a-number\n' }, request), /non-numeric/);
});
