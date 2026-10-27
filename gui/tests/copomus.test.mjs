import { test } from 'node:test';
import assert from 'node:assert/strict';
import { mutationAnnotations, mutationEnergyPlot } from '../copomus.js';

const candidate = {
  values: { mutation: 'G-3C&C5G', qIndex: -3, tIndex: 5, bpMutated: 'CG' },
  interactions: {
    ww: { E: -12.5 }, wm: { E: -6.25 }, mw: null, mm: { E: -11 },
  },
};

test('CopomuS annotations preserve signed query-first coordinates in vaRRI target-first order', () => {
  assert.equal(mutationAnnotations(candidate), '1:5G:006400,2:-3C:006400');
  for (const values of [
    { ...candidate.values, qIndex: 0 },
    { ...candidate.values, tIndex: 1.5 },
    { ...candidate.values, bpMutated: 'NN' },
  ]) assert.throws(() => mutationAnnotations({ ...candidate, values }));
});

test('CopomuS energy plot keeps requested state order and wild-type baseline', () => {
  const plot = mutationEnergyPlot(candidate);
  assert.deepEqual(plot.data[0].x, ['ww', 'wm', 'mw', 'mm']);
  assert.deepEqual(plot.data[0].y, [-12.5, -6.25, null, -11]);
  assert.equal(plot.layout.xaxis.title, undefined);
  assert.equal(plot.layout.yaxis.title.text, 'E [kcal/mol]');
  assert.deepEqual(plot.layout.shapes[0], {
    type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: -12.5, y1: -12.5,
    line: { color: '#145a96', width: 2, dash: 'dot' },
  });
  assert.deepEqual(mutationEnergyPlot({ ...candidate, interactions: { ...candidate.interactions, ww: null } }).layout.shapes, []);
});
