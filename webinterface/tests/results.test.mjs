import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseResults, selectRows } from '../results.js';
const header = 'id1;start1;end1;id2;start2;end2;E\n';

test('native CSV orientation, numeric coordinates, energy and relative suboptimality', () => {
  const rows = parseResults(header + 'target;2;10;query;1;9;-19.6\ntarget;1;10;query;1;10;-21.4\n');
  assert.equal(rows[0].start1, 2);
  assert.equal(rows[0].start2, 1);
  assert.equal(rows[0].deltaE, 1.8);
  assert.equal(rows[1].deltaE, 0);
  assert.deepEqual(selectRows(rows, '', 'E', 1).map((row) => row.E), [-21.4, -19.6]);
  assert.deepEqual(selectRows(rows, '', 'start1', -1).map((row) => row.start1), [2, 1]);
  assert.equal(selectRows(rows, 'QUERY', 'id2', 1).length, 2);
  assert.equal(selectRows(rows, 'missing', 'id2', 1).length, 0);
});

test('empty output and literal native ID quotes; malformed output fails clearly', () => {
  assert.deepEqual(parseResults(''), []);
  assert.deepEqual(parseResults(header), []);
  assert.equal(parseResults(header + '"name";1;7;query;1;7;-3\n')[0].id1, '"name"');
  for (const csv of ['unexpected\n', header + 'too;few\n', header + 't;0;1;q;1;2;-3\n', header + 't;2;1;q;1;2;-3\n', header + 't;1;7;q;1;7;NAN\n', header + 't;1;7;q;1;7;\n']) {
    assert.throws(() => parseResults(csv));
  }
});
