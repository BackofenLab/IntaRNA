import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFastaFile } from '../uploads.js';
const file = (text) => ({ size: new TextEncoder().encode(text).length, text: async () => text });
test('FASTA uploads accept text/CRLF and validate size, format and task before replacement', async () => {
  const input = '>q\r\ngggaggaggg\r\n';
  assert.equal(await readFastaFile(file(input), 'Query', 'prediction'), input);
  const multi = '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG';
  assert.equal(await readFastaFile(file(multi), 'Query', 'prediction'), multi);
  await assert.rejects(readFastaFile(file(multi), 'Query', 'evaluation'), /exactly one/);
  for (const invalid of ['', '\u0000bin', '>q\n--!', '>q\n']) await assert.rejects(readFastaFile(file(invalid), 'Query', 'prediction'));
  const batch = Array.from({ length: 50 }, (_, i) => `>q${i}\n${'A'.repeat(5000)}\n`).join('');
  assert.equal(await readFastaFile(file(batch), 'Query', 'prediction'), batch);
  await assert.rejects(readFastaFile({ size: 1000001, text: () => { throw Error('must not read'); } }, 'Query', 'prediction'), /1 MB/);
});
