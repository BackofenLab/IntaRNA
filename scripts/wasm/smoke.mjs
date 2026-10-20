// Exercise the actual CLI and RNAlib in a fresh WASM instance.
import assert from 'node:assert/strict';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
const asset = resolve(process.argv[2] || 'build-wasm/site/wasm/intarna.js');
const { default: createIntaRNA } = await import(pathToFileURL(asset));
const diagnostics = [];
const module = await createIntaRNA({
  noInitialRun: true,
  locateFile: () => new URL('./intarna.wasm', pathToFileURL(asset)).href,
  print: (line) => diagnostics.push(line),
  printErr: (line) => diagnostics.push(line),
});
const code = module.callMain([
  '--query=GGGAGGAGGG', '--target=CCCUCCUCCC', '--outMode=C',
  '--outNumber=10', '--outOverlap=B',
  '--outCsvCols=id1,start1,end1,id2,start2,end2,E', '--out=/results.csv',
]);
assert.equal(code, 0, diagnostics.join('\n'));
const csv = module.FS.readFile('/results.csv', { encoding: 'utf8' });
assert.match(csv, /^id1;start1;end1;id2;start2;end2;E\n/);
// Reference from the native CLI with the same Turner/default-accessibility model.
assert.equal(csv.trim().split('\n')[1], 'target;1;10;query;1;10;-21.4');
assert.equal(csv.trim().split('\n').length, 11);
console.log('WASM CLI: native-reference energy/coordinates and 10 suboptimals passed');
