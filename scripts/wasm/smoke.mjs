// Exercise the actual CLI and RNAlib in fresh WASM instances.
import assert from 'node:assert/strict';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
const asset = resolve(process.argv[2] || 'build-wasm/site/wasm/intarna.js');
const { default: createIntaRNA } = await import(pathToFileURL(asset));
async function predict(args) {
  const diagnostics = [];
  const module = await createIntaRNA({
    noInitialRun: true,
    locateFile: () => new URL('./intarna.wasm', pathToFileURL(asset)).href,
    print: (line) => diagnostics.push(line),
    printErr: (line) => diagnostics.push(line),
  });
  const code = module.callMain([...args, '--outMode=C', '--outOverlap=B', '--out=/results.csv']);
  return { module, code, diagnostics: diagnostics.join('\n') };
}
const columns = '--outCsvCols=id1,start1,end1,id2,start2,end2,E';
const first = await predict([
  '--query=GGGAGGAGGG', '--target=CCCUCCUCCC', '--outNumber=10', columns,
]);
assert.equal(first.code, 0, first.diagnostics);
const csv = first.module.FS.readFile('/results.csv', { encoding: 'utf8' });
assert.match(csv, /^id1;start1;end1;id2;start2;end2;E\n/);
// Reference from the native CLI with the same Turner/default-accessibility model.
assert.equal(csv.trim().split('\n')[1], 'target;1;10;query;1;10;-21.4');
assert.equal(csv.trim().split('\n').length, 11);

// Hairpin-forming strands check nonzero accessibility costs as well as pairing.
const structured = await predict([
  '--query=GGGAAACCCGGGAAACCC', '--target=GGGUUUCCCGGGUUUCCC',
  '--outNumber=3', `${columns},ED1,ED2`,
]);
assert.equal(structured.code, 0, structured.diagnostics);
const structuredCsv = structured.module.FS.readFile('/results.csv', { encoding: 'utf8' });
assert.equal(structuredCsv.trim().split('\n')[1], 'target;1;18;query;1;18;-28.04;4.6;4.86');
assert.equal(structuredCsv.trim().split('\n').length, 4);

// Input errors must travel through the C++ exception handler, not trap WASM.
const invalid = await predict(['--query=ACGU!ACGU', '--target=ACGUACGU', columns]);
assert.notEqual(invalid.code, 0);
assert.match(invalid.diagnostics, /error/i);
// The CLI's nonzero exit status is local to the deliberately invalid instance.
process.exitCode = 0;
console.log('WASM CLI: native-reference coordinates/energy, structured accessibility, suboptimals and C++ error handling passed');
