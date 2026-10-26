// Integration coverage of browser argument translation against the real CLI.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { DEFAULTS, CSV_COLUMNS, EVALUATION_EXAMPLE, prepareRequest, makeParameterDocument } from '../../gui/protocol.js';
import { CSV_CATALOG } from '../../gui/columns.js';
import { parseResults, projectCsv } from '../../gui/results.js';
import { archiveFiles } from '../../gui/archive.js';
const asset = resolve(process.argv[2] || 'build-wasm/site/wasm/intarna.js');
const { default: createIntaRNA } = await import(pathToFileURL(asset));
const compiledVersion = readFileSync(new URL('../build-info.txt', pathToFileURL(asset)), 'utf8');
async function run(input, format = 'C') {
  const request = prepareRequest({ ...DEFAULTS, ...input });
  const diagnostics = [];
  const module = await createIntaRNA({ noInitialRun: true, locateFile: () => new URL('./intarna.wasm', pathToFileURL(asset)).href, print: (line) => diagnostics.push(line), printErr: (line) => diagnostics.push(line) });
  module.FS.writeFile('/parameter.json', JSON.stringify(makeParameterDocument(request, format)));
  module.FS.mkdir('/output');
  const args = ['--parameterFile=/parameter.json'];
  assert.equal(module.callMain(args), 0, `${args.join(' ')}\n${diagnostics.join('\n')}`);
  const text = module.FS.readFile(format === 'C' ? '/results.csv' : '/report.txt', { encoding: 'utf8' });
  return { rows: format === 'C' ? parseResults(text) : [], text, fs: module.FS, request };
}
const defaults = await run({ parameters: { outNumber: 10 } });
assert.equal(defaults.rows.length, 10);
assert.equal(defaults.rows[0].E, -21.4);
const shifted = await run({ qIdxPos0: -3, tIdxPos0: -5 });
assert.deepEqual(['start1', 'end1', 'start2', 'end2'].map((key) => shifted.rows[0][key]), [-5, 5, -3, 7]);
const multi = await run({ query: '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG', target: '>t1\nCCCUCCUCCC\n>t2\nCCCUCCUCCC', outputs: ['qMinE', 'tMinE'] });
assert.equal(multi.rows.length, 4);
assert.ok(multi.rows.every((row) => row.deltaE === 0));
assert.equal(multi.fs.readdir('/output').filter((name) => name.endsWith('.csv')).length, 8);
const pairwise = await run({ query: '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG', target: '>t1\nCCCUCCUCCC\n>t2\nCCCUCCUCCC', parameters: { outPairwise: true } });
assert.equal(pairwise.rows.length, 2);
const evaluation = await run({ task: 'evaluation', rri: EVALUATION_EXAMPLE + ':1|||||||||&2|||||||||', outputs: ['qAcc', 'tPu'], csvColumns: CSV_CATALOG.map((c) => c.key) });
assert.equal(evaluation.rows.length, 2);
assert.equal(evaluation.rows[0].E, -21.4);
assert.equal(evaluation.rows[0].seedE, 'NAN');
assert.ok((await run({ task: 'evaluation', rri: EVALUATION_EXAMPLE }, 'D')).text.includes('-21.4'));
assert.ok((await run({ task: 'evaluation', rri: EVALUATION_EXAMPLE }, 'E')).text.length > 0);
const probabilities = await run({ parameters: { model: 'P' }, csvColumns: CSV_CATALOG.map((c) => c.key), outputs: ['bpProb', 'bpsvg', 'qSpotProb', 'tSpotProb', 'spotProb', 'pMinE'] });
assert.equal(probabilities.rows.length, 1);
assert.ok(probabilities.fs.readFile('/output/bpsvg.svg', { encoding: 'utf8' }).includes('<svg'));
const accessibility = await run({ parameters: { outNumber: 0 }, outputs: ['qPu'] });
const imported = await run({ parameters: { qAcc: 'P', qAccFile: accessibility.fs.readFile('/output/qPu.csv', { encoding: 'utf8' }) } });
assert.equal(imported.rows[0].E, (await run({ parameters: { qIntLenMax: 9 } })).rows[0].E);
for (const parameters of [{ energyVRNA: 'Turner99', temperature: 20 }, { energyVRNA: 'Andronescu07' }, { model: 'B' }, { mode: 'K' }, { noSeed: true, mode: 'M' }, { seedMaxUP: 1 }, { qAccConstr: '..........', qShape: '1 0.2\n2 0.3' }, { outNumberPerSeed: 1 }, { outOverlap: 'N', outNumber: 3 }, { outCsvSort: 'E', outNumber: 3 }]) assert.ok((await run({ parameters })).rows.length >= 1);
console.log('Browser protocol: native-reference prediction/evaluation, signed coordinates, multi-FASTA, all columns, imported accessibility, energy/seed/helix/SHAPE options and output files passed');

// Execute exported parameter files in fresh engines with only archived inputs.
// Byte comparisons cover selected CSV, per-record extra files and report modes.
async function replayArchive(input) {
  const original = await run(input), { request, fs } = original;
  const outputs = fs.readdir('/output').filter((name) => name !== '.' && name !== '..').map((name) => ({ name, text: fs.readFile(`/output/${name}`, { encoding: 'utf8' }) }));
  const reports = [['normal', 'N'], ['detailed', 'D'], ['ensemble', 'E']].filter(([key]) => request.outputs.includes(key));
  for (const [key, format] of reports) outputs.push({ name: `intarna-${key}.txt`, text: (await run(input, format)).text });
  const files = archiveFiles(request, { csv: original.text, files: outputs }, compiledVersion);
  for (const config of ['parameters', ...reports.map(([key]) => key)]) {
    const diagnostics = [];
    const module = await createIntaRNA({ noInitialRun: true, locateFile: () => new URL('./intarna.wasm', pathToFileURL(asset)).href, print: (line) => diagnostics.push(line), printErr: (line) => diagnostics.push(line) });
    module.FS.writeFile('/parameter.json', files.find((file) => file.name === 'parameter.json').text);
    const arguments_ = ['--parameterFile=/parameter.json', `--out=${config === 'parameters' ? 'intarna-interactions.csv' : `intarna-${config}.txt`}`];
    if (config !== 'parameters') arguments_.push(`--outMode=${reports.find(([key]) => key === config)[1]}`, `--outCsvCols=${CSV_COLUMNS}`, "--outCsvSort", "");
    assert.equal(module.callMain(arguments_), 0, diagnostics.join('\n'));
    const expected = config === 'parameters' ? [{ name: 'intarna-interactions.csv', text: projectCsv(original.text, request.csvColumns) }, ...outputs.filter((file) => !file.name.startsWith('intarna-'))] : outputs.filter((file) => file.name === `intarna-${config}.txt`);
    for (const file of expected) assert.equal(module.FS.readFile(`/${file.name}`, { encoding: 'utf8' }), file.text, file.name);
  }
}
await replayArchive({ query: '>q1\nGGGAGGAGGG\n>q2\nGGGAGGAGGG', qIdxPos0: -3, tIdxPos0: -5, parameters: { temperature: 20, outNumber: 3, noSeed: true }, outputs: ['qMinE', 'tPu', 'normal', 'detailed'], csvColumns: ['id2', 'E'] });
await replayArchive({ task: 'evaluation', rri: EVALUATION_EXAMPLE, parameters: { qShape: '1 0.2\n2 0.3' }, outputs: ['detailed', 'ensemble'] });
await replayArchive({ parameters: { qAcc: 'P', qAccFile: accessibility.fs.readFile('/output/qPu.csv', { encoding: 'utf8' }), model: 'P' }, outputs: ['bpProb', 'bpsvg'], csvColumns: ['E'] });
console.log('ZIP parameter replay: CSV, multi-FASTA/signed indices, non-default energy/seed options, SHAPE/accessibility imports, per-record outputs, probability SVG and ASCII/detailed/ensemble reports are byte-identical');
