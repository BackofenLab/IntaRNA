import { prepareRequest, makeParameterDocument } from './protocol.js';

const PYTHON_FILES = [
  'intarna_tools.py',
  'copomus/__init__.py', 'copomus/browser.py', 'copomus/workflow.py',
  'copomus/IntaRNA.py', 'copomus/indexing.py', 'copomus/candidate_filters.py',
  'copomus/measures.py', 'copomus/mutation.py', 'copomus/mutation_generators.py',
];
const CALL_COLUMNS = new Set(['id1', 'start1', 'end1', 'id2', 'start2', 'end2', 'subseqDP', 'hybridDP', 'E', 'bpList']);
const CALL_OUTPUTS = new Set(['qAcc', 'qMinE', 'tMinE']);
let started = false;

self.onmessage = async ({ data }) => {
  if (started) return;
  started = true;
  const diagnostics = [];
  const log = (line) => {
    diagnostics.push(String(line).slice(0, 2000));
    if (diagnostics.length > 50) diagnostics.shift();
  };
  try {
    const request = prepareRequest(data);
    if (request.task !== 'mutation') throw new Error('The mutation worker only accepts Mutation Design requests.');
    self.postMessage({ type: 'status', message: 'Loading Python and the prediction engine…' });
    const [{ loadPyodide }, { default: createIntaRNA }] = await Promise.all([
      import('./vendor/pyodide/pyodide.mjs'), import('./wasm/intarna.js'),
    ]);
    const pyodide = await loadPyodide({
      indexURL: new URL('./vendor/pyodide/', import.meta.url).href,
      stdout: log, stderr: log,
    });
    pyodide.FS.mkdirTree('/app/copomus');
    await Promise.all(PYTHON_FILES.map(async (name) => {
      const response = await fetch(new URL(`./python/${name}`, import.meta.url));
      if (!response.ok) throw new Error(`Could not load CopomuS Python module ${name}.`);
      pyodide.FS.writeFile(`/app/${name}`, await response.text(), { encoding: 'utf8' });
    }));
    pyodide.runPython("import sys; sys.path.insert(0, '/app')");

    const predict = async (payloadJson) => {
      const call = JSON.parse(String(payloadJson));
      if (typeof call.query !== 'string' || typeof call.target !== 'string'
          || call.qIdxPos0 !== request.parameters.qIdxPos0 || call.tIdxPos0 !== request.parameters.tIdxPos0
          || !Array.isArray(call.columns) || !call.columns.length || call.columns.some((column) => !CALL_COLUMNS.has(column))
          || !Number.isInteger(call.outNumber) || call.outNumber < 1 || call.outNumber > 1000
          || !Array.isArray(call.outputs) || call.outputs.some((output) => !CALL_OUTPUTS.has(output))
          || call.parameters === null || typeof call.parameters !== 'object' || Array.isArray(call.parameters)
          || Object.keys(call.parameters).some((key) => !['qRegion', 'tRegion'].includes(key))) {
        return JSON.stringify({ error: 'CopomuS sent an invalid IntaRNA request.' });
      }
      const module = await createIntaRNA({
        noInitialRun: true,
        locateFile: () => new URL('./wasm/intarna.wasm', import.meta.url).href,
        print: log, printErr: log,
      });
      const document = makeParameterDocument(request);
      const parameters = document.parameters;
      parameters.query = { text: `>${request.query.records[0].id}\n${call.query}\n` };
      parameters.target = { text: `>${request.target.records[0].id}\n${call.target}\n` };
      parameters.outMode = 'C'; parameters.outCsvCols = call.columns.join(','); parameters.outNumber = call.outNumber;
      delete parameters.outCsvSort;
      Object.assign(parameters, call.parameters);
      parameters.out = ['/results.csv', ...call.outputs.map((output) => `${output}:/output/${output}.txt`)];
      module.FS.mkdir('/output');
      module.FS.writeFile('/parameter.json', JSON.stringify(document));
      const code = module.callMain(['--parameterFile=/parameter.json']);
      if (code !== 0) return JSON.stringify({ error: `IntaRNA exited with status ${code}. Check the diagnostics and settings.` });
      const files = Object.fromEntries(call.outputs.map((output) => [output, module.FS.readFile(`/output/${output}.txt`, { encoding: 'utf8' })]));
      return JSON.stringify({ csv: module.FS.readFile('/results.csv', { encoding: 'utf8' }), files });
    };
    const progress = (message) => self.postMessage({ type: 'status', message: String(message) });
    const payload = JSON.stringify({
      query: request.query.records[0].sequence, target: request.target.records[0].sequence,
      qIdxPos0: request.parameters.qIdxPos0, tIdxPos0: request.parameters.tIdxPos0,
      copomus: request.copomus,
    });
    pyodide.globals.set('browser_request', payload);
    pyodide.globals.set('browser_predict', predict);
    pyodide.globals.set('browser_progress', progress);
    const resultJson = await pyodide.runPythonAsync(`
from copomus.browser import run
await run(browser_request, browser_predict, browser_progress)
`);
    const result = JSON.parse(String(resultJson));
    self.postMessage({ type: 'result', mutation: result, files: [], diagnostics: diagnostics.join('\n') });
  } catch (error) {
    self.postMessage({ type: 'error', message: error instanceof Error ? error.message : String(error), diagnostics: diagnostics.join('\n') });
  }
};
