import { prepareRequest, makeParameterDocument } from './protocol.js';

// Each task owns a worker; each CLI invocation owns a fresh WASM instance.
let started = false;
self.onmessage = async ({ data }) => {
  if (started) return;
  started = true;
  const diagnostics = [];
  const log = (line) => {
    diagnostics.push(String(line).slice(0, 2000));
    if (diagnostics.length > 25) diagnostics.shift();
  };
  try {
    const request = prepareRequest(data);
    self.postMessage({ type: 'status', message: 'Loading the prediction engine…' });
    const { default: createIntaRNA } = await import('./wasm/intarna.js');
    const run = async (format) => {
      const module = await createIntaRNA({
        noInitialRun: true,
        locateFile: () => new URL('./wasm/intarna.wasm', import.meta.url).href,
        print: log, printErr: log,
      });
      module.FS.writeFile('/parameter.json', JSON.stringify(makeParameterDocument(request, format)));
      module.FS.mkdir('/output');
      self.postMessage({ type: 'status', message: request.task === 'evaluation' ? 'Evaluating supplied RNA interactions…' : 'Predicting RNA interactions…' });
      const code = module.callMain(['--parameterFile=/parameter.json']);
      if (code !== 0) throw new Error(`IntaRNA exited with status ${code}. Check the diagnostics and input settings.`);
      return module.FS;
    };
    const fs = await run('C');
    const csv = fs.readFile('/results.csv', { encoding: 'utf8' });
    const files = fs.readdir('/output').filter((name) => name !== '.' && name !== '..').map((name) => ({ name, text: fs.readFile(`/output/${name}`, { encoding: 'utf8' }) }));
    for (const [key, format] of [['normal', 'N'], ['detailed', 'D'], ['ensemble', 'E']]) {
      if (!request.outputs.includes(key)) continue;
      const reportFS = await run(format);
      files.push({ name: `intarna-${key}.txt`, text: reportFS.readFile('/report.txt', { encoding: 'utf8' }) });
    }
    self.postMessage({ type: 'result', csv, csvColumns: request.csvColumns, files, diagnostics: diagnostics.join('\n') });
  } catch (error) {
    self.postMessage({ type: 'error', message: error instanceof Error ? error.message : String(error), diagnostics: diagnostics.join('\n') });
  }
};
