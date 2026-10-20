import { prepareRequest, makeArguments } from './protocol.js';

// One main() call per worker: C++ globals and MEMFS never cross predictions.
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
    // Validate again at the worker boundary before allocating the WASM runtime.
    const request = prepareRequest(data);
    self.postMessage({ type: 'status', message: 'Loading the prediction engine…' });
    const { default: createIntaRNA } = await import('./wasm/intarna.js');
    const module = await createIntaRNA({ noInitialRun: true, print: log, printErr: log });
    module.FS.writeFile('/query.fa', request.query.fasta);
    module.FS.writeFile('/target.fa', request.target.fasta);
    self.postMessage({ type: 'status', message: 'Predicting RNA interactions…' });
    const code = module.callMain(makeArguments(request.count));
    if (code !== 0) throw new Error(`IntaRNA exited with status ${code}.`);
    const csv = module.FS.readFile('/results.csv', { encoding: 'utf8' });
    self.postMessage({ type: 'result', csv, diagnostics: diagnostics.join('\n') });
  } catch (error) {
    self.postMessage({
      type: 'error',
      message: error instanceof Error ? error.message : String(error),
      diagnostics: diagnostics.join('\n'),
    });
  }
};
