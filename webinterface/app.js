import { DEFAULTS, parseSequence, prepareRequest } from './protocol.js';
import { PredictionRunner } from './runner.js';

const byId = (id) => document.getElementById(id);
const form = byId('prediction-form');
const runner = new PredictionRunner();
const dialog = byId('loading');
const status = byId('result-status');
const errorBox = byId('input-error');

function clearResults(message = 'Inputs changed. Compute to update the interactions.') {
  byId('csv-preview').hidden = true;
  byId('csv-preview').textContent = '';
  byId('diagnostics').hidden = true;
  status.textContent = message;
}

function updateLength(name) {
  try {
    const rna = parseSequence(byId(name).value, name === 'query' ? 'Query' : 'Target');
    byId(`${name}-length`).textContent = `${rna.sequence.length.toLocaleString()} nucleotides`;
  } catch {
    byId(`${name}-length`).textContent = 'Check sequence format and length';
  }
}

function setBusy(busy) {
  byId('inputs').disabled = busy;
  byId('results').setAttribute('aria-busy', String(busy));
  if (busy) dialog.showModal();
  else {
    dialog.close();
    byId('compute').focus();
  }
}

for (const name of ['query', 'target', 'count']) {
  byId(name).addEventListener('input', () => {
    if (name !== 'count') updateLength(name);
    errorBox.hidden = true;
    clearResults();
  });
}
byId('example').addEventListener('click', () => {
  for (const name of ['query', 'target', 'count']) byId(name).value = DEFAULTS[name];
  updateLength('query');
  updateLength('target');
  errorBox.hidden = true;
  clearResults('Example loaded. Compute to predict its interactions.');
});
byId('abort').addEventListener('click', () => runner.abort());
dialog.addEventListener('cancel', (event) => { event.preventDefault(); runner.abort(); });
window.addEventListener('pagehide', () => runner.abort());

form.addEventListener('submit', async (event) => {
  event.preventDefault();
  const request = { query: byId('query').value, target: byId('target').value, count: Number(byId('count').value) };
  errorBox.hidden = true;
  try {
    prepareRequest(request);
  } catch (error) {
    errorBox.textContent = error.message;
    errorBox.hidden = false;
    return;
  }
  clearResults('Computing interactions…');
  byId('loading-status').textContent = 'Loading the prediction engine…';
  setBusy(true);
  try {
    const result = await runner.run(request, (message) => { byId('loading-status').textContent = message; });
    byId('csv-preview').textContent = result.csv;
    byId('csv-preview').hidden = false;
    status.textContent = 'Prediction complete. Coordinates are 1-based and inclusive.';
    byId('diagnostic-text').textContent = result.diagnostics;
    byId('diagnostics').hidden = !result.diagnostics;
  } catch (error) {
    if (error.name === 'AbortError') status.textContent = 'Computation aborted. You can edit the sequences and try again.';
    else {
      status.textContent = 'Prediction failed.';
      errorBox.textContent = error.message;
      errorBox.hidden = false;
      byId('diagnostic-text').textContent = error.diagnostics || '';
      byId('diagnostics').hidden = !error.diagnostics;
    }
  } finally {
    setBusy(false);
  }
});

if (!globalThis.Worker || !globalThis.WebAssembly) {
  byId('compute').disabled = true;
  errorBox.textContent = 'This browser needs WebAssembly and Web Workers. Use a current browser or the command line version.';
  errorBox.hidden = false;
}
