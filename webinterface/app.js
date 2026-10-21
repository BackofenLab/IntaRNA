import { DEFAULTS, EVALUATION_EXAMPLE, parseFasta, prepareRequest, CSV_COLUMNS } from './protocol.js';
import { GROUPS, FIELDS, EXTRA_OUTPUTS } from './options.js';
import { CSV_CATALOG } from './columns.js';
import { PredictionRunner } from './runner.js';
import { ResultsView, downloadText } from './results.js';

const byId = (id) => document.getElementById(id);
const form = byId('prediction-form'), runner = new PredictionRunner();
const dialog = byId('loading'), status = byId('result-status'), errorBox = byId('input-error');
const results = new ResultsView(byId('results-view'));
const csvInputs = [], outputInputs = [];
function element(tag, text, className) {
  const node = document.createElement(tag);
  if (text) node.textContent = text;
  if (className) node.className = className;
  return node;
}
function group(key, label) {
  const details = element('details', '', 'parameter-group');
  details.id = key; details.name = 'parameters'; details.append(element('summary', label));
  const body = element('div', '', 'group-body'); details.append(body);
  byId('parameter-groups').append(details); return body;
}
for (const definition of GROUPS) {
  const body = group(definition.key, definition.label), grid = element('div', '', 'field-grid'); body.append(grid);
  for (const field of definition.fields) {
    const wrapper = element('div', '', `field${field.type === 'checkbox' ? ' check' : ''}`);
    const label = element('label', field.label); label.htmlFor = field.key; label.title = `--${field.key}`;
    const input = element(field.type === 'select' ? 'select' : field.type === 'textarea' ? 'textarea' : 'input'); input.id = field.key;
    if (input.tagName === 'INPUT') input.type = field.type;
    if (field.type === 'select') for (const [value, name] of Object.entries(field.choices)) { const option = element('option', name); option.value = value; input.append(option); }
    if (field.type === 'number') { input.min = field.min; input.max = field.max; input.step = field.step ?? 1; }
    if (field.type === 'checkbox') { input.checked = field.value; wrapper.append(input, label); }
    else { input.value = field.value; wrapper.append(label, input); }
    if (field.hint) { const hint = element('p', field.hint, 'field-hint'); hint.id = `${field.key}-help`; input.setAttribute('aria-describedby', hint.id); wrapper.append(hint); }
    if (field.type === 'textarea') {
      input.spellcheck = false;
      const upload = element('input'); upload.type = 'file'; upload.accept = '.txt,.csv,.shape,.dat'; upload.setAttribute('aria-label', `Load ${field.label} from a local text file`);
      upload.dataset.for = field.key;
      upload.addEventListener('change', async () => {
        const file = upload.files[0]; if (!file) return;
        if (file.size > 1000000) { errorBox.textContent = 'Choose a text file of at most 1 MB.'; errorBox.hidden = false; return; }
        input.value = await file.text(); input.dispatchEvent(new Event('input', { bubbles: true })); upload.value = '';
      }); wrapper.append(upload);
    }
    grid.append(wrapper);
  }
}
function checkbox(container, key, labelText, description) {
  const label = element('label'), input = element('input'); input.type = 'checkbox'; input.value = key;
  if (description) label.title = description;
  label.append(input, document.createTextNode(labelText)); container.append(label); return input;
}
const outputBody = byId('output-options').querySelector('.group-body');
outputBody.append(element('h3', 'CSV output columns', 'subheading'));
outputBody.append(element('p', 'seq1 is the target; seq2 is the query. Select columns for both the table and CSV download. You can hide table columns after computation without changing the CSV.', 'input-help'));
const selectionTools = element('div', '', 'selection-tools'); outputBody.append(selectionTools);
const csvPicker = element('div', '', 'checkbox-grid'); outputBody.append(csvPicker);
for (const column of [...CSV_COLUMNS.split(',').map((key) => CSV_CATALOG.find((column) => column.key === key)), ...CSV_CATALOG.filter((column) => !CSV_COLUMNS.split(',').includes(column.key))]) csvInputs.push(checkbox(csvPicker, column.key, `${column.key} — ${column.description}`, column.description));
for (const [name, selection] of [['Select all', CSV_CATALOG.map((c) => c.key)], ['CLI defaults', CSV_COLUMNS.split(',')], ['Clear selection', []]]) {
  const button = element('button', name); button.type = 'button';
  button.addEventListener('click', () => { for (const input of csvInputs) input.checked = selection.includes(input.value); clearResults(); }); selectionTools.append(button);
}
outputBody.append(element('h3', 'Additional output files', 'subheading'));
outputBody.append(element('p', 'Reports are computed with the same settings. ASCII, detailed and ensemble reports each run a fresh engine instance and can take extra time.', 'input-help'));
const outputPicker = element('div', '', 'checkbox-grid'); outputBody.append(outputPicker);
for (const output of EXTRA_OUTPUTS) outputInputs.push(checkbox(outputPicker, output.key, output.label));
const generalBody = group('general-options', 'General');
generalBody.append(element('p', 'Tasks run locally in one worker. Abort or Escape stops computation immediately. Process-level settings (threads, log paths, parameter-file paths and terminal help) are available in the CLI.'));
const docs = element('a', 'Full CLI parameter documentation'); docs.href = 'https://github.com/BackofenLab/IntaRNA#readme'; generalBody.append(docs);

// Use explicit grouping as a fallback for browsers without details[name].
for (const details of document.querySelectorAll('.parameter-group')) {
  details.addEventListener('toggle', () => { if (details.open) for (const other of document.querySelectorAll('.parameter-group')) if (other !== details) other.open = false; });
}
function clearResults(message = 'Inputs changed. Compute to update the interactions.') {
  results.clear(); byId('output-files').hidden = true; byId('file-list').replaceChildren();
  byId('diagnostics').hidden = true; byId('diagnostic-text').textContent = ''; status.textContent = message;
}
function updateLength(name) {
  try {
    const parsed = parseFasta(byId(name).value, name === 'query' ? 'Query' : 'Target', byId('task').value === 'prediction');
    byId(`${name}-length`).textContent = `${parsed.records.length} record${parsed.records.length === 1 ? '' : 's'}, ${parsed.records.reduce((sum, r) => sum + r.sequence.length, 0).toLocaleString()} nucleotides`;
  } catch { byId(`${name}-length`).textContent = 'Check sequence format and length'; }
}
function updateAvailability() {
  const evaluation = byId('task').value === 'evaluation';
  byId('evaluation-input').hidden = !evaluation;
  byId('task-help').textContent = evaluation ? 'Evaluate supplied structures for one query and one target. Prediction constraints are ignored.' : 'Predict interactions for all query–target sequence pairs.';
  byId('compute').textContent = evaluation ? 'Evaluate interactions' : 'Compute interactions';
  for (const definition of GROUPS) {
    const details = byId(definition.key); details.hidden = evaluation && !!definition.prediction;
    for (const field of definition.fields) {
      const input = byId(field.key);
      input.disabled = evaluation && !!field.prediction || !evaluation && definition.key === 'helix-options' && byId('model').value !== 'B' || field.key === 'kineticScore' && byId('mode').value !== 'K';
    }
  }
  for (const input of outputInputs) { input.disabled = evaluation && EXTRA_OUTPUTS.find((entry) => entry.key === input.value).prediction; if (input.disabled) input.checked = false; }
  updateLength('query'); updateLength('target');
}
function resetParameters() {
  for (const field of FIELDS) { const input = byId(field.key); if (field.type === 'checkbox') input.checked = field.value; else input.value = field.value; }
  for (const input of csvInputs) input.checked = CSV_COLUMNS.split(',').includes(input.value);
  for (const input of outputInputs) input.checked = false;
  byId('qIdxPos0').value = 1; byId('tIdxPos0').value = 1; byId('rri').value = EVALUATION_EXAMPLE;
  updateAvailability(); errorBox.hidden = true; clearResults('Parameters reset to CLI defaults. Compute to update the interactions.');
}
byId('task').addEventListener('change', resetParameters);
byId('reset').addEventListener('click', resetParameters);
form.addEventListener('input', (event) => {
  if (event.target.type === 'file') return;
  updateAvailability(); errorBox.hidden = true; clearResults();
});
byId('example').addEventListener('click', () => {
  byId('query').value = DEFAULTS.query; byId('target').value = DEFAULTS.target;
  resetParameters(); byId('sequence-input').open = true; clearResults('Example loaded. Compute to analyze its interactions.');
});
byId('abort').addEventListener('click', () => runner.abort());
dialog.addEventListener('cancel', (event) => { event.preventDefault(); runner.abort(); });
window.addEventListener('pagehide', () => runner.abort());
function showFiles(files) {
  for (const file of files) {
    const details = element('details'), summary = element('summary', file.name), preview = element('pre', file.text.slice(0, 100000));
    const download = element('button', 'Download file', 'secondary'); download.type = 'button'; download.addEventListener('click', () => downloadText(file.name, file.text));
    details.append(summary, download, preview);
    if (file.text.length > 100000) details.append(element('p', 'Preview limited to 100,000 characters. Download contains the complete file.'));
    byId('file-list').append(details);
  }
  byId('output-files').hidden = !files.length;
}
form.addEventListener('submit', async (event) => {
  event.preventDefault();
  const request = {
    task: byId('task').value, query: byId('query').value, target: byId('target').value,
    qIdxPos0: byId('qIdxPos0').value, tIdxPos0: byId('tIdxPos0').value, rri: byId('rri').value,
    parameters: Object.fromEntries(FIELDS.map((field) => [field.key, field.type === 'checkbox' ? byId(field.key).checked : byId(field.key).value])),
    csvColumns: csvInputs.filter((input) => input.checked).map((input) => input.value), outputs: outputInputs.filter((input) => input.checked && !input.disabled).map((input) => input.value),
  };
  errorBox.hidden = true;
  try { prepareRequest(request); } catch (error) { errorBox.textContent = error.message; errorBox.hidden = false; return; }
  clearResults('Computing interactions…'); byId('loading-status').textContent = 'Loading the prediction engine…';
  byId('inputs').disabled = true; byId('results').setAttribute('aria-busy', 'true'); dialog.showModal();
  try {
    const result = await runner.run(request, (message) => { byId('loading-status').textContent = message; });
    const count = results.show(result.csv, result.csvColumns);
    status.textContent = count ? `${request.task === 'evaluation' ? 'Evaluation' : 'Prediction'} complete: ${count} interaction${count === 1 ? '' : 's'}.` : request.task === 'prediction' && Number(request.parameters.outNumber) === 0 ? 'Accessibility computation complete. No interaction rows requested.' : 'No favorable interactions found with these sequences and settings.';
    showFiles(result.files); byId('diagnostic-text').textContent = result.diagnostics; byId('diagnostics').hidden = !result.diagnostics;
  } catch (error) {
    if (error.name === 'AbortError') status.textContent = 'Computation aborted. You can edit the sequences and try again.';
    else { status.textContent = 'Computation failed.'; errorBox.textContent = error.message; errorBox.hidden = false; byId('diagnostic-text').textContent = error.diagnostics || ''; byId('diagnostics').hidden = !error.diagnostics; }
  } finally {
    byId('inputs').disabled = false; updateAvailability(); byId('results').setAttribute('aria-busy', 'false'); dialog.close(); byId('compute').focus();
  }
});
resetParameters();
fetch(new URL('./build-info.json', import.meta.url)).then((response) => {
  if (!response.ok) throw new Error('Missing build provenance'); return response.json();
}).then((info) => { byId('version-info').textContent = `IntaRNA ${info.intarna} (${info.commit.slice(0, 12)}${info.dirty ? ', modified sources' : ''}) · ViennaRNA ${info.viennarna} · Emscripten ${info.emscripten}`; }).catch(() => { byId('version-info').textContent = 'Build versions unavailable. Serve the staged interface with its compiled engine and build-info files.'; });
if (!globalThis.Worker || !globalThis.WebAssembly) { byId('compute').disabled = true; errorBox.textContent = 'This browser needs WebAssembly and Web Workers. Use a current browser or the command line version.'; errorBox.hidden = false; }
