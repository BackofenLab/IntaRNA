import { parseFasta, prepareRequest, CSV_COLUMNS } from './protocol.js';
import { GROUPS, FIELDS, EXTRA_OUTPUTS } from './options.js';
import { CSV_CATALOG } from './columns.js';
import { PredictionRunner } from './runner.js';
import { ResultsView, downloadText, downloadBlob } from './results.js';
import { createArchive } from './archive.js';
import { importParameters, serializeParameters, MAX_PARAMETER_BYTES } from './parameter-file.js';
import { readFastaFile } from './uploads.js';

const byId = (id) => document.getElementById(id);
const form = byId('prediction-form'), runner = new PredictionRunner();
const dialog = byId('loading'), status = byId('result-status'), errorBox = byId('input-error');
const results = new ResultsView(byId('results-view'));
const csvInputs = [], outputInputs = [];
let completedRun = null;
let inputRevision = 0;
let csvOrder = CSV_COLUMNS.split(',');
function element(tag, text, className) {
  const node = document.createElement(tag);
  if (text) node.textContent = text;
  if (className) node.className = className;
  return node;
}
function activateGroup(key, focus = false) {
  for (const tab of byId('parameter-tabs').children) {
    const active = tab.getAttribute('aria-controls') === key;
    tab.setAttribute('aria-selected', String(active)); tab.tabIndex = active ? 0 : -1;
    byId(tab.getAttribute('aria-controls')).hidden = !active;
    if (active && focus) tab.focus();
  }
}
function addTab(key, label) {
  const tab = element('button', label); tab.type = 'button'; tab.id = `tab-${key}`;
  tab.setAttribute('role', 'tab'); tab.setAttribute('aria-controls', key);
  tab.addEventListener('click', () => activateGroup(key));
  byId('parameter-tabs').append(tab);
}
function group(key, label) {
  const panel = element('section', '', 'parameter-group'); panel.id = key;
  panel.setAttribute('role', 'tabpanel'); panel.setAttribute('aria-labelledby', `tab-${key}`);
  const body = element('div', '', 'group-body'); panel.append(body);
  byId('parameter-groups').append(panel); addTab(key, label); return body;
}
addTab('sequence-input', 'Sequence input');
byId('parameter-tabs').addEventListener('keydown', (event) => {
  if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
  const tabs = [...byId('parameter-tabs').children].filter((tab) => !tab.hidden);
  const current = tabs.indexOf(document.activeElement);
  const next = event.key === 'Home' ? 0 : event.key === 'End' ? tabs.length - 1 : (current + (event.key === 'ArrowRight' ? 1 : -1) + tabs.length) % tabs.length;
  event.preventDefault(); activateGroup(tabs[next].getAttribute('aria-controls'), true);
});
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
  button.addEventListener('click', () => { csvOrder = [...selection]; for (const input of csvInputs) input.checked = selection.includes(input.value); clearResults(); }); selectionTools.append(button);
}
outputBody.append(element('h3', 'Additional output files', 'subheading'));
outputBody.append(element('p', 'Reports are computed with the same settings. ASCII, detailed and ensemble reports each run a fresh engine instance and can take extra time.', 'input-help'));
const outputPicker = element('div', '', 'checkbox-grid'); outputBody.append(outputPicker);
for (const output of EXTRA_OUTPUTS) outputInputs.push(checkbox(outputPicker, output.key, output.label));
const generalBody = group('general-options', 'General');
generalBody.append(element('p', 'Tasks run locally in one worker. Abort or Escape stops computation immediately. Process-level settings (threads, log paths, parameter-file paths and terminal help) are available in the CLI.'));
const docs = element('a', 'Full CLI parameter documentation'); docs.href = 'https://github.com/BackofenLab/IntaRNA#readme'; generalBody.append(docs);

function clearResults(message = 'Inputs changed. Compute to update the interactions.') {
  inputRevision++;
  completedRun = null; results.clear(); byId('output-files').hidden = true; byId('file-list').replaceChildren();
  byId('diagnostics').hidden = true; byId('diagnostic-text').textContent = ''; status.textContent = message;
}
function updateLength(name) {
  if (!byId(name).value.trim()) { byId(`${name}-length`).textContent = 'No sequence entered'; return; }
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
    const tab = byId(`tab-${definition.key}`); tab.hidden = evaluation && !!definition.prediction;
    if (tab.hidden) {
      byId(definition.key).hidden = true;
      if (tab.getAttribute('aria-selected') === 'true') activateGroup('sequence-input');
    }
    for (const field of definition.fields) {
      const input = byId(field.key);
      input.disabled = evaluation && !!field.prediction || !evaluation && definition.key === 'helix-options' && byId('model').value !== 'B' || field.key === 'kineticScore' && byId('mode').value !== 'K';
    }
  }
  for (const input of outputInputs) { input.disabled = evaluation && EXTRA_OUTPUTS.find((entry) => entry.key === input.value).prediction; if (input.disabled) input.checked = false; }
  updateLength('query'); updateLength('target');
}
function resetParameters({ clearInputs = false } = {}) {
  if (clearInputs) {
    byId('task').value = 'prediction';
    byId('query').value = ''; byId('target').value = '';
    for (const input of form.querySelectorAll('input[type=file]')) input.value = '';
  }
  csvOrder = CSV_COLUMNS.split(',');
  for (const field of FIELDS) { const input = byId(field.key); if (field.type === 'checkbox') input.checked = field.value; else input.value = field.value; }
  for (const input of csvInputs) input.checked = CSV_COLUMNS.split(',').includes(input.value);
  for (const input of outputInputs) input.checked = false;
  byId('qIdxPos0').value = 1; byId('tIdxPos0').value = 1; byId('rri').value = '';
  updateAvailability(); errorBox.hidden = true; clearResults(clearInputs ? 'Inputs cleared. Enter sequences or choose an example.' : 'Parameters reset to CLI defaults. Compute to update the interactions.');
}
byId('task').addEventListener('change', resetParameters);
byId('reset').addEventListener('click', () => resetParameters({ clearInputs: true }));
form.addEventListener('input', (event) => {
  if (event.target.type === 'file') return;
  updateAvailability(); errorBox.hidden = true; clearResults();
});
function showInputError(error) {
  errorBox.textContent = error.message; errorBox.hidden = false;
}
function applyParameters(request, message = 'Parameters loaded. Compute to analyze the inputs.') {
  // Import validates every setting before this atomic UI update.
  byId('task').value = request.task;
  resetParameters();
  for (const key of ['query', 'target', 'qIdxPos0', 'tIdxPos0', 'rri']) byId(key).value = request[key];
  for (const field of FIELDS) {
    const value = request.parameters[field.key] ?? field.value;
    if (field.type === 'checkbox') byId(field.key).checked = value;
    else byId(field.key).value = value;
  }
  csvOrder = [...request.csvColumns];
  for (const input of csvInputs) input.checked = request.csvColumns.includes(input.value);
  for (const input of outputInputs) input.checked = request.outputs.includes(input.value);
  updateAvailability(); activateGroup('sequence-input'); clearResults(message);
}
for (const strand of ['query', 'target']) {
  const input = byId(strand), picker = byId(`${strand}-upload`);
  const loadFasta = async (files) => {
    if (byId('inputs').disabled) return;
    const revision = ++inputRevision;
    try {
      if (files.length !== 1) throw new Error('Choose or drop one FASTA file per field.');
      const content = await readFastaFile(files[0], strand === 'query' ? 'Query' : 'Target', byId('task').value);
      if (inputRevision !== revision) return;
      input.value = content;
      updateAvailability(); errorBox.hidden = true;
      clearResults(`${strand === 'query' ? 'Query' : 'Target'} FASTA loaded. Compute to analyze the inputs.`);
    } catch (error) { if (inputRevision === revision) showInputError(error); }
  };
  byId(`${strand}-upload-button`).addEventListener('click', () => picker.click());
  picker.addEventListener('change', () => { const files = [...picker.files]; picker.value = ''; if (files.length) loadFasta(files); });
  input.addEventListener('dragover', (event) => {
    if (event.dataTransfer.types.includes('Files')) { event.preventDefault(); event.dataTransfer.dropEffect = 'copy'; input.classList.add('drag-over'); }
  });
  input.addEventListener('dragleave', () => input.classList.remove('drag-over'));
  input.addEventListener('drop', (event) => {
    input.classList.remove('drag-over');
    if (!event.dataTransfer.types.includes('Files')) return;
    event.preventDefault(); loadFasta([...event.dataTransfer.files]);
  });
}
byId('load-parameters').addEventListener('click', () => byId('parameter-upload').click());
byId('parameter-upload').addEventListener('change', async (event) => {
  const file = event.target.files[0]; event.target.value = ''; if (!file) return;
  const revision = ++inputRevision;
  try {
    if (file.size > MAX_PARAMETER_BYTES) throw new Error('Choose a JSON file of at most 5 MB.');
    const imported = importParameters(await file.text());
    if (inputRevision !== revision) return;
    applyParameters(imported);
  } catch (error) { if (inputRevision === revision) showInputError(error); }
});
byId('store-parameters').addEventListener('click', async () => {
  // Snapshot before waiting for build metadata so edits cannot change the download.
  const request = currentRequest();
  try { downloadText('parameter.json', serializeParameters(request, await buildInfo)); }
  catch (error) { showInputError(error); }
});
async function loadExamples() {
  const container = byId('example-options');
  try {
    const response = await fetch(new URL('./examples/index.json', import.meta.url));
    if (!response.ok) throw new Error('Example list is unavailable.');
    const files = await response.json();
    if (!Array.isArray(files) || !files.length || files.some((name) => typeof name !== 'string' || !/^[a-z0-9-]+\.json$/.test(name))) throw new Error('Invalid example list.');
    const examples = await Promise.all(files.map(async (filename) => {
      const response = await fetch(new URL(`./examples/${filename}`, import.meta.url));
      if (!response.ok) throw new Error('An example is unavailable.');
      const example = await response.json();
      if (typeof example.name !== 'string' || typeof example.description !== 'string') throw new Error('Invalid example description.');
      const request = importParameters(example.configuration);
      prepareRequest(request);
      return { ...example, request };
    }));
    container.replaceChildren();
    for (const example of examples) {
      const button = element('button'); button.type = 'button';
      button.append(element('strong', example.name), element('span', example.description));
      button.addEventListener('click', () => {
        if (byId('inputs').disabled) return;
        applyParameters(example.request, `${example.name} loaded. Compute to analyze the example.`);
        byId('examples-menu').open = false;
      });
      container.append(button);
    }
  } catch (error) { container.replaceChildren(element('p', `${error.message} Reload to try again.`, 'input-help')); }
}
loadExamples();
byId('examples-menu').addEventListener('keydown', (event) => {
  if (event.key === 'Escape') { byId('examples-menu').open = false; byId('examples-menu').querySelector('summary').focus(); }
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
function currentRequest() {
  const selected = csvInputs.filter((input) => input.checked).map((input) => input.value);
  const columns = [...csvOrder.filter((key) => selected.includes(key)), ...selected.filter((key) => !csvOrder.includes(key))];
  return {
    task: byId('task').value, query: byId('query').value, target: byId('target').value,
    qIdxPos0: byId('qIdxPos0').value, tIdxPos0: byId('tIdxPos0').value, rri: byId('rri').value,
    parameters: Object.fromEntries(FIELDS.map((field) => [field.key, field.type === 'checkbox' ? byId(field.key).checked : byId(field.key).value])),
    csvColumns: columns, outputs: outputInputs.filter((input) => input.checked && !input.disabled).map((input) => input.value),
  };
}
form.addEventListener('submit', async (event) => {
  event.preventDefault();
  const request = currentRequest();
  errorBox.hidden = true;
  let prepared;
  try { prepared = prepareRequest(request); } catch (error) { errorBox.textContent = error.message; errorBox.hidden = false; return; }
  clearResults('Computing interactions…'); byId('loading-status').textContent = 'Loading the prediction engine…';
  byId('inputs').disabled = true; byId('results').setAttribute('aria-busy', 'true'); dialog.showModal();
  try {
    const result = await runner.run(request, (message) => { byId('loading-status').textContent = message; });
    const count = results.show(result.csv, result.csvColumns, prepared);
    status.textContent = count ? `${request.task === 'evaluation' ? 'Evaluation' : 'Prediction'} complete: ${count} interaction${count === 1 ? '' : 's'}.` : request.task === 'prediction' && Number(request.parameters.outNumber) === 0 ? 'Accessibility computation complete. No interaction rows requested.' : 'No favorable interactions found with these sequences and settings.';
    completedRun = { request: prepared, result };
    showFiles(result.files); byId('diagnostic-text').textContent = result.diagnostics; byId('diagnostics').hidden = !result.diagnostics;
  } catch (error) {
    if (error.name === 'AbortError') status.textContent = 'Computation aborted. You can edit the sequences and try again.';
    else { status.textContent = 'Computation failed.'; errorBox.textContent = error.message; errorBox.hidden = false; byId('diagnostic-text').textContent = error.diagnostics || ''; byId('diagnostics').hidden = !error.diagnostics; }
  } finally {
    byId('inputs').disabled = false; updateAvailability(); byId('results').setAttribute('aria-busy', 'false'); dialog.close(); byId('compute').focus();
  }
});
activateGroup('sequence-input');
resetParameters({ clearInputs: true });
const buildInfo = fetch(new URL('./build-info.txt', import.meta.url)).then((response) => {
  if (!response.ok) throw new Error('Missing build provenance'); return response.text();
});
// Keep the rejection handled even if no archive is requested.
buildInfo.catch(() => {});
fetch(new URL('./build-info.json', import.meta.url)).then((response) => {
  if (!response.ok) throw new Error('Missing build provenance'); return response.json();
}).then((info) => { byId('version-info').textContent = `IntaRNA ${info.intarna} (${info.commit.slice(0, 12)}${info.dirty ? ', modified sources' : ''}) · ViennaRNA ${info.viennarna} · Emscripten ${info.emscripten}`; }).catch(() => { byId('version-info').textContent = 'Build versions unavailable. Serve the staged interface with its compiled engine and build-info files.'; });
if (!globalThis.Worker || !globalThis.WebAssembly) { byId('compute').disabled = true; errorBox.textContent = 'This browser needs WebAssembly and Web Workers. Use a current browser or the command line version.'; errorBox.hidden = false; }

byId('download-zip').addEventListener('click', async () => {
  const run = completedRun; if (!run) return;
  const button = byId('download-zip'); button.disabled = true;
  try {
    const provenance = await buildInfo;
    if (completedRun !== run) return;
    downloadBlob('intarna-results.zip', new Blob([createArchive(run.request, run.result, provenance)], { type: 'application/zip' }));
  } catch (error) {
    errorBox.textContent = `ZIP download failed: ${error.message}. Reload the staged interface with its compiled version files.`; errorBox.hidden = false;
  } finally { button.disabled = false; }
});
