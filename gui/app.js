import { parseFasta, prepareRequest, CSV_COLUMNS } from './protocol.js';
import { GROUPS, FIELDS, EXTRA_OUTPUTS, UI_DEFAULTS } from './options.js';
import { CSV_CATALOG, CSV_GROUPS, defaultColumns } from './columns.js';
import { PredictionRunner } from './runner.js';
import { downloadText } from './results.js';
import { SessionsView } from './sessions.js';
import { importParameters, serializeParameters, MAX_PARAMETER_BYTES } from './parameter-file.js';
import { readFastaFile } from './uploads.js';

const byId = (id) => document.getElementById(id);
const form = byId('prediction-form'), runner = new PredictionRunner();
const getTask = () => form.querySelector('input[name=task]:checked').value;
const setTask = (task) => { for (const input of form.querySelectorAll('input[name=task]')) input.checked = input.value === task; };
const dialog = byId('loading'), status = byId('result-status'), errorBox = byId('input-error');
const sessions = new SessionsView(byId('sessions'), byId('session-template'), (request, name) => {
  applyParameters(request, `${name} selected. Its computed inputs have been restored.`);
}, showInputError, (message) => { status.textContent = message; });
const csvInputs = [], outputInputs = [];
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
addTab('sequence-input', 'Sequences');
byId('parameter-tabs').addEventListener('keydown', (event) => {
  if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
  const tabs = [...byId('parameter-tabs').children].filter((tab) => !tab.hidden);
  const current = tabs.indexOf(document.activeElement);
  const next = event.key === 'Home' ? 0 : event.key === 'End' ? tabs.length - 1 : (current + (event.key === 'ArrowRight' ? 1 : -1) + tabs.length) % tabs.length;
  event.preventDefault(); activateGroup(tabs[next].getAttribute('aria-controls'), true);
});
// Render each requested row explicitly; labels share a grid row so wrapped
// labels do not push neighboring controls out of alignment.
function renderField(field, id = field.key) {
  const wrapper = element('div', '', `field${field.type === 'checkbox' ? ' check' : ''}`);
  const label = element('label', field.label); label.htmlFor = id; label.title = `--${field.key}`;
  const input = element(field.type === 'select' ? 'select' : field.type === 'textarea' ? 'textarea' : 'input'); input.id = id;
  if (input.tagName === 'INPUT') input.type = field.type;
  if (field.type === 'select') for (const [value, name] of Object.entries(field.choices)) { const option = element('option', name); option.value = value; input.append(option); }
  if (field.type === 'number') { input.min = field.min; input.max = field.max; input.step = field.step ?? 1; }
  if (field.type === 'checkbox') { input.checked = field.value; wrapper.append(input, label); }
  else { input.value = field.value; wrapper.append(label, input); }
  if (field.hint) { const hint = element('p', field.hint, 'field-hint'); hint.id = `${id}-help`; input.setAttribute('aria-describedby', hint.id); wrapper.append(hint); }
  if (field.type === 'textarea') {
    input.spellcheck = false;
    const upload = element('input'); upload.type = 'file'; upload.accept = '.txt,.csv,.shape,.dat'; upload.setAttribute('aria-label', `Load ${field.label} from a local text file`);
    upload.dataset.for = field.key;
    upload.addEventListener('change', async () => {
      const file = upload.files[0]; upload.value = ''; if (!file) return;
      const revision = ++inputRevision;
      try {
        if (file.size > 1000000) throw new Error('Choose a text file of at most 1 MB.');
        const content = await file.text();
        if (inputRevision !== revision || input.disabled) return;
        input.value = content; input.dispatchEvent(new Event('input', { bubbles: true }));
      } catch (error) { if (inputRevision === revision) showInputError(error); }
    }); wrapper.append(upload);
  }
  if (field.key === 'noSeed') wrapper.hidden = true;
  return wrapper;
}
function fieldRow(body, fields, suffix = '') {
  const grid = element('div', '', `field-grid columns-${fields.length}`);
  for (const field of fields) grid.append(renderField(field, field.key + suffix));
  body.append(grid); return grid;
}
for (const definition of GROUPS) {
  const body = group(definition.key, definition.label);
  if (definition.key === 'seed-options') {
    const primary = element('div', '', 'seed-primary-row');
    const wrapper = element('div', '', 'field seed-mode');
    const label = element('label', 'Seed constraint'); label.htmlFor = 'seed-mode';
    const select = element('select'); select.id = 'seed-mode';
    for (const [value, name] of [['none', 'Disable seed constraint'], ['explicit', 'Explicit seed'], ['identify', 'Identify seeds']]) {
      const option = element('option', name); option.value = value; select.append(option);
    }
    select.value = 'identify'; wrapper.append(label, select);
    primary.append(wrapper, renderField(FIELDS.find((field) => field.key === 'seedTQ')));
    body.append(primary);
  }
  if (definition.rows) {
    for (const row of definition.rows) {
      if (definition.key === 'seed-options' && row[0].key === 'seedTQ') continue;
      fieldRow(body, row);
    }
  } else fieldRow(body, definition.fields);
  if (definition.strand) {
    const p = definition.strand;
    const shared = FIELDS.filter((field) => ['accNoLP', 'accNoGUend'].includes(field.key));
    fieldRow(body, shared, p === 't' ? '-target' : '');
    body.append(element('p', 'These two accessibility settings apply to both strands.', 'input-help'));
    const columns = element('div', '', 'accessibility-columns'); body.append(columns);
    const accessibility = element('div'), shape = element('div'); columns.append(accessibility, shape);
    const field = (key) => FIELDS.find((entry) => entry.key === p + key);
    fieldRow(accessibility, [field('AccConstr')]);
    fieldRow(accessibility, [field('AccFile')]);
    fieldRow(shape, [field('Shape')]);
    fieldRow(shape, [field('ShapeMethod'), field('ShapeConversion')]);
  }
}
function checkbox(container, key, labelText, description) {
  const label = element('label'), input = element('input'); input.type = 'checkbox'; input.value = key;
  if (description) label.title = description;
  label.append(input, document.createTextNode(labelText)); container.append(label); return input;
}
const csvBody = byId('csv-options').querySelector('.group-body');
csvBody.append(element('h3', 'CSV output columns', 'subheading'));
csvBody.append(element('p', 'seq1 is the target; seq2 is the query. Hover over a column name for details. Select columns for the table and CSV download; table visibility can be changed after computation.', 'input-help'));
const selectionTools = element('div', '', 'selection-tools'); csvBody.append(selectionTools);
const selectColumns = (selection) => {
  csvOrder = [...selection]; for (const input of csvInputs) input.checked = selection.includes(input.value); inputsChanged();
};
for (const [name, selection] of [['Select all', () => CSV_CATALOG.map((c) => c.key)], ['Task defaults', () => defaultColumns(getTask())], ['Clear selection', () => []]]) {
  const button = element('button', name); button.type = 'button';
  button.addEventListener('click', () => selectColumns(selection())); selectionTools.append(button);
}
for (const group of CSV_GROUPS) {
  const section = element('fieldset', '', 'csv-group'), legend = element('legend', group.label);
  section.dataset.csvGroup = group.label;
  for (const [name, checked] of [['Select all', true], ['Deselect all', false]]) {
    const button = element('button', name, 'text-button'); button.type = 'button';
    button.setAttribute('aria-label', `${name} ${group.label} columns`);
    button.addEventListener('click', () => {
      for (const input of csvInputs.filter((input) => group.keys.includes(input.value))) input.checked = checked;
      inputsChanged();
    }); legend.append(button);
  }
  const picker = element('div', '', 'checkbox-grid csv-columns'); section.append(legend, picker); csvBody.append(section);
  for (const key of group.keys) {
    const column = CSV_CATALOG.find((column) => column.key === key);
    csvInputs.push(checkbox(picker, key, key, column.description));
  }
}
const outputBody = byId('output-options').querySelector('.group-body');
outputBody.append(element('h3', 'Additional output files', 'subheading'));
outputBody.append(element('p', 'Reports are computed with the same settings. ASCII, detailed and ensemble reports each run a fresh engine instance and can take extra time.', 'input-help'));
const outputPicker = element('div', '', 'checkbox-grid'); outputBody.append(outputPicker);
for (const output of EXTRA_OUTPUTS) outputInputs.push(checkbox(outputPicker, output.key, output.label, output.description));
const generalBody = group('general-options', 'CLI-only');
generalBody.append(element('p', 'Process-level settings (threads, log paths, parameter-file paths and terminal help) are available in the CLI.'));
const docs = element('a', 'Full CLI parameter documentation'); docs.href = 'https://github.com/BackofenLab/IntaRNA#readme'; generalBody.append(docs);

function inputsChanged(message = 'Inputs changed. Compute to update the interactions.') {
  inputRevision++;
  sessions.deselect();
  byId('diagnostics').hidden = true; byId('diagnostic-text').textContent = ''; status.textContent = message;
}
function updateLength(name) {
  const warning = byId(`${name}-warning`); warning.hidden = true; warning.textContent = '';
  if (!byId(name).value.trim()) { byId(`${name}-length`).textContent = 'No sequence entered'; return; }
  try {
    const parsed = parseFasta(byId(name).value, name === 'query' ? 'Query' : 'Target', getTask() === 'prediction');
    warning.textContent = parsed.warnings.join(' '); warning.hidden = !parsed.warnings.length;
    byId(`${name}-length`).textContent = `${parsed.records.length} record${parsed.records.length === 1 ? '' : 's'}, ${parsed.records.reduce((sum, r) => sum + r.sequence.length, 0).toLocaleString()} nucleotides`;
  } catch { byId(`${name}-length`).textContent = 'Check sequence format and length'; }
}
function updateAvailability() {
  const evaluation = getTask() === 'evaluation';
  for (const button of byId('example-options').querySelectorAll('[data-example-task]')) {
    button.hidden = button.dataset.exampleTask !== getTask();
  }
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
      input.disabled = evaluation && !!field.prediction || field.key.startsWith('helix') && byId('model').value !== 'B' || field.key === 'kineticScore' && byId('mode').value !== 'K';
    }
  }
  byId('noSeed').checked = byId('seed-mode').value === 'none';
  for (const field of FIELDS.filter((field) => field.key.startsWith('seed'))) {
    byId(field.key).disabled = evaluation || (field.key === 'seedTQ' ? byId('seed-mode').value !== 'explicit' : byId('seed-mode').value !== 'identify');
  }
  for (const p of ['q', 't']) {
    const mode = byId(`${p}Acc`).value;
    for (const key of ['AccW', 'AccL', 'PfScale', 'AccConstr', 'AccFile', 'Shape', 'ShapeMethod', 'ShapeConversion']) {
      byId(p + key).disabled = mode === 'N' || ['AccW', 'AccL'].includes(key) && mode !== 'C' || key === 'AccFile' && !['P', 'E'].includes(mode);
    }
    for (const key of ['accNoLP', 'accNoGUend']) byId(key + (p === 't' ? '-target' : '')).disabled = mode === 'N';
  }
  // Upload controls follow the text area they populate, including async reads.
  for (const upload of form.querySelectorAll('input[data-for]')) upload.disabled = byId(upload.dataset.for).disabled;
  for (const input of outputInputs) { input.disabled = evaluation && EXTRA_OUTPUTS.find((entry) => entry.key === input.value).prediction; if (input.disabled) input.checked = false; }
  updateLength('query'); updateLength('target');
}
function resetParameters({ clearInputs = false } = {}) {
  if (clearInputs) {
    setTask('prediction');
    byId('query').value = ''; byId('target').value = '';
    for (const input of form.querySelectorAll('input[type=file]')) input.value = '';
    byId('examples-menu').open = false; activateGroup('sequence-input');
  }
  csvOrder = defaultColumns(getTask());
  byId('seed-mode').value = 'identify';
  for (const field of FIELDS) { const input = byId(field.key); if (field.type === 'checkbox') input.checked = UI_DEFAULTS[field.key]; else input.value = UI_DEFAULTS[field.key]; }
  syncAccessibilityFlags();
  for (const input of csvInputs) input.checked = csvOrder.includes(input.value);
  for (const input of outputInputs) input.checked = false;
  byId('qIdxPos0').value = 1; byId('tIdxPos0').value = 1; byId('rri').value = '';
  updateAvailability(); errorBox.hidden = true; inputsChanged(clearInputs ? 'Inputs cleared. Enter sequences or choose an example.' : 'Parameters reset to task defaults. Compute to update the interactions.');
}
byId('task').addEventListener('change', resetParameters);
byId('reset').addEventListener('click', () => resetParameters({ clearInputs: true }));
function syncAccessibilityFlags(source) {
  for (const key of ['accNoLP', 'accNoGUend']) {
    if (source?.id === key + '-target') byId(key).checked = source.checked;
    byId(key + '-target').checked = byId(key).checked;
  }
}
form.addEventListener('input', (event) => {
  if (event.target.type === 'file') return;
  syncAccessibilityFlags(event.target);
  updateAvailability(); errorBox.hidden = true; inputsChanged();
});
function showInputError(error) {
  errorBox.textContent = error.message; errorBox.hidden = false;
}
function applyParameters(request, message = 'Parameters loaded. Compute to analyze the inputs.') {
  // Import validates every setting before this atomic UI update.
  setTask(request.task);
  resetParameters();
  for (const key of ['query', 'target', 'qIdxPos0', 'tIdxPos0', 'rri']) byId(key).value = request[key];
  for (const field of FIELDS) {
    const value = request.parameters[field.key] ?? field.value;
    if (field.type === 'checkbox') byId(field.key).checked = value;
    else byId(field.key).value = value;
  }
  byId('seed-mode').value = request.parameters.noSeed ? 'none' : request.parameters.seedTQ ? 'explicit' : 'identify';
  syncAccessibilityFlags();
  csvOrder = [...request.csvColumns];
  for (const input of csvInputs) input.checked = request.csvColumns.includes(input.value);
  for (const input of outputInputs) input.checked = request.outputs.includes(input.value);
  updateAvailability(); activateGroup('sequence-input'); inputsChanged(message);
}
for (const strand of ['query', 'target']) {
  const input = byId(strand), picker = byId(`${strand}-upload`);
  const loadFasta = async (files) => {
    if (byId('inputs').disabled) return;
    const revision = ++inputRevision;
    try {
      if (files.length !== 1) throw new Error('Choose or drop one FASTA file per field.');
      const content = await readFastaFile(files[0], strand === 'query' ? 'Query' : 'Target', getTask());
      if (inputRevision !== revision) return;
      input.value = content;
      updateAvailability(); errorBox.hidden = true;
      inputsChanged(`${strand === 'query' ? 'Query' : 'Target'} FASTA loaded. Compute to analyze the inputs.`);
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
      button.dataset.exampleTask = example.request.task;
      button.append(element('strong', example.name), element('span', example.description));
      button.addEventListener('click', () => {
        if (byId('inputs').disabled) return;
        applyParameters(example.request, `${example.name} loaded. Compute to analyze the example.`);
        byId('examples-menu').open = false;
      });
      container.append(button);
    }
    updateAvailability();
  } catch (error) { container.replaceChildren(element('p', `${error.message} Reload to try again.`, 'input-help')); }
}
loadExamples();
byId('examples-menu').addEventListener('keydown', (event) => {
  if (event.key === 'Escape') { byId('examples-menu').open = false; byId('examples-menu').querySelector('summary').focus(); }
});
byId('abort').addEventListener('click', () => runner.abort());
dialog.addEventListener('cancel', (event) => { event.preventDefault(); runner.abort(); });
window.addEventListener('pagehide', () => runner.abort());
function currentParameters() {
  // Keep inactive edits in the form for toggling back, but exclude them from
  // validation, saved JSON and native computation. Shared flags stay active
  // whenever either strand uses accessibility.
  return Object.fromEntries(FIELDS.map((field) => {
    const input = byId(field.key);
    const inactive = ['accNoLP', 'accNoGUend'].includes(field.key)
      ? input.disabled && byId(field.key + '-target').disabled : input.disabled;
    return [field.key, inactive ? field.value : field.type === 'checkbox' ? input.checked : input.value];
  }));
}
function currentRequest() {
  const selected = csvInputs.filter((input) => input.checked).map((input) => input.value);
  const columns = [...csvOrder.filter((key) => selected.includes(key)), ...selected.filter((key) => !csvOrder.includes(key))];
  return {
    task: getTask(), query: byId('query').value, target: byId('target').value,
    qIdxPos0: byId('qIdxPos0').value, tIdxPos0: byId('tIdxPos0').value, rri: byId('rri').value,
    parameters: currentParameters(),
    csvColumns: columns, outputs: outputInputs.filter((input) => input.checked && !input.disabled).map((input) => input.value),
  };
}
form.addEventListener('submit', async (event) => {
  event.preventDefault();
  const request = currentRequest();
  errorBox.hidden = true;
  let prepared;
  try {
    if (getTask() === 'prediction' && byId('seed-mode').value === 'explicit' && !request.parameters.seedTQ.trim()) throw new Error('Explicit seed: enter target&query dot-bar encodings.');
    prepared = prepareRequest(request); } catch (error) { errorBox.textContent = error.message; errorBox.hidden = false; return; }
  inputsChanged('Computing interactions…'); byId('loading-status').textContent = 'Loading the prediction engine…';
  byId('inputs').disabled = true; byId('results').setAttribute('aria-busy', 'true'); dialog.showModal();
  try {
    const result = await runner.run(request, (message) => { byId('loading-status').textContent = message; });
    // Keep the executed inputs and outputs together, independent of later edits.
    const provenance = await buildInfo.catch(() => 'Compiled version information is unavailable.');
    status.textContent = sessions.add(prepared, result, provenance);
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
