import { JSON_SCHEMA, CSV_COLUMNS, normalizeParameters } from './protocol.js';
import { FIELDS, PARAMETER_DEFAULTS, EXTRA_OUTPUTS } from './options.js';
import { CSV_CATALOG } from './columns.js';

export const MAX_PARAMETER_BYTES = 5000000;
const reports = { normal: 'N', detailed: 'D', ensemble: 'E' };
const embeddedKeys = ['qShape', 'tShape', 'qAccFile', 'tAccFile'];
const object = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
function text(value, label, limit = 1000000) {
  if (typeof value !== 'string' || value.length > limit || value.includes('\0')) throw new Error(`${label}: invalid or oversized text.`);
  return value;
}
function selections(request) {
  const columns = request.csvColumns ?? CSV_COLUMNS.split(','), outputs = request.outputs ?? [];
  for (const [values, catalog, label] of [[columns, CSV_CATALOG, 'CSV columns'], [outputs, EXTRA_OUTPUTS.filter((entry) => request.task !== 'evaluation' || !entry.prediction), 'output files']]) {
    if (!Array.isArray(values) || new Set(values).size !== values.length || values.some((key) => !catalog.some((entry) => entry.key === key))) throw new Error(`Choose available ${label}, without duplicates.`);
  }
  return { columns, outputs };
}

// Save editable settings, even before sequences have been entered. Computation
// performs full sequence-dependent validation later through prepareRequest.
export function exportParameters(request, provenance) {
  if (!provenance?.trim()) throw new Error('Compiled version information is unavailable');
  const normalized = normalizeParameters(request);
  const { columns, outputs } = selections(request);
  const parameters = {
    query: { text: text(typeof request.query === 'string' ? request.query : request.query.fasta, 'Query', 60000) },
    target: { text: text(typeof request.target === 'string' ? request.target : request.target.fasta, 'Target', 60000) },
  };
  if (request.task === 'evaluation') parameters.personality = 'IntaRNAeval';
  for (const [key, value] of Object.entries(normalized)) {
    if (key === 'outSpotProbSpots' || value === '' || value === (PARAMETER_DEFAULTS[key] ?? 1)) continue;
    parameters[key] = embeddedKeys.includes(key) ? { text: value } : value;
  }
  if (request.task === 'evaluation' && request.rri) parameters.rri = text(request.rri, 'Interactions', 10000);
  parameters.outMode = 'C';
  if (columns.join(',') !== CSV_COLUMNS) parameters.outCsvCols = columns.join(',');
  const nativeOutputs = outputs.filter((key) => !Object.hasOwn(reports, key)).map((key) => {
    const spots = key === 'spotProb' && normalized.outSpotProbSpots ? `${normalized.outSpotProbSpots}:` : '';
    return `${key}:${spots}${key}${key === 'bpsvg' ? '.svg' : '.csv'}`;
  });
  if (nativeOutputs.length) parameters.out = nativeOutputs;
  const selectedReports = outputs.filter((key) => Object.hasOwn(reports, key));
  const comments = [
    ...provenance.trimEnd().split(/\r?\n/),
    'Reproduce: IntaRNA --parameterFile=parameter.json --out=intarna-interactions.csv',
    ...selectedReports.map((key) => `Additional report: IntaRNA --parameterFile=parameter.json --outMode=${reports[key]} --outCsvCols=${CSV_COLUMNS} --outCsvSort '' --out=intarna-${key}.txt`),
  ];
  const gui = {};
  if (selectedReports.length) gui.reports = selectedReports;
  if (normalized.outSpotProbSpots) gui.outSpotProbSpots = normalized.outSpotProbSpots;
  return { schema: JSON_SCHEMA, _comment: comments, parameters, ...(Object.keys(gui).length ? { gui } : {}) };
}

export function serializeParameters(request, provenance) {
  return JSON.stringify(exportParameters(request, provenance), null, 2) + '\n';
}

// Import only controls the browser can represent. Native-only settings produce
// an actionable error instead of being silently discarded or reaching WASM.
export function importParameters(source) {
  if (typeof source === 'string' && source.length > MAX_PARAMETER_BYTES) throw new Error('Choose a JSON file of at most 5 MB.');
  const document = typeof source === 'string' ? JSON.parse(source.replace(/^\uFEFF/, '')) : source;
  if (!object(document) || document.schema !== JSON_SCHEMA || !object(document.parameters)) throw new Error(`Expected a ${JSON_SCHEMA} parameter document.`);
  for (const key of Object.keys(document)) if (!['schema', '_comment', 'parameters', 'gui'].includes(key)) throw new Error(`Unknown JSON document field: ${key}`);
  const request = { task: 'prediction', query: '', target: '', qIdxPos0: 1, tIdxPos0: 1, rri: '', parameters: {}, outputs: [], csvColumns: CSV_COLUMNS.split(',') };
  const known = new Set(FIELDS.map((field) => field.key));
  for (const [key, value] of Object.entries(document.parameters)) {
    if (value === null) throw new Error(`${key}: null is not a parameter value.`);
    if (['query', 'target', ...embeddedKeys].includes(key)) {
      const content = object(value) && Object.keys(value).length === 1 && typeof value.text === 'string' ? value.text : value;
      // Plain nucleotide strings are native inline sequences; file paths need
      // explicit embedding because the browser cannot open local paths.
      if (typeof content !== 'string' || (content === value && (embeddedKeys.includes(key) || !/^[ACGUTNacgutn]*$/.test(value)))) throw new Error(`${key}: embed the input using {"text":"file contents"}.`);
      if (key === 'query' || key === 'target') request[key] = text(content, key, 60000);
      else request.parameters[key] = text(content, key);
    } else if (key === 'personality') {
      if (!['IntaRNA', 'IntaRNAeval'].includes(value)) throw new Error('The browser supports IntaRNA prediction and IntaRNAeval evaluation.');
      request.task = value === 'IntaRNAeval' ? 'evaluation' : 'prediction';
    } else if (key === 'qIdxPos0' || key === 'tIdxPos0') request[key] = value;
    else if (key === 'rri') request.rri = text(value, key, 10000);
    else if (key === 'outMode') { if (value !== 'C') throw new Error('Use CSV outMode C and gui.reports for additional browser reports.'); }
    else if (key === 'outCsvCols') request.csvColumns = text(value, key).split(',').filter(Boolean);
    else if (key === 'out') {
      const entries = Array.isArray(value) ? value : [value];
      for (const entry of entries) {
        const match = /^([A-Za-z]+):(?:(-?[\d&,\-]+):)?([^:]+)$/.exec(text(entry, 'out'));
        if (!match || !EXTRA_OUTPUTS.some((output) => output.key === match[1] && !Object.hasOwn(reports, output.key))) throw new Error('The browser supports named additional outputs; omit main output paths from the parameter file.');
        request.outputs.push(match[1]);
        if (match[2]) {
          if (match[1] !== 'spotProb') throw new Error('Only spotProb accepts explicit probability spots.');
          request.parameters.outSpotProbSpots = match[2];
        }
      }
    } else if (known.has(key)) request.parameters[key] = value;
    else throw new Error(`The browser cannot represent parameter ${key}; use the native CLI.`);
  }
  if (document.gui !== undefined) {
    if (!object(document.gui)) throw new Error('gui must be an object.');
    for (const [key, value] of Object.entries(document.gui)) {
      if (key === 'reports' && Array.isArray(value) && value.every((entry) => Object.hasOwn(reports, entry))) request.outputs.push(...value);
      else if (key === 'outSpotProbSpots') {
        if (request.parameters.outSpotProbSpots && request.parameters.outSpotProbSpots !== value) throw new Error('Conflicting spotProb selections.');
        request.parameters.outSpotProbSpots = text(value, key, 10000);
      } else throw new Error(`Unsupported GUI setting: ${key}`);
    }
  }
  selections(request);
  // Keep prediction-only controls from being silently ignored on import.
  if (request.task === 'evaluation' && FIELDS.some((field) => field.prediction && Object.hasOwn(request.parameters, field.key))) throw new Error('Evaluation documents must omit prediction-only parameters.');
  const normalized = normalizeParameters(request);
  request.qIdxPos0 = normalized.qIdxPos0; request.tIdxPos0 = normalized.tIdxPos0;
  request.parameters = normalized;
  return request;
}
