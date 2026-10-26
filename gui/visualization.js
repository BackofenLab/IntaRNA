import { sequenceOffset } from './protocol.js';

const SEED_HIGHLIGHT = '800080:0.3';

function seedHighlights(row) {
  const values = ['seedStart1', 'seedEnd1', 'seedStart2', 'seedEnd2'].map((key) => row[key]);
  if (values.some((value) => value === undefined || /^(?:NA|NAN)$/i.test(String(value).trim()))) return '';
  const coordinates = values.map((value) => String(value).split(':').map(Number));
  const count = coordinates[0].length;
  if (!count || coordinates.some((list) => list.length !== count || list.some((value) => !Number.isInteger(value) || value === 0))) return '';
  const highlights = [];
  for (let index = 0; index < count; index++) {
    const [start1, end1, start2, end2] = coordinates.map((list) => list[index]);
    if (start1 > end1 || start2 > end2) return '';
    highlights.push(`${start1}-${end1}&${start2}-${end2}:${SEED_HIGHLIGHT}`);
  }
  return highlights.join(',');
}

export function visualizationUrl(row, request) {
  // Native seq1 is target and seq2 is query; both structures run 5′ to 3′.
  const structures = row.hybridDP?.split('&');
  if (structures?.length !== 2) throw new Error('The interaction has no dot-bracket structure.');
  const sequences = [], fullStructures = [], starts = [];
  for (const [index, strand] of ['target', 'query'].entries()) {
    const number = index + 1, start = request.parameters[index ? 'qIdxPos0' : 'tIdxPos0'];
    const record = request[strand].records.find((record) => record.id.replaceAll(';', '_') === row[`id${number}`]);
    if (!record) throw new Error('The interaction sequence was not found in this run.');
    if (start === 0) throw new Error('vaRRI requires nonzero sequence start indices. Choose a nonzero start index to visualize this pair.');
    const offset = sequenceOffset(row[`start${number}`], start), structure = structures[index];
    const end = sequenceOffset(row[`end${number}`], start);
    if (offset < 0 || end >= record.sequence.length || structure.length !== end - offset + 1 || /[^.()]/.test(structure)) throw new Error('The interaction structure does not match its sequence coordinates.');
    sequences.push(record.sequence); starts.push(start);
    fullStructures.push('.'.repeat(offset) + structure + '.'.repeat(record.sequence.length - end - 1));
  }
  const url = new URL('https://backofenlab.github.io/vaRRI/'), parameters = new URLSearchParams({
    sequence: sequences.join('&'), structure: fullStructures.join('&'),
    seqName1: row.id1, seqName2: row.id2, startIndex1: starts[0], startIndex2: starts[1],
    showRenderingOnly: 'true', forceLayout: 'false', forceLayoutLinearRRI: 'true',
    cropping: 2
  });
  const highlights = seedHighlights(row);
  if (highlights) parameters.set('regionHighlights', highlights);
  url.search = parameters;
  return url.href;
}

export class InteractionVisualization {
  constructor(root) {
    this.frame = root.querySelector('[data-role=varri-frame]'); this.link = root.querySelector('[data-role=varri-link]');
    this.status = root.querySelector('[data-role=visualization-status]');
  }
  clear() {
    this.frame.removeAttribute('src'); this.frame.hidden = true;
    this.link.removeAttribute('href'); this.link.hidden = true;
    this.status.textContent = 'Select View on an interaction to visualize it.';
  }
  show(row, request) {
    this.clear();
    try {
      const url = visualizationUrl(row, request);
      this.frame.src = url; this.frame.hidden = false;
      const fullUrl = new URL(url); fullUrl.searchParams.delete('showRenderingOnly');
      this.link.href = fullUrl.href; this.link.hidden = false;
      this.status.textContent = `${row.id1} (${row.start1}–${row.end1}) & ${row.id2} (${row.start2}–${row.end2}), ${row.E.toFixed(2)} kcal/mol. If the viewer is unavailable, use Open in vaRRI to retry.`;
    } catch (error) { this.status.textContent = error.message; }
  }
}
