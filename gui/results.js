import { CSV_CATALOG } from './columns.js';
import { CONTEXT_COLUMNS } from './protocol.js';
import { InteractionVisualization } from './visualization.js';

export const COLUMNS = Object.freeze([...CSV_CATALOG, { key: 'deltaE', label: 'mfe Δ (kcal/mol)', numeric: true }]);
const catalog = new Map(COLUMNS.map((column) => [column.key, column]));

// Native CSV is unquoted, semicolon-separated; unavailable optional fields use
// NAN and seed fields may contain colon-separated lists. Preserve those strings.
export function parseResults(csv) {
  if (typeof csv !== 'string') throw new Error('The engine did not return CSV text.');
  if (!csv.trim()) return [];
  const lines = csv.trimEnd().split(/\r?\n/), keys = lines.shift().split(';');
  if (new Set(keys).size !== keys.length || keys.some((key) => !catalog.has(key)) || CONTEXT_COLUMNS.some((key) => !keys.includes(key))) throw new Error('Unexpected IntaRNA CSV columns. Serve the matching interface and engine.');
  const rows = lines.filter((line) => line.trim()).map((line) => {
    const values = line.split(';');
    if (values.length !== keys.length) throw new Error('Malformed IntaRNA CSV row.');
    const row = Object.fromEntries(keys.map((key, index) => [key, values[index]]));
    for (const key of CONTEXT_COLUMNS.filter((key) => key !== 'id1' && key !== 'id2')) {
      if (!row[key].trim() || !Number.isFinite(Number(row[key]))) throw new Error(`Invalid ${key} in interaction output.`);
      row[key] = Number(row[key]);
    }
    for (const key of ['start1', 'end1', 'start2', 'end2']) {
      if (!Number.isInteger(row[key])) throw new Error('Invalid interaction coordinates.');
    }
    if (row.start1 > row.end1 || row.start2 > row.end2) throw new Error('Reversed interaction coordinates.');
    return row;
  });
  // Each sequence pair has its own optimum in a multi-FASTA screen.
  const best = new Map();
  for (const row of rows) {
    const pair = JSON.stringify([row.id1, row.id2]);
    best.set(pair, Math.min(best.get(pair) ?? Infinity, row.E));
  }
  return rows.map((row) => ({ ...row, deltaE: Number((row.E - best.get(JSON.stringify([row.id1, row.id2]))).toFixed(2)) }));
}
function numeric(value) { return String(value).trim() && Number.isFinite(Number(value)) ? Number(value) : null; }
export function selectRows(rows, search, key, direction) {
  const query = search.toLowerCase().trim(), column = catalog.get(key);
  return rows.filter((row) => Object.values(row).some((value) => String(value).toLowerCase().includes(query))).toSorted((a, b) => {
    if (column?.numeric) {
      const x = numeric(a[key]), y = numeric(b[key]);
      if (x !== null && y !== null) return direction * (x - y);
      if (x === null && y !== null) return 1;
      if (x !== null && y === null) return -1;
    }
    return direction * String(a[key] ?? '').localeCompare(String(b[key] ?? ''));
  });
}
export function projectCsv(csv, columns) {
  const lines = csv.trimEnd().split(/\r?\n/), keys = lines[0].split(';');
  const indices = columns.map((key) => keys.indexOf(key));
  if (indices.some((index) => index < 0)) throw new Error('Selected columns are absent from the engine output.');
  return lines.map((line) => { const cells = line.split(';'); return indices.map((index) => cells[index]).join(';'); }).join('\n') + '\n';
}
export function downloadText(name, text) {
  downloadBlob(name, new Blob([text], { type: 'text/plain;charset=utf-8' }));
}
export function downloadBlob(name, blob) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url; link.download = name; link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
export class ResultsView {
  constructor(root) {
    this.root = root; this.rows = []; this.key = 'E'; this.direction = 1;
    this.visualization = new InteractionVisualization(root);
    this.search = root.querySelector('[data-role=search]');
    this.table = root.querySelector('table');
    this.body = this.table.querySelector('tbody');
    this.header = this.table.querySelector('thead tr');
    this.count = root.querySelector('[data-role=visible-count]');
    this.empty = root.querySelector('[data-role=no-matches]');
    this.picker = root.querySelector('[data-role=visible-columns]');
    this.search.addEventListener('input', () => this.render());
    root.querySelector('[data-role=download]').addEventListener('click', () => downloadText('intarna-interactions.csv', this.csv));
  }
  clear() {
    this.rows = []; this.csv = ''; this.body.replaceChildren(); this.root.hidden = true;
    this.selected = null; this.visualization.clear();
  }
  show(csv, columns = csv.trim().split('\n')[0].split(';'), request) {
    this.rows = parseResults(csv); this.csv = projectCsv(csv, columns);
    this.request = request; this.selected = null; this.visualization.clear();
    this.columns = columns.map((key) => catalog.get(key));
    if (columns.includes('E')) this.columns.push(catalog.get('deltaE'));
    this.visible = new Set(this.columns.map((c) => c.key));
    this.picker.replaceChildren();
    for (const column of this.columns) {
      const label = document.createElement('label'), input = document.createElement('input');
      input.type = 'checkbox'; input.checked = true; input.value = column.key;
      input.addEventListener('change', () => {
        if (input.checked) this.visible.add(column.key); else this.visible.delete(column.key);
        this.render();
      });
      label.append(input, document.createTextNode(column.label)); this.picker.append(label);
    }
    this.search.value = ''; this.key = columns.includes('E') ? 'E' : columns[0]; this.direction = 1;
    this.root.hidden = false; this.render(true); return this.rows.length;
  }
  render(selectFirst = false) {
    const rows = selectRows(this.rows, this.search.value, this.key, this.direction);
    const columns = this.columns.filter((column) => this.visible.has(column.key));
    // The completed run opens on the same interaction shown at the top of the
    // initial table; later sorting and filtering retain that selection.
    if (selectFirst && rows.length) {
      this.selected = rows[0];
      this.visualization.show(this.selected, this.request);
    }
    this.header.replaceChildren();
    const viewHeading = document.createElement('th'); viewHeading.scope = 'col'; viewHeading.setAttribute('aria-label', 'Interaction visualization');
    this.header.append(viewHeading);
    for (const column of columns) {
      const cell = document.createElement('th'), button = document.createElement('button');
      cell.scope = 'col'; button.type = 'button';
      const split = column.numeric ? column.label.indexOf(' ') : -1;
      if (split > 0) button.append(document.createTextNode(column.label.slice(0, split)), document.createElement('br'), document.createTextNode(column.label.slice(split + 1)));
      else button.textContent = column.label;
      cell.title = button.title = column.description ?? column.label;
      const active = column.key === this.key, direction = this.direction === 1 ? 'ascending' : 'descending';
      cell.setAttribute('aria-sort', active ? direction : 'none');
      button.setAttribute('aria-label', `Sort by ${column.label}${active ? `, currently ${direction}` : ''}`);
      button.addEventListener('click', () => { this.direction = active ? -this.direction : 1; this.key = column.key; this.render(); });
      cell.append(button); this.header.append(cell);
    }
    this.body.replaceChildren();
    for (const row of rows) {
      const tr = document.createElement('tr');
      tr.classList.toggle('selected', row === this.selected);
      const action = document.createElement('td'), view = document.createElement('button');
      view.type = 'button'; view.textContent = 'View'; view.className = 'view-interaction';
      view.setAttribute('aria-label', `Visualize ${row.id1} ${row.start1}–${row.end1} with ${row.id2} ${row.start2}–${row.end2}`);
      view.setAttribute('aria-pressed', String(row === this.selected));
      tr.addEventListener('click', () => {
        this.selected = row; this.visualization.show(row, this.request);
        for (const other of this.body.children) {
          other.classList.toggle('selected', other === tr);
          other.querySelector('.view-interaction').setAttribute('aria-pressed', String(other === tr));
        }
      });
      action.append(view); tr.append(action);
      for (const column of columns) {
        const cell = document.createElement('td');
        cell.textContent = ['E', 'deltaE'].includes(column.key) ? Number(row[column.key]).toFixed(2) : String(row[column.key]);
        if (column.numeric) cell.className = 'number'; tr.append(cell);
      }
      this.body.append(tr);
    }
    this.count.textContent = `${rows.length} of ${this.rows.length} interactions`;
    this.empty.hidden = rows.length > 0;
  }
}
