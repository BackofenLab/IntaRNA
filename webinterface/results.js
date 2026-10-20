import { CSV_COLUMNS } from './protocol.js';

export const COLUMNS = Object.freeze([
  { key: 'id2', label: 'Query', numeric: false },
  { key: 'start2', label: 'Query start', numeric: true },
  { key: 'end2', label: 'Query end', numeric: true },
  { key: 'id1', label: 'Target', numeric: false },
  { key: 'start1', label: 'Target start', numeric: true },
  { key: 'end1', label: 'Target end', numeric: true },
  { key: 'E', label: 'Energy (kcal/mol)', numeric: true },
  { key: 'deltaE', label: 'ΔE (kcal/mol)', numeric: true },
]);

// IntaRNA emits unquoted semicolon records and replaces delimiters in IDs.
// Header control characters are rejected before prediction; quotes are literal.
export function parseResults(csv) {
  if (typeof csv !== 'string') throw new Error('The prediction did not return CSV text.');
  const lines = csv.trim().split(/\r?\n/);
  if (!csv.trim()) return [];
  if (lines.shift() !== CSV_COLUMNS.replaceAll(',', ';')) {
    throw new Error('Unexpected IntaRNA CSV columns. Rebuild and serve the matching web interface and engine.');
  }
  const keys = CSV_COLUMNS.split(',');
  const rows = lines.filter((line) => line.trim()).map((line) => {
    const values = line.split(';');
    if (values.length !== keys.length) throw new Error('Malformed IntaRNA CSV row.');
    const row = Object.fromEntries(keys.map((key, index) => [key, values[index]]));
    for (const key of ['start1', 'end1', 'start2', 'end2', 'E']) {
      if (!row[key].trim() || !Number.isFinite(Number(row[key]))) throw new Error(`Invalid ${key} in prediction output.`);
      row[key] = Number(row[key]);
    }
    for (const key of ['start1', 'end1', 'start2', 'end2']) {
      if (!Number.isInteger(row[key]) || row[key] < 1) throw new Error('Invalid interaction coordinates.');
    }
    if (row.start1 > row.end1 || row.start2 > row.end2) throw new Error('Reversed interaction coordinates.');
    return row;
  });
  const best = Math.min(...rows.map((row) => row.E));
  return rows.map((row) => ({ ...row, deltaE: Number((row.E - best).toFixed(2)) }));
}

export function selectRows(rows, search, key, direction) {
  const query = search.toLowerCase().trim();
  const column = COLUMNS.find((entry) => entry.key === key);
  return rows.filter((row) => COLUMNS.some((entry) => String(row[entry.key]).toLowerCase().includes(query)))
    .toSorted((a, b) => direction * (column.numeric ? a[key] - b[key] : a[key].localeCompare(b[key])));
}

export class ResultsView {
  constructor(root) {
    this.root = root;
    this.rows = [];
    this.key = 'E';
    this.direction = 1;
    this.search = root.querySelector('#search');
    this.table = root.querySelector('table');
    this.body = this.table.querySelector('tbody');
    this.count = root.querySelector('#visible-count');
    this.empty = root.querySelector('#no-matches');
    this.search.addEventListener('input', () => this.render());
    const header = this.table.querySelector('thead tr');
    for (const column of COLUMNS) {
      const cell = document.createElement('th');
      cell.scope = 'col';
      const button = document.createElement('button');
      button.type = 'button';
      button.textContent = column.label;
      button.addEventListener('click', () => {
        this.direction = this.key === column.key ? -this.direction : 1;
        this.key = column.key;
        this.render();
      });
      cell.append(button);
      header.append(cell);
    }
    root.querySelector('#download').addEventListener('click', () => {
      const url = URL.createObjectURL(new Blob([this.csv], { type: 'text/csv;charset=utf-8' }));
      const link = document.createElement('a');
      link.href = url;
      link.download = 'intarna-interactions.csv';
      link.click();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    });
  }

  clear() {
    this.rows = [];
    this.csv = '';
    this.body.replaceChildren();
    this.root.hidden = true;
  }

  show(csv) {
    this.rows = parseResults(csv);
    this.csv = csv;
    this.search.value = '';
    this.key = 'E';
    this.direction = 1;
    this.root.hidden = this.rows.length === 0;
    this.render();
    return this.rows.length;
  }

  render() {
    const visible = selectRows(this.rows, this.search.value, this.key, this.direction);
    this.body.replaceChildren();
    for (const row of visible) {
      const tr = document.createElement('tr');
      for (const column of COLUMNS) {
        const cell = document.createElement('td');
        cell.textContent = ['E', 'deltaE'].includes(column.key) ? row[column.key].toFixed(2) : String(row[column.key]);
        if (column.numeric) cell.className = 'number';
        tr.append(cell);
      }
      this.body.append(tr);
    }
    this.table.querySelectorAll('th').forEach((cell, index) => {
      const active = COLUMNS[index].key === this.key;
      cell.setAttribute('aria-sort', active ? (this.direction === 1 ? 'ascending' : 'descending') : 'none');
      cell.firstChild.setAttribute('aria-label', `Sort by ${COLUMNS[index].label}${active ? `, currently ${this.direction === 1 ? 'ascending' : 'descending'}` : ''}`);
    });
    this.count.textContent = `${visible.length} of ${this.rows.length} interactions`;
    this.empty.hidden = visible.length > 0;
  }
}
