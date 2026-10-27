import { downloadText } from './results.js';
import { InteractionVisualization } from './visualization.js';
export { COPOMUS_MEASURES, COPOMUS_FILTERS, COPOMUS_DEFAULTS, normalizeCopomus } from './copomus-options.js';

const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
const MUTATION_COLOR = '006400';
const ENERGY_STATES = Object.freeze(['ww', 'wm', 'mw', 'mm']);

function numeric(value) {
  return typeof value === 'number' && Number.isFinite(value) ? value : null;
}

export function mutationAnnotations(candidate) {
  const queryPosition = Number(candidate.values.qIndex), targetPosition = Number(candidate.values.tIndex);
  const mutated = String(candidate.values.bpMutated ?? '');
  if (!Number.isInteger(queryPosition) || queryPosition === 0 || !Number.isInteger(targetPosition) || targetPosition === 0 || !/^[ACGU]{2}$/.test(mutated)) throw new Error('The selected mutation has invalid sequence coordinates or bases.');
  // vaRRI receives native sequence 1 (target) before sequence 2 (query), while
  // CopomuS stores each mutation and base pair in query-first order.
  return `1:${targetPosition}${mutated[1]}:${MUTATION_COLOR},2:${queryPosition}${mutated[0]}:${MUTATION_COLOR}`;
}

export function mutationEnergyPlot(candidate) {
  const energies = ENERGY_STATES.map((state) => numeric(candidate.interactions[state]?.E));
  const wildtype = energies[0];
  return {
    data: [{
      type: 'bar', x: [...ENERGY_STATES], y: energies,
      marker: { color: '#79a9cc', line: { color: '#145a96', width: 1 } },
      hovertemplate: '%{x}: %{y:.2f} kcal/mol<extra></extra>',
    }],
    layout: {
      autosize: true, height: 398, margin: { l: 58, r: 14, t: 18, b: 48 },
      showlegend: false, paper_bgcolor: '#ffffff', plot_bgcolor: '#ffffff',
      font: { family: '"Segoe UI", Tahoma, Geneva, sans-serif', color: '#1c2b39', size: 12 },
      xaxis: { type: 'category', fixedrange: true },
      yaxis: { title: { text: 'E [kcal/mol]' }, fixedrange: true, zerolinecolor: '#ccd9e5', gridcolor: '#e0e8f0' },
      shapes: wildtype === null ? [] : [{
        type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: wildtype, y1: wildtype,
        line: { color: '#145a96', width: 2, dash: 'dot' },
      }],
    },
    config: { displayModeBar: false, responsive: true, staticPlot: false },
  };
}

class MutationVisualization {
  constructor(root) {
    this.root = root;
    this.viewer = new InteractionVisualization(root, 'mutation-');
    this.title = root.querySelector('[data-role=mutation-visualization-title]');
    this.plot = root.querySelector('[data-role=mutation-energy-plot]');
  }

  clear() {
    this.candidate = null; this.request = null;
    if (globalThis.Plotly?.purge) globalThis.Plotly.purge(this.plot);
    this.plot.hidden = true;
    this.viewer.clear(); this.title.textContent = 'Compensatory-mutation visualization';
  }

  show(candidate, request) {
    this.candidate = candidate; this.request = request;
    this.title.textContent = candidate.values.mutation;
    const interaction = candidate.interactions.ww;
    if (!interaction) {
      this.viewer.clear();
      this.viewer.status.textContent = 'The wild type has no favorable reported interaction to visualize.';
    } else {
      const visualRequest = {
        ...request,
        query: { ...request.query, records: [{ ...request.query.records[0], sequence: candidate.queryWildtype }] },
        target: { ...request.target, records: [{ ...request.target.records[0], sequence: candidate.targetWildtype }] },
      };
      try {
        this.viewer.show(interaction, visualRequest, { mutations: mutationAnnotations(candidate) });
        if (!this.viewer.frame.hidden) this.viewer.status.textContent += ' Designed point mutations are annotated in dark green.';
      } catch (error) { this.viewer.clear(); this.viewer.status.textContent = error.message; }
    }
    const plot = mutationEnergyPlot(candidate);
    this.plot.hidden = false;
    this.plot.setAttribute('aria-label', `Interaction energies for ${candidate.values.mutation}: ${ENERGY_STATES.map((state, index) => `${state} ${plot.data[0].y[index] ?? 'not available'}`).join(', ')} kcal/mol`);
    if (globalThis.Plotly?.newPlot) globalThis.Plotly.newPlot(this.plot, plot.data, plot.layout, plot.config);
  }
}

export class CopomusResultsView {
  constructor(root) {
    this.root = root; this.candidates = []; this.key = 'rank'; this.direction = 1;
    this.visualization = new MutationVisualization(root);
    this.search = root.querySelector('[data-role=mutation-search]');
    this.table = root.querySelector('table'); this.header = this.table.querySelector('thead tr'); this.body = this.table.querySelector('tbody');
    this.count = root.querySelector('[data-role=mutation-visible-count]'); this.empty = root.querySelector('[data-role=mutation-no-matches]');
    this.picker = root.querySelector('[data-role=mutation-visible-columns]');
    this.search.addEventListener('input', () => this.render());
    root.querySelector('[data-role=mutation-download]').addEventListener('click', () => downloadText('copomus-mutations.tsv', this.csv));
  }

  clear() {
    this.candidates = []; this.csv = ''; this.selected = null; this.body.replaceChildren(); this.root.hidden = true; this.visualization.clear();
  }

  show(result, request) {
    if (!isObject(result) || typeof result.tsv !== 'string' || !Array.isArray(result.columns) || !Array.isArray(result.candidates) || result.columns.length < 6 || result.candidates.some((candidate) => !isObject(candidate.values) || !isObject(candidate.interactions))) throw new Error('The mutation-design engine returned malformed results.');
    this.csv = result.tsv; this.columns = [...result.columns]; this.candidates = result.candidates; this.request = request;
    this.selected = null; this.visualization.clear(); this.visible = new Set(this.columns);
    this.picker.replaceChildren();
    for (const key of this.columns) {
      const label = document.createElement('label'), input = document.createElement('input');
      input.type = 'checkbox'; input.checked = true; input.value = key;
      input.addEventListener('change', () => { if (input.checked) this.visible.add(key); else this.visible.delete(key); this.render(); });
      label.append(input, document.createTextNode(key)); this.picker.append(label);
    }
    this.search.value = ''; this.key = this.columns.includes('rank') ? 'rank' : this.columns[0]; this.direction = 1;
    this.root.hidden = false; this.render(); return this.candidates.length;
  }

  selectedCandidates() {
    const query = this.search.value.trim().toLowerCase();
    return this.candidates.filter((candidate) => Object.values(candidate.values).some((value) => String(value).toLowerCase().includes(query))).toSorted((left, right) => {
      const a = left.values[this.key], b = right.values[this.key], x = numeric(a), y = numeric(b);
      if (x !== null && y !== null) return this.direction * (x - y);
      if (x === null && y !== null) return 1;
      if (x !== null && y === null) return -1;
      return this.direction * String(a ?? '').localeCompare(String(b ?? ''));
    });
  }

  render() {
    const rows = this.selectedCandidates(), columns = this.columns.filter((key) => this.visible.has(key));
    this.header.replaceChildren();
    const viewHeading = document.createElement('th'); viewHeading.scope = 'col'; viewHeading.setAttribute('aria-label', 'Mutation visualization'); this.header.append(viewHeading);
    for (const key of columns) {
      const cell = document.createElement('th'), button = document.createElement('button'); cell.scope = 'col'; button.type = 'button'; button.textContent = key;
      const active = key === this.key, direction = this.direction === 1 ? 'ascending' : 'descending';
      cell.setAttribute('aria-sort', active ? direction : 'none'); button.setAttribute('aria-label', `Sort by ${key}${active ? `, currently ${direction}` : ''}`);
      button.addEventListener('click', () => { this.direction = active ? -this.direction : 1; this.key = key; this.render(); });
      cell.append(button); this.header.append(cell);
    }
    this.body.replaceChildren();
    for (const candidate of rows) {
      const row = document.createElement('tr'); row.classList.toggle('selected', candidate === this.selected);
      const action = document.createElement('td'), view = document.createElement('button'); view.type = 'button'; view.textContent = 'View'; view.className = 'view-interaction';
      view.setAttribute('aria-label', `Visualize mutation ${candidate.values.mutation}`); view.setAttribute('aria-pressed', String(candidate === this.selected));
      row.addEventListener('click', () => {
        this.selected = candidate; this.visualization.show(candidate, this.request);
        for (const other of this.body.children) { other.classList.toggle('selected', other === row); other.querySelector('.view-interaction').setAttribute('aria-pressed', String(other === row)); }
      });
      action.append(view); row.append(action);
      for (const key of columns) { const cell = document.createElement('td'); cell.textContent = String(candidate.values[key]); if (numeric(candidate.values[key]) !== null) cell.className = 'number'; row.append(cell); }
      this.body.append(row);
    }
    this.count.textContent = `${rows.length} of ${this.candidates.length} mutations`; this.empty.hidden = rows.length > 0;
  }
}
