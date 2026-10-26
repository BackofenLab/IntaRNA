const PROFILE_TYPES = new Map([
  ['qMinE', { strand: 'query', value: 'minE', title: 'Minimum-energy profile', yTitle: 'Minimum interaction energy (kcal/mol)' }],
  ['tMinE', { strand: 'target', value: 'minE', title: 'Minimum-energy profile', yTitle: 'Minimum interaction energy (kcal/mol)' }],
  ['qSpotProb', { strand: 'query', value: 'spotProb', title: 'Interaction-probability profile', yTitle: 'Interaction probability' }],
  ['tSpotProb', { strand: 'target', value: 'spotProb', title: 'Interaction-probability profile', yTitle: 'Interaction probability' }],
]);
const MATRIX_TYPES = new Map([
  ['pMinE', { value: 'minE', title: 'Pair minimum energies', colorTitle: 'kcal/mol', colorscale: [[0, '#0f4676'], [0.5, '#8bbbd8'], [1, '#f8fbff']] }],
  ['spotProb', { value: 'spotProb', title: 'Interaction spot probabilities', colorTitle: 'Probability', colorscale: [[0, '#f8fbff'], [0.5, '#80b6d8'], [1, '#145a96']] }],
]);

function outputType(name) {
  const match = /^(qMinE|tMinE|qSpotProb|tSpotProb|pMinE|spotProb)(?:-t\d+q\d+)?\.csv$/.exec(name);
  return match?.[1] ?? null;
}

function table(text) {
  if (typeof text !== 'string' || !text.trim()) throw new Error('The output file is empty.');
  const lines = text.trimEnd().split(/\r?\n/), header = lines.shift().split(';');
  if (header.length < 2 || new Set(header).size !== header.length) throw new Error('The output file has an invalid header.');
  const rows = lines.filter((line) => line.trim()).map((line) => {
    const cells = line.split(';');
    if (cells.length !== header.length) throw new Error('The output file contains a malformed row.');
    return cells;
  });
  if (!rows.length) throw new Error('The output file contains no values to plot.');
  return { header, rows };
}

function number(value) {
  if (/^(NA|NAN)$/i.test(value.trim())) return null;
  const parsed = Number(value);
  if (!Number.isFinite(parsed)) throw new Error(`The output contains a non-numeric value: ${value}`);
  return parsed;
}

function plotText(value) {
  return String(value).replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;');
}

function record(request, strand, index) {
  return request?.[strand]?.records?.[index - 1]?.id ?? '';
}

export function outputFileContext(name, request) {
  const type = outputType(name) ?? /^(qAcc|tAcc|qPu|tPu|bpProb|bpsvg)/.exec(name)?.[1] ?? '';
  const pair = /-t(\d+)q(\d+)(?=\.[^.]+$)/.exec(name);
  if (pair) {
    const target = record(request, 'target', Number(pair[1])), query = record(request, 'query', Number(pair[2]));
    return { target, query, label: [target, query].filter(Boolean).join(' : ') };
  }
  const single = /-s(\d+)(?=\.[^.]+$)/.exec(name);
  if (single && /^[qt](Acc|Pu)/.test(type)) {
    const strand = type.startsWith('q') ? 'query' : 'target', value = record(request, strand, Number(single[1]));
    return { [strand]: value, label: value };
  }
  const target = request?.target?.records?.length === 1 ? request.target.records[0].id : '';
  const query = request?.query?.records?.length === 1 ? request.query.records[0].id : '';
  if (/^[qt](Acc|Pu)$/.test(type)) {
    const value = type.startsWith('q') ? query : target;
    return { [type.startsWith('q') ? 'query' : 'target']: value, label: value };
  }
  return { target, query, label: [target, query].filter(Boolean).join(' : ') };
}

function baseLayout(title) {
  return {
    title: { text: title, x: 0.04, xanchor: 'left', font: { size: 17 } },
    autosize: true,
    height: 420,
    margin: { l: 72, r: 28, t: 62, b: 78 },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: '#fff',
    font: { family: 'Segoe UI, Tahoma, Geneva, sans-serif', color: '#1c2b39', size: 12 },
  };
}

function categoricalAxis(title, values) {
  return { title: { text: title }, type: 'category', categoryorder: 'array', categoryarray: values, automargin: true, gridcolor: '#e0e8f0' };
}

export function outputVisualization(file, request) {
  const type = outputType(file.name);
  if (!type) return null;
  const parsed = table(file.text), context = outputFileContext(file.name, request);
  if (PROFILE_TYPES.has(type)) {
    const profile = PROFILE_TYPES.get(type);
    if (parsed.header.length !== 3 || parsed.header[0] !== 'idx' || parsed.header[2] !== profile.value) throw new Error(`Unexpected ${type} profile columns.`);
    const x = parsed.rows.map((row) => String(row[0])), y = parsed.rows.map((row) => number(row[2]));
    const sequence = parsed.header[1] || context[profile.strand] || (profile.strand === 'query' ? 'Query' : 'Target');
    return {
      kind: 'profile', context,
      traces: [{
        type: 'scatter', mode: 'lines+markers', x, y, text: parsed.rows.map((row) => row[1]),
        name: plotText(sequence), line: { color: '#145a96', width: 2 }, marker: { color: '#145a96', size: 6 },
        connectgaps: false, hovertemplate: 'Index %{x}<br>Nucleotide %{text}<br>Value %{y}<extra></extra>',
      }],
      layout: {
        ...baseLayout(`${profile.title}: ${plotText(sequence)}`),
        xaxis: categoricalAxis(`${profile.strand === 'query' ? 'Query' : 'Target'} index`, x),
        yaxis: { title: { text: profile.yTitle }, automargin: true, gridcolor: '#e0e8f0', rangemode: profile.value === 'spotProb' ? 'tozero' : 'normal' },
      },
    };
  }

  const matrix = MATRIX_TYPES.get(type);
  if (type === 'spotProb' && parsed.header[0] === 'spot') {
    if (parsed.header.length !== 2 || parsed.header[1] !== 'probability') throw new Error('Unexpected selected-spot probability columns.');
    const x = parsed.rows.map((row) => String(row[0])), y = parsed.rows.map((row) => number(row[1]));
    return {
      kind: 'spots', context,
      traces: [{ type: 'bar', x, y, marker: { color: '#145a96' }, hovertemplate: 'Spot %{x}<br>Probability %{y}<extra></extra>' }],
      layout: {
        ...baseLayout(`${matrix.title}${context.label ? `: ${plotText(context.label)}` : ''}`),
        xaxis: categoricalAxis('Target & query index', x),
        yaxis: { title: { text: matrix.colorTitle }, automargin: true, gridcolor: '#e0e8f0', rangemode: 'tozero' },
      },
    };
  }
  if (parsed.header[0] !== matrix.value || parsed.header.length < 2) throw new Error(`Unexpected ${type} matrix columns.`);
  const x = parsed.header.slice(1).map(String), y = parsed.rows.map((row) => String(row[0]));
  const z = parsed.rows.map((row) => row.slice(1).map(number));
  const target = context.target || 'Target', query = context.query || 'Query';
  return {
    kind: 'matrix', context,
    traces: [{
      type: 'heatmap', x, y, z, colorscale: matrix.colorscale, hoverongaps: false,
      colorbar: { title: { text: matrix.colorTitle, side: 'right' } },
      hovertemplate: `${plotText(target)} %{y}<br>${plotText(query)} %{x}<br>Value %{z}<extra></extra>`,
    }],
    layout: {
      ...baseLayout(`${matrix.title}: ${plotText(target)} : ${plotText(query)}`),
      xaxis: categoricalAxis(`${plotText(query)} index`, x),
      yaxis: categoricalAxis(`${plotText(target)} index`, y),
    },
  };
}

function button(text, className = '') {
  const node = document.createElement('button'); node.type = 'button'; node.textContent = text; node.className = className; return node;
}

export class AdditionalOutputView {
  constructor(file, request, id, download) {
    this.file = file; this.id = id; this.root = document.createElement('details'); this.root.className = 'output-file';
    const summary = document.createElement('summary'); summary.className = 'output-file-summary';
    const name = document.createElement('span'); name.className = 'output-file-name'; name.textContent = file.name;
    const downloadButton = button('Download', 'text-button output-download');
    downloadButton.setAttribute('aria-label', `Download ${file.name}`);
    downloadButton.addEventListener('click', (event) => { event.preventDefault(); event.stopPropagation(); download(); });
    const context = outputFileContext(file.name, request), description = document.createElement('span');
    description.className = 'output-file-context'; description.textContent = context.label ? `— ${context.label}` : '';
    summary.append(name, downloadButton, description); this.root.append(summary);

    const body = document.createElement('div'); body.className = 'output-file-body'; this.root.append(body);
    const preview = document.createElement('pre'); preview.textContent = file.text.slice(0, 100000);
    const contentPanel = document.createElement('section'); contentPanel.className = 'output-tab-panel'; contentPanel.append(preview);
    if (file.text.length > 100000) {
      const note = document.createElement('p'); note.textContent = 'Preview limited to 100,000 characters. Download contains the complete file.'; contentPanel.append(note);
    }

    let spec = null, plotError = null;
    try { spec = outputVisualization(file, request); } catch (error) { plotError = error; }
    if (!spec && !plotError) { body.append(contentPanel); return; }

    const tabs = document.createElement('div'); tabs.className = 'output-tabs'; tabs.setAttribute('role', 'tablist'); tabs.setAttribute('aria-label', `${file.name} views`);
    const visualizationTab = button('Visualization', 'output-tab'), contentTab = button('File Content', 'output-tab');
    const visualizationPanel = document.createElement('section'); visualizationPanel.className = 'output-tab-panel';
    visualizationTab.id = `${id}-visualization-tab`; visualizationPanel.id = `${id}-visualization-panel`;
    contentTab.id = `${id}-content-tab`; contentPanel.id = `${id}-content-panel`;
    for (const [tab, panel] of [[visualizationTab, visualizationPanel], [contentTab, contentPanel]]) {
      tab.setAttribute('role', 'tab'); tab.setAttribute('aria-controls', panel.id);
      panel.setAttribute('role', 'tabpanel'); panel.setAttribute('aria-labelledby', tab.id);
    }
    tabs.append(visualizationTab, contentTab); body.append(tabs, visualizationPanel, contentPanel);
    this.tabs = [visualizationTab, contentTab]; this.panels = [visualizationPanel, contentPanel];
    this.selectTab = (index, focus = false) => {
      this.tabs.forEach((tab, current) => { tab.setAttribute('aria-selected', String(current === index)); tab.tabIndex = current === index ? 0 : -1; this.panels[current].hidden = current !== index; });
      if (focus) this.tabs[index].focus();
      if (index === 0) this.render();
    };
    visualizationTab.addEventListener('click', () => this.selectTab(0)); contentTab.addEventListener('click', () => this.selectTab(1));
    tabs.addEventListener('keydown', (event) => {
      if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
      const current = this.tabs.indexOf(document.activeElement);
      const next = event.key === 'Home' ? 0 : event.key === 'End' ? 1 : (current + (event.key === 'ArrowRight' ? 1 : -1) + 2) % 2;
      event.preventDefault(); this.selectTab(next, true);
    });

    const tools = document.createElement('div'); tools.className = 'plot-tools';
    const status = document.createElement('p'); status.className = 'output-plot-status'; status.setAttribute('role', 'status');
    visualizationPanel.append(tools, status); this.status = status;
    if (plotError) { status.textContent = `Visualization unavailable: ${plotError.message}`; this.selectTab(0); return; }
    const exportButton = button('Export PNG', 'secondary'); exportButton.disabled = true;
    exportButton.addEventListener('click', async () => {
      try {
        await this.ready;
        const image = await globalThis.Plotly.toImage(this.plot, { format: 'png', width: 1400, height: 850 });
        const encoded = image.slice(image.indexOf(',') + 1), bytes = Uint8Array.from(atob(encoded), (character) => character.charCodeAt(0));
        const url = URL.createObjectURL(new Blob([bytes], { type: 'image/png' }));
        const link = document.createElement('a'); link.href = url; link.download = file.name.replace(/\.csv$/i, '.png');
        document.body.append(link); link.click(); link.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
      } catch (error) { status.textContent = `PNG export failed: ${error.message}`; }
    });
    tools.append(exportButton); this.exportButton = exportButton;
    this.plot = document.createElement('div'); this.plot.className = 'output-plot'; visualizationPanel.append(this.plot); this.spec = spec;
    this.root.addEventListener('toggle', () => { if (this.root.open && !visualizationPanel.hidden) this.render(); });
    this.selectTab(0);
  }

  render() {
    if (!this.spec || this.rendered || !this.root.open) return;
    if (!globalThis.Plotly?.newPlot) { this.status.textContent = 'Visualization unavailable: Plotly did not load.'; return; }
    this.rendered = true; this.status.textContent = 'Rendering visualization…'; this.plot.setAttribute('aria-busy', 'true');
    this.ready = Promise.resolve(globalThis.Plotly.newPlot(this.plot, this.spec.traces, this.spec.layout, {
      responsive: true, displaylogo: false, scrollZoom: false,
      modeBarButtonsToRemove: ['sendDataToCloud', 'lasso2d', 'select2d'],
    })).then(() => {
      this.plot.removeAttribute('aria-busy'); this.status.textContent = ''; this.exportButton.disabled = false;
    }).catch((error) => { this.plot.removeAttribute('aria-busy'); this.status.textContent = `Visualization unavailable: ${error.message}`; });
  }

  resize() {
    if (this.rendered && this.root.open && !this.panels?.[0].hidden) requestAnimationFrame(() => globalThis.Plotly?.Plots?.resize(this.plot));
  }

  destroy() {
    if (this.rendered) globalThis.Plotly?.purge?.(this.plot);
  }
}
