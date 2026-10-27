import { ResultsView, downloadBlob } from './results.js';
import { serializeParameters, importParameters } from './parameter-file.js';
import { zipFiles } from './archive.js';
import { AdditionalOutputView } from './output-visualization.js';
import { CopomusResultsView } from './copomus.js';

// Sessions live only in this page. Each owns a JSON snapshot and its own view;
// editing the form or exploring another session cannot change an earlier run.
export class SessionsView {
  constructor(root, template, restore, reportError, announce) {
    this.root = root; this.template = template; this.restore = restore; this.reportError = reportError; this.announce = announce;
    this.entries = new Map(); this.nextNumber = 1; this.selected = null;
  }

  deselect() {
    this.selected = null;
    for (const session of this.entries.values()) {
      session.root.dataset.selected = 'false';
      session.toggle.textContent = 'Show';
    }
  }

  setExpanded(session, expanded) {
    session.panel.hidden = !expanded;
    session.toggle.setAttribute('aria-expanded', String(expanded));
    session.toggle.textContent = expanded ? 'Hide' : 'Show';
    // Unload third-party viewers while folded. The table selection is retained
    // and can be visualized again without recomputing the prediction. A newly
    // rendered result already has its initial viewer URL and must not reload it.
    if (!expanded) session.view.visualization.clear();
    else {
      if (session.view.selected && !session.view.visualization.frame.hasAttribute('src')) session.view.visualization.show(session.view.selected, session.view.request);
      for (const fileView of session.fileViews) fileView.resize();
    }
  }

  select(session, restore = true) {
    if (restore) this.restore(importParameters(session.parameterJson), session.name);
    this.selected = session.id;
    for (const entry of this.entries.values()) {
      entry.root.dataset.selected = String(entry === session);
      this.setExpanded(entry, entry === session);
    }
  }

  add(request, result, provenance) {
    // Prepared indices live in parameters, whereas the JSON serializer accepts
    // editable requests with top-level indices. Preserve signed coordinates.
    const parameterJson = serializeParameters({ ...request, qIdxPos0: request.parameters.qIdxPos0, tIdxPos0: request.parameters.tIdxPos0 }, provenance);
    const root = this.template.content.firstElementChild.cloneNode(true);
    const find = (role) => root.querySelector(`[data-role=${role}]`);
    const id = `session-${this.nextNumber}`, task = request.task === 'evaluation' ? 'Evaluation' : request.task === 'mutation' ? 'Mutation Design' : 'Prediction';
    const mutation = request.task === 'mutation';
    const session = {
      id, name: `${task} ${this.nextNumber}`, parameterJson, root,
      panel: find('panel'), toggle: find('toggle'),
      view: mutation ? new CopomusResultsView(find('mutation-results-view')) : new ResultsView(find('results-view')),
      fileViews: [],
    };
    const count = mutation ? session.view.show(result.mutation, request) : session.view.show(result.csv, result.csvColumns, request);
    const summary = mutation ? `Mutation Design complete: ${count} ranked mutation${count === 1 ? '' : 's'}.` : count ? `${task} complete: ${count} interaction${count === 1 ? '' : 's'}.` : request.task === 'prediction' && request.parameters.outNumber === 0 ? 'Accessibility computation complete. No interaction rows requested.' : 'No favorable interactions found with these sequences and settings.';
    find('summary').textContent = summary;
    root.dataset.sessionId = id;
    session.panel.id = `${id}-panel`; session.toggle.setAttribute('aria-controls', session.panel.id);
    const name = find('name'); name.value = session.name;
    const rename = () => {
      session.name = name.value.trim() || `${task} ${id.slice('session-'.length)}`;
      root.setAttribute('aria-label', session.name);
      session.toggle.setAttribute('aria-label', `Show or hide results and restore inputs for ${session.name}`);
      find('delete').setAttribute('aria-label', `Delete session ${session.name}`);
      find('download-zip').setAttribute('aria-label', `Download ZIP for ${session.name}`);
    };
    name.addEventListener('input', rename);
    name.addEventListener('change', () => { rename(); name.value = session.name; });
    name.addEventListener('keydown', (event) => { if (event.key === 'Enter') name.blur(); });
    rename();
    // The name and header background also select a run. Action buttons retain
    // their own behavior, so exporting/deleting cannot overwrite edited inputs.
    const header = root.querySelector('.session-header');
    header.tabIndex = 0; header.setAttribute('role', 'group');
    header.setAttribute('aria-label', 'Select result and restore its inputs');
    for (const event of ['click', 'focusin']) header.addEventListener(event, (event) => {
      if (!event.target.closest('button') && (this.selected !== id || session.panel.hidden)) this.select(session);
    });
    header.addEventListener('keydown', (event) => {
      if (event.target === header && ['Enter', ' '].includes(event.key)) { event.preventDefault(); this.select(session); }
    });
    // Clicking or keyboard-focusing an output area selects its computed inputs.
    for (const event of ['click', 'focusin']) session.panel.addEventListener(event, () => {
      if (this.selected !== id) this.select(session);
    });
    session.toggle.addEventListener('click', () => {
      if (this.selected === id && !session.panel.hidden) this.setExpanded(session, false);
      else this.select(session);
    });

    // Each archive uses the saved JSON and full outputs, never the current form,
    // filtered table, or the session selected at download time.
    find('download-zip').addEventListener('click', () => {
      try {
        const files = [{ name: 'parameter.json', text: parameterJson }, { name: mutation ? 'copomus-mutations.tsv' : 'intarna-interactions.csv', text: session.view.csv }, ...result.files];
        const filename = session.name.replace(/[^\p{L}\p{N}._-]+/gu, '-').replace(/^[.-]+|[.-]+$/g, '').slice(0, 80) || id;
        downloadBlob(`${filename}.zip`, new Blob([zipFiles(files)], { type: 'application/zip' }));
      } catch (error) { this.reportError(new Error(`ZIP download failed: ${error.message}`)); }
    });
    for (const [index, file] of result.files.entries()) {
      const fileView = new AdditionalOutputView(file, request, `${id}-file-${index}`, () => {
        downloadBlob(file.name, new Blob([file.text], { type: 'text/plain;charset=utf-8' }));
      });
      session.fileViews.push(fileView); find('file-list').append(fileView.root);
    }
    find('output-files').hidden = !result.files.length;
    find('diagnostic-text').textContent = result.diagnostics; find('diagnostics').hidden = !result.diagnostics;

    find('delete').addEventListener('click', () => {
      const entries = [...this.entries.values()], index = entries.indexOf(session);
      const next = entries[index + 1] ?? entries[index - 1];
      session.view.clear();
      for (const fileView of session.fileViews) fileView.destroy();
      root.remove(); this.entries.delete(id);
      // Deleting output does not silently overwrite an edited form.
      if (this.selected === id) this.selected = null;
      this.announce(`${session.name} deleted. ${this.entries.size} session${this.entries.size === 1 ? '' : 's'} remaining.`);
      (next?.toggle ?? document.getElementById('compute')).focus();
    });
    this.entries.set(id, session); this.nextNumber++;
    this.root.prepend(root); this.select(session, false);
    return summary;
  }
}
