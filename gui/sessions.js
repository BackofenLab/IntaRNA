import { ResultsView, downloadBlob } from './results.js';
import { serializeParameters, importParameters } from './parameter-file.js';
import { zipFiles } from './archive.js';

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
      session.toggle.textContent = session.panel.hidden ? 'Show results' : 'Restore inputs';
    }
  }

  setExpanded(session, expanded) {
    session.panel.hidden = !expanded;
    session.toggle.setAttribute('aria-expanded', String(expanded));
    session.toggle.textContent = expanded ? 'Hide results' : 'Show results';
    // Unload third-party viewers while folded. The table selection is retained
    // and can be visualized again without recomputing the prediction.
    if (!expanded) session.view.visualization.clear();
    else if (session.view.selected) session.view.visualization.show(session.view.selected, session.view.request);
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
    const id = `session-${this.nextNumber}`, task = request.task === 'evaluation' ? 'Evaluation' : 'Prediction';
    const session = {
      id, name: `${task} ${this.nextNumber}`, parameterJson, root,
      panel: find('panel'), toggle: find('toggle'), view: new ResultsView(find('results-view')),
    };
    const count = session.view.show(result.csv, result.csvColumns, request);
    const summary = count ? `${task} complete: ${count} interaction${count === 1 ? '' : 's'}.` : request.task === 'prediction' && request.parameters.outNumber === 0 ? 'Accessibility computation complete. No interaction rows requested.' : 'No favorable interactions found with these sequences and settings.';
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
        const files = [{ name: 'parameter.json', text: parameterJson }, { name: 'intarna-interactions.csv', text: session.view.csv }, ...result.files];
        const filename = session.name.replace(/[^\p{L}\p{N}._-]+/gu, '-').replace(/^[.-]+|[.-]+$/g, '').slice(0, 80) || id;
        downloadBlob(`${filename}.zip`, new Blob([zipFiles(files)], { type: 'application/zip' }));
      } catch (error) { this.reportError(new Error(`ZIP download failed: ${error.message}`)); }
    });
    for (const file of result.files) {
      const details = document.createElement('details'), label = document.createElement('summary'); label.textContent = file.name;
      const preview = document.createElement('pre'); preview.textContent = file.text.slice(0, 100000);
      const download = document.createElement('button'); download.type = 'button'; download.className = 'secondary'; download.textContent = 'Download file';
      download.addEventListener('click', () => downloadBlob(file.name, new Blob([file.text], { type: 'text/plain;charset=utf-8' })));
      details.append(label, download, preview);
      if (file.text.length > 100000) {
        const note = document.createElement('p'); note.textContent = 'Preview limited to 100,000 characters. Download contains the complete file.'; details.append(note);
      }
      find('file-list').append(details);
    }
    find('output-files').hidden = !result.files.length;
    find('diagnostic-text').textContent = result.diagnostics; find('diagnostics').hidden = !result.diagnostics;

    find('delete').addEventListener('click', () => {
      const entries = [...this.entries.values()], index = entries.indexOf(session);
      const next = entries[index + 1] ?? entries[index - 1];
      session.view.clear(); root.remove(); this.entries.delete(id);
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
