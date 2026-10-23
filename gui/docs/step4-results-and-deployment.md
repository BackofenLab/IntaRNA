# Step 4: results and Pages deployment (#286)

## Results contract

The CSV checklist covers all `OutputHandlerCsv` column IDs. The worker requests
those columns plus internal target/query IDs, coordinates, energy and `hybridDP`
needed for pair grouping, result validation and visualization. Downloads project the untouched native CSV
cells onto the user's selected columns, in selection order, retaining every row
regardless of table filtering/visibility. No computed mfe difference is inserted
into native CSV. The selected-column catalog is checked against the C++ mapping.

Native output is unquoted and semicolon-delimited. Semicolons in IDs become
underscores; quotes remain literal. The parser checks column IDs, finite core
energy/coordinate values, integer inclusive coordinates and ascending bounds.
Signed coordinates are allowed. Optional `NAN` and colon-separated seed lists
remain strings; numeric sorting puts unavailable/list values after numbers.

Sequence 1 is **target**, sequence 2 is **query**. Energy uses kcal/mol. The
computed **mfe Δ** is `E - min(reported E for the same query–target pair)`, rounded
to hundredths. Evaluation uses the best supplied structure as its reference.
This is relative suboptimality, not a free-energy change or probability.

Every heading is a keyboard-operable sort button with `aria-sort`. Search covers
all returned values. A separate checklist shows/hides table columns after the
run; it does not change selected CSV output. All identifiers, output previews
and diagnostics render with `textContent`. Selected extra files appear below the
table with bounded text previews and complete downloads, including SVG source.
Empty output is successful and still permits CSV/accessibility downloads.

## Selected interaction in vaRRI

Each row has a keyboard-accessible View button; clicking a row also selects it.
Selection survives table sorting, filtering and visibility changes. Input edits,
new runs, errors and cancellation clear the viewer and archive snapshot.
`visualization.js` uses native seq1 = target, seq2 = query and maps the native
CSV identifiers to the completed run's normalized FASTA records (including
semicolon-to-underscore ID conversion). It pads the interaction dot-bracket
structure with unpaired positions to the full input lengths, using each start
index and the native negative-to-positive transition that skips zero.

The URL follows [vaRRI's embedding API](https://github.com/BackofenLab/vaRRI/blob/main-v2/docs/sharing-and-input.md):
`sequence`, `structure`, `seqName1/2`, `startIndex1/2`, and
`showRenderingOnly=true`. Names use the main-v2 API; older deployed vaRRI versions
may retain generic strand names. A separate Open in vaRRI link opens the full
viewer. The local status always identifies the selected pair and coordinates.
Zero-based CLI start indices remain accepted for computation/export, but vaRRI
requires nonzero start indices; choosing View for such a run explains this limit.

No external viewer loads until a row is selected. The interface explains that
selection shares that pair and structure through vaRRI's URL. Computation and
ZIP generation remain local. CSP permits frames only at the vaRRI deployment;
the iframe sends no referrer. Viewer availability depends on that external site;
the full-viewer link permits retry, and scientific outputs remain downloadable.

## Reproducible ZIP download

Download ZIP sits beside Download CSV. `archive.js` creates a standard UTF-8,
uncompressed ZIP with CRC-32 checksums locally, without a library or CDN. It uses
the completed request and result, independent of table sorting/filtering and
column visibility. Editing inputs discards the completed snapshot.

The archive contains:

- `query.fa` and `target.fa`: the exact normalized FASTA consumed by the worker.
- `parameters.cfg`: native parameter-file syntax, non-default effective analysis
  settings, explicit input/output routing and the selected CSV columns. Comments
  include the exact compiled version/provenance text and the replay command
  `IntaRNA --parameterFile=parameters.cfg` from the extracted directory.
- `build-info.txt`: original compiled engine/library versions and asset hashes.
- `intarna-interactions.csv`: all results with the selected CSV columns, identical
  to Download CSV even if the table has been sorted, filtered or columns hidden.
- Every generated additional output, with its original filename and contents.
- Local SHAPE/accessibility input files when supplied, referenced by relative paths.
- Separate `normal.cfg`, `detailed.cfg` or `ensemble.cfg` when those reports were
  requested, with their replay commands also listed in `parameters.cfg`. Each
  report requires a separate native invocation, matching worker behavior.

Replay with the recorded IntaRNA and library versions. Export translates MEMFS
paths to portable relative paths and omits defaults using the worker's own
argument builder. ZIP download reports a clear error if compiled version text
is unavailable, rather than claiming incomplete provenance.

## Packaging and deployment

```sh
bash scripts/wasm/build.sh
node scripts/wasm/smoke.mjs
node scripts/wasm/smoke-protocol.mjs
bash scripts/wasm/stage.sh
npm test --prefix gui
python3 -m pip install -r scripts/wasm/requirements-browser.txt
python3 -m playwright install chromium firefox
python3 -m http.server 8765 --directory build-wasm/site
# In another terminal:
python3 scripts/wasm/test-browser.py http://127.0.0.1:8765
python3 scripts/wasm/test-browser.py http://127.0.0.1:8765 --browser firefox
```

The staging script includes matched WASM assets, only the application source
files needed at runtime, compiled version metadata and dependency notices. Build-time `build-info.json`
and `build-info.txt` record actual CLI IntaRNA/ViennaRNA/Boost versions, source
commit and tracked modifications, SDK version and asset SHA-256 digests. Staging
checks these digests and preserves the original metadata; it never replaces the
compiled revision with a later staging checkout. The footer displays versions.
It needs the same activated SDK and dependency build directory as the build.
No compiled artifacts are committed. HTTP(S), module workers and WebAssembly
are required; `file://` and legacy browsers are unsupported. No backend,
SharedArrayBuffer, credentials, CDN or cross-origin headers are needed.

The existing repository already publishes its user guide and API with Pages.
Publishing a second independent artifact would overwrite them. Therefore
`documentation.yml` is now a reusable build (also manually dispatchable), and
`build-wasm.yml` is the **single publisher**. Both builds run on PRs (including dependent feature-branch PRs) and on
the actual default branch, `master`. After JS tests, the real WASM smoke test and
Chromium/Firefox integration tests pass, the package job combines:

- Existing Jekyll guide at `/IntaRNA/` and API at `/IntaRNA/api/`.
- Tested application at `/IntaRNA/gui/`.

The artifact `intarna-pages-preview` permits PR review of the complete site.
Pull requests and manual runs on feature branches do not publish. Only a push
or manual run on `BackofenLab/IntaRNA`'s `master` may upload the Pages artifact
and use the `github-pages` deployment environment. Pages permissions are scoped
to that deploy job. Concurrent deployments are serialized without cancellation.

The repository's current Pages source is already **GitHub Actions**. A fork
would need that source enabled, its repository guard adjusted deliberately, and
an environment allowing its default branch. Existing environment protection
rules still apply. After merging this PR, check the workflow's deploy job and
its returned environment URL. This implementation PR does not merge itself or
claim that deployment happened before merge.

## Validation scope

Node tests cover input boundaries, lifecycle races, parser failures, numeric
sorting and suboptimal differences. The WASM smoke test compares a real
Turner/default-accessibility prediction with the native reference. Browser
integration runs real WASM (no mocked scientific results), checking coordinates,
energy, sorting, filtering, CSV download, mobile overflow, abort during download
and active C++, clean restart, empty predictions and asset-load recovery.

### Local verification, 2026-10-20

- Emscripten 4.0.23 with the exact SHA-256-pinned official archives: successful
  dependency probe, cross-build and real WASM smoke test. The structured-strand
  reference reports E = -28.04, ED1 = 4.60, ED2 = 4.86 kcal/mol; the short example
  reports E = -21.40 with target/query coordinates 1–10 and ten suboptimals.
  Invalid CLI input returns a nonzero status through C++ error handling.
- `npm test --prefix gui`: all five input/lifecycle/CSV test groups pass.
- Playwright 1.64.0, Chromium and Firefox: full integration suite passes with the
  official-dependency engine, including a nested `/site/` URL as used by Pages,
  active-C++ cancellation, literal HTML-like identifiers and retry after failure.
  Desktop screenshot visually reviewed; 390 px viewport has no page overflow.
- Native release build with OpenMP, conda GCC 16.2.0 and ViennaRNA 2.7.2:
  `make tests -j2` passes all eight API/CLI groups, including 462,689 assertions
  across 119 API cases. For the separate build directory, tracked `tests/data`
  fixtures were linked into the build's data directory because the existing
  shell harness discovers fixtures relative to its working directory.
- Shell syntax, workflow YAML parsing and `git diff --check`: pass.

GitHub's platform checks independently cover GCC 14 release/debug and Apple
Clang. Deployment itself occurs only after merge; local browser tests do not
claim a live Pages deployment.

### PR #287 review verification, 2026-10-21

- Emscripten 4.0.23 rebuild using the pinned official dependency archives passes
  with `-Werror=inconsistent-missing-override`. An all-public-header audit using
  `-Wsuggest-override -Wsuggest-destructor-override` reports no project override
  diagnostics (third-party Easylogging++ is intentionally unchanged).
- `node scripts/wasm/smoke.mjs` retains the original native-reference energy,
  coordinates, structured accessibility and C++ error-handling checks.
- `node scripts/wasm/smoke-protocol.mjs` passes real CLI translations for default
  prediction, shifted negative coordinates, all-vs-all/pairwise multi-FASTA,
  evaluation with all CSV columns, extra reports/probability SVGs, imported
  accessibility and non-default energy/seed/helix/SHAPE/overlap settings.
  Native imported accessibility reserves one length column for dangling-end
  treatment; its energy is compared with the corresponding restricted native
  interaction length rather than the unrestricted prediction.
- `npm test --prefix gui`: ten protocol/results/lifecycle groups pass.
- Real Chromium/Firefox integration covers all parameter groups, CLI defaults,
  compiled versions/logos, task changes, multiple records, signed coordinates,
  selected/visible CSV columns, additional previews/downloads and the original
  cancellation, empty/error, literal-ID and mobile checks. Desktop screenshot
  reviewed; the 390 px viewport has no document overflow.
- Native OpenMP release build with GCC 16.2.0/ViennaRNA 2.7.2:
  `make tests -j3` passes all eight groups, including 462,689 assertions across
  119 API cases. Tracked CLI fixtures are linked into the separate build data
  directory as in the previous verification; expected results are unchanged.
- Staging preserves compiled metadata byte-for-byte; a deliberately changed
  copied JS engine is rejected by its asset digest check. Shell syntax,
  workflow YAML, Python syntax and `git diff --check` pass.

### PR #287 follow-up verification, 2026-10-22

- All 12 Node test groups pass, including coordinate/strand mapping, non-default
  export settings, auxiliary inputs and report configurations.
- `node scripts/wasm/smoke.mjs` and `node scripts/wasm/smoke-protocol.mjs` pass
  using the previously compiled Emscripten 4.0.23 engine (revision `fc62a23`).
  The latter executes exported parameter files in fresh engines with only
  archived inputs and verifies byte-identical CSV, per-record outputs,
  probability SVGs and ASCII/detailed/ensemble reports. Cases include signed
  multi-FASTA, changed energy/seed settings, SHAPE and accessibility imports.
- Chromium 156 and Firefox 157 with Playwright 1.64.0 pass the real-WASM browser
  suite at a nested URL. It covers tab keyboard behavior, selection/encoded
  viewer URLs, Python `zipfile` validation of downloaded ZIP contents, hidden
  CSV columns, evaluation reports, empty results and the existing lifecycle
  checks. Only the third-party iframe response is stubbed in CI; computation
  runs the real engine. Desktop/parameter/mobile screenshots were reviewed.
- A separate live Chromium check renders RNA base circles in the actual hosted
  vaRRI iframe. An independent Python `zipfile` check verifies CRCs, Unicode
  filenames/data and empty files in the ZIP writer's output.
- Native `make tests -j2` passes all eight API/CLI groups (462,689 assertions in
  119 API cases). No C++ or compiled engine changes were needed for this follow-up.
- Staging, shell/Python syntax and `git diff --check` pass.
