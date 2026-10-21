# Step 4: results and Pages deployment (#286)

## Results contract

The CSV checklist covers all `OutputHandlerCsv` column IDs. The worker requests
those columns plus internal target/query IDs, coordinates and energy needed for
pair grouping and result validation. Downloads project the untouched native CSV
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

## Packaging and deployment

```sh
bash scripts/wasm/build.sh
node scripts/wasm/smoke.mjs
node scripts/wasm/smoke-protocol.mjs
bash scripts/wasm/stage.sh
npm test --prefix webinterface
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
- Tested application at `/IntaRNA/webinterface/`.

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
- `npm test --prefix webinterface`: all five input/lifecycle/CSV test groups pass.
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
- `npm test --prefix webinterface`: ten protocol/results/lifecycle groups pass.
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
