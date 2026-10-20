# Step 4: results and Pages deployment (#286)

## Results contract

The worker requests exactly `id1,start1,end1,id2,start2,end2,E`. IntaRNA's native
CSV is unquoted and semicolon-delimited; its output handler replaces semicolons
inside sequence identifiers with underscores. The browser parser treats quotes
as literal ID characters, checks the header, validates finite numeric values
and positive inclusive coordinates, and rejects malformed rows. Header control
characters/multirecord input are rejected before the CLI runs. Empty CSV or a
header-only result is a successful prediction with no favorable interactions.

Sequence 1 is **target** and sequence 2 is **query**. The table deliberately
shows query first, with each strand's original 1-based inclusive coordinates.
Energy uses kcal/mol. ΔE is derived as `E - min(reported E)` rounded to hundredths;
it is not a probability. Suboptimal interactions have distinct site boundaries
and may overlap. Scientific prediction defaults are unchanged.

Every heading is a keyboard-operable sort button with `aria-sort`. Numerical
columns compare numbers, IDs compare text, and the search matches all columns
without changing the raw prediction. All values are rendered through
`textContent`. CSV export returns the original native output, including all
reported rows regardless of filtering. Input edits clear the old result.

## Packaging and deployment

```sh
bash scripts/wasm/build.sh
node scripts/wasm/smoke.mjs
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
files needed at runtime, source commit/version metadata and dependency notices.
It needs the same activated SDK and dependency build directory as the build.
No compiled artifacts are committed. HTTP(S), module workers and WebAssembly
are required; `file://` and legacy browsers are unsupported. No backend,
SharedArrayBuffer, credentials, CDN or cross-origin headers are needed.

The existing repository already publishes its user guide and API with Pages.
Publishing a second independent artifact would overwrite them. Therefore
`documentation.yml` is now a reusable build (also manually dispatchable), and
`build-wasm.yml` is the **single publisher**. Both builds run on PRs and on the
actual default branch, `master`. After JS tests, the real WASM smoke test and
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
