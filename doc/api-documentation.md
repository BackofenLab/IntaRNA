# Building and publishing the API reference

The API reference is generated from `src/IntaRNA/*.h` and
`doc/doxygen/mainpage.md`. Its HTML header, footer, and extra stylesheet live in
`doc/doxygen/`; the common configuration is `doc/doxygen.cfg`.

## Local HTML build

Install Doxygen and Graphviz (on Ubuntu: `sudo apt-get install doxygen graphviz`).
The HTML templates are based on Doxygen 1.9.8, the version provided by the
Ubuntu 24.04 CI runner. When updating Doxygen, compare them with freshly generated
templates (`doxygen -w html header.html footer.html default.css`) and verify
navigation and search before publishing.

From the repository root:

```bash
bash doc/build-api.sh
python3 doc/check-api.py doxygen-doc/html
```

Open `doxygen-doc/html/index.html`, or serve the directory with
`python3 -m http.server --directory doxygen-doc/html` and browse
`http://localhost:8000/`. An optional output directory can be passed to the build
script; relative paths are resolved from the calling directory. The script only
builds HTML and does not require a compiler, Boost, ViennaRNA, or an IntaRNA build.
It uses the package version from `configure.ac`.

In a configured Autotools checkout, the existing `make doxygen-doc` target uses
the same content and styling. Configure with `--disable-doxygen-pdf` for an
HTML-only build; PDF generation additionally requires the configured LaTeX tools.

The link check verifies the entry page's local links (including class/member
anchors), its styles and scripts, and the class/namespace/search entry points.
Review Doxygen's warnings as well: existing header documentation can emit
warnings, so the workflow does not treat every legacy documentation warning as
fatal. New unresolved references should be fixed before merging.

## GitHub Pages

`.github/workflows/documentation.yml` builds documentation on pull requests to
`master`, on pushes to `master`, and on manual dispatch. Every successful build
uploads a `github-pages` artifact for review. Only a push or manual dispatch on
`BackofenLab/IntaRNA`'s `master` can deploy; pull requests and forks cannot publish.
The deployment job alone receives Pages and OIDC write permissions.

The workflow first builds the existing README-based Jekyll site with its Slate
theme, then adds Doxygen output under `api/`. It publishes both together:

- User guide: <https://backofenlab.github.io/IntaRNA/>
- Development API: <https://backofenlab.github.io/IntaRNA/api/>

### One-time activation after merge

The repository currently publishes from the root of the `master` branch. After
merging the documentation workflow, a repository administrator must select
**Settings → Pages → Build and deployment → Source → GitHub Actions**, then run
the **Documentation** workflow on `master` (or push a subsequent commit).
Retain any required approval rules on the `github-pages` environment. The PR
does not change the live publishing settings. See GitHub's
[custom Pages workflow documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).

No generated HTML, access token, or separate documentation branch is needed.
Both the existing guide and the API are included in the deployment artifact;
deploying only the API would replace the guide at the site's root.
