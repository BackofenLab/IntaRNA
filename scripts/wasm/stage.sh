#!/usr/bin/env bash
# Package only the static application, matched WASM assets and distribution notices.
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
build=${WASM_BUILD_DIR:-"$repo/build-wasm"}
site="$build/site"
test -s "$site/wasm/intarna.js"
test -s "$site/wasm/intarna.wasm"
# Select notices from the archive used by these assets, including dlib upgrades.
read -r vienna_version vienna_sha256 dlib_directory < <(python3 - "$site" "$build" <<'PYSOURCE'
import json, re, sys
from pathlib import Path
site, build = map(Path, sys.argv[1:])
metadata = json.loads((site / "build-info.json").read_text())
archive = metadata["viennarnaArchive"]
version, digest = archive["version"], archive["sha256"]
if version != metadata["viennarna"] or not re.fullmatch(r"2\.\d+\.\d+", version) or not re.fullmatch(r"[0-9a-f]{64}", digest):
    raise SystemExit("Invalid compiled ViennaRNA archive identity")
source = build / "sources" / f"vienna-{version}-{digest}" / f"ViennaRNA-{version}"
dlib = list((source / "src").glob("dlib-*/dlib/LICENSE.txt"))
if len(dlib) != 1 or not re.fullmatch(r"dlib-[0-9.]+", dlib[0].parent.parent.name):
    raise SystemExit("Expected one bundled dlib license in the compiled ViennaRNA source")
print(version, digest, dlib[0].parent.parent.name)
PYSOURCE
)
vienna_source="$build/sources/vienna-$vienna_version-$vienna_sha256/ViennaRNA-$vienna_version"
mkdir -p "$site/licenses"
cp -R "$repo/webinterface/assets" "$site/"
for file in index.html styles.css app.js mark.svg protocol.js runner.js results.js options.js columns.js archive.js visualization.js; do
  cp "$repo/webinterface/$file" "$site/$file"
done
cp "$repo/webinterface/intarna-worker.js" "$site/intarna-worker.js"
touch "$site/.nojekyll"
cp "$repo/LICENSE" "$site/licenses/IntaRNA.txt"
cp "$repo/src/easylogging++.LICENCE.txt" "$site/licenses/easylogging.txt"
cp "$repo/src/mdspan/LICENSE" "$site/licenses/mdspan.txt"
cp "$build/sources/boost_1_87_0/LICENSE_1_0.txt" "$site/licenses/Boost.txt"
cp "$build/sources/zlib-1.3.1/zlib.h" "$site/licenses/zlib.h"
cp "$vienna_source/COPYING" "$site/licenses/ViennaRNA.txt"
cp "$vienna_source/src/$dlib_directory/dlib/LICENSE.txt" "$site/licenses/dlib.txt"
cp "$vienna_source/src/ViennaRNA/plotting/naview/naview.c" "$site/licenses/naview.c"
emscripten_dir=$(dirname "$(command -v emcc)")
cp "$emscripten_dir/LICENSE" "$site/licenses/Emscripten.txt"
cp "$emscripten_dir/system/lib/libc/musl/COPYRIGHT" "$site/licenses/musl.txt"
for library in libcxx libcxxabi libunwind compiler-rt; do
  cp "$emscripten_dir/system/lib/$library/LICENSE.TXT" "$site/licenses/$library.txt"
done
cat > "$site/licenses/index.html" <<HTML
<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>IntaRNA licenses</title>
<h1>IntaRNA and third-party notices</h1>
<p>IntaRNA: Bioinformatics Group of Prof. Backofen, University of Freiburg.<br>ViennaRNA: authors and Institute for Theoretical Chemistry, University of Vienna.</p>
<p>This distribution uses ViennaRNA $vienna_version, Boost 1.87.0, zlib 1.3.1, bundled dlib ${dlib_directory#dlib-}, Kokkos mdspan, Easylogging++ and the Emscripten 4.0.23 runtime.</p>
<ul>
<li><a href="IntaRNA.txt">IntaRNA</a></li><li><a href="ViennaRNA.txt">ViennaRNA</a></li>
<li><a href="Boost.txt">Boost</a></li><li><a href="zlib.h">zlib</a></li>
<li><a href="dlib.txt">dlib</a></li><li><a href="naview.c">Naview</a></li>
<li><a href="mdspan.txt">Kokkos mdspan</a></li><li><a href="easylogging.txt">Easylogging++</a></li>
<li><a href="Emscripten.txt">Emscripten</a></li><li><a href="musl.txt">musl</a></li>
<li><a href="libcxx.txt">libc++</a></li><li><a href="libcxxabi.txt">libc++abi</a></li>
<li><a href="libunwind.txt">libunwind</a></li><li><a href="compiler-rt.txt">compiler-rt</a></li>
</ul><p><a href="../build-info.json">Build provenance</a> · <a href="../">Back to IntaRNA</a></p>
</html>
HTML
# Provenance belongs to the compiled assets, never the later staging revision.
python3 - "$site" <<'PYINFO'
import hashlib, json, sys
from pathlib import Path
site = Path(sys.argv[1])
metadata = json.loads((site / "build-info.json").read_text())
if not (site / "build-info.txt").is_file():
    raise SystemExit("Missing compiled version text")
for name, digest in metadata["sha256"].items():
    if hashlib.sha256((site / "wasm" / name).read_bytes()).hexdigest() != digest:
        raise SystemExit(f"Mismatched engine: {name}")
PYINFO
echo "Static site: $site"
