#!/usr/bin/env bash
set -euo pipefail
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python3 "$script_dir/testBasePairProbabilities.py" "$INTARNABINPATH/src/bin/IntaRNA"
python3 "$script_dir/testBasePairProbabilitySvg.py" "$INTARNABINPATH/src/bin/IntaRNA"
