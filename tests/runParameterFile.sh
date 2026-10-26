#!/bin/sh
set -eu
command -v python3 >/dev/null 2>&1 || exit 77
exec python3 "${srcdir:-.}/testParameterFile.py" "$INTARNABINPATH/src/bin/IntaRNA"
