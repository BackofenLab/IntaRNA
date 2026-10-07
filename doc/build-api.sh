#!/usr/bin/env bash
# Build the HTML API without configuring or compiling the C++ library.
set -euo pipefail

if (( $# > 1 )); then
  echo "Usage: $0 [output-directory]" >&2
  exit 2
fi

source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
output_dir="${1:-$source_dir/doxygen-doc}"
mkdir -p -- "$output_dir"
output_dir="$(cd -- "$output_dir" && pwd)"
cd -- "$source_dir"

version="$(sed -n 's/^AC_INIT(\[IntaRNA\], \[\([^]]*\)\].*/\1/p' configure.ac)"
test -n "$version"
for tool in doxygen dot; do
  if ! command -v "$tool" >/dev/null; then
    echo "Missing documentation tool: $tool (install Doxygen and Graphviz)" >&2
    exit 1
  fi
done
dot_path="$(dirname -- "$(command -v dot)")"

# These are the same substitutions provided by m4/doxygen.m4 for make doxygen-doc.
env SRCDIR=. PROJECT=IntaRNA VERSION="$version" DOCDIR="$output_dir" \
  GENERATE_HTML=YES GENERATE_LATEX=NO GENERATE_PDF=NO \
  GENERATE_HTMLHELP=NO GENERATE_CHI=NO GENERATE_RTF=NO \
  GENERATE_MAN=NO GENERATE_XML=NO PAPER_SIZE=a4 \
  HAVE_DOT=YES DOT_PATH="$dot_path" \
  doxygen doc/doxygen.cfg

test -s "$output_dir/html/index.html"
echo "API documentation: $output_dir/html/index.html"
