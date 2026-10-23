#!/usr/bin/env bash
# Cache only verified dependency archives; recover from invalid cached responses.
set -euo pipefail
url=$1
destination=$2
checksum=$3

# Verify existing archives too, including caches from earlier build scripts.
if [[ -f "$destination" ]]; then
  if printf '%s  %s\n' "$checksum" "$destination" | sha256sum --check --status; then
    echo "$destination: OK (cached)"
    exit 0
  fi
  echo "Discarding invalid cached archive: $destination" >&2
  rm -- "$destination"
fi

# Never promote HTTP error pages or incomplete transfers to the download cache.
partial="$destination.part"
trap 'rm -f -- "$partial"' EXIT
curl --fail --location --retry 3 "$url" -o "$partial"
printf '%s  %s\n' "$checksum" "$partial" | sha256sum --check -
mv -- "$partial" "$destination"
