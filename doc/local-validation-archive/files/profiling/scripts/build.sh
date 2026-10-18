#!/usr/bin/env bash
set -euo pipefail
profile_root=$(cd "$(dirname "$0")/.." && pwd)
repo_root=${INTARNA_REPO:-"$profile_root/../intarna_optimization_supervised"}
dependency_root=${INTARNA_DEPS:-/home/alex/CursorProjects/IntaRNA_optimization/intaRNA_legacy/.conda-env}
revision=${INTARNA_REVISION:-1b8377989ba79a0f1b8a450163b71bd8ed7c3e45}
if [[ -e "$profile_root/source" ]]; then
    echo 'source already exists; use a new profiling directory to preserve evidence' >&2
    exit 1
fi
mkdir -p "$profile_root/source" "$profile_root/logs"
git -C "$repo_root" archive "$revision" | tar -x -C "$profile_root/source"
cd "$profile_root/source"
autoreconf -fi > ../logs/autoreconf.log 2>&1
./configure --with-boost="$dependency_root" --with-boost-libdir="$dependency_root/lib" \
    --with-vrna="$dependency_root" --prefix="$profile_root/install" \
    CXX="$dependency_root/bin/x86_64-conda-linux-gnu-g++" \
    CC="$dependency_root/bin/x86_64-conda-linux-gnu-gcc" \
    CXXFLAGS='-g -fno-omit-frame-pointer' CPPFLAGS="-I$dependency_root/include" \
    LDFLAGS="-L$dependency_root/lib -Wl,-rpath,$dependency_root/lib" \
    PKG_CONFIG_PATH="$dependency_root/lib/pkgconfig" > ../logs/configure.log 2>&1
make -j4 V=1 > ../logs/build.log 2>&1
make -j4 check > ../logs/check.log 2>&1
