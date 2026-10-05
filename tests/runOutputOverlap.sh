#!/usr/bin/env bash
# Regional overlap guards, including automatic decomposition.
set -euo pipefail
bin="$INTARNABINPATH/src/bin/IntaRNA"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
common=(--energy=B --acc=C --seedBP=2 --threads=1 --outMode=C
        --outCsvCols=start1,end1,start2,end2,E --outNumber=10 --default-log-file="$tmp/info.log")

check_regions() {
    local expect=$1
    shift
    : > "$tmp/info.log"
    if "$bin" "${common[@]}" "$@" > "$tmp/result" 2> "$tmp/error"; then
        if test "$expect" != ok; then
            echo "Expected regional-overlap rejection: $*" >&2
            exit 1
        fi
    else
        if test "$expect" != reject; then cat "$tmp/error" "$tmp/info.log" "$tmp/result" >&2; exit 1; fi
        grep -q 'regions cannot be merged with --outOverlap=' "$tmp/error" "$tmp/info.log" "$tmp/result"
    fi
}

for overlap in B N T Q; do
    for per_region in false true; do
        for target_regions in 1 2; do
            for query_regions in 1 2; do
                expect=ok
                if test "$per_region" = false; then
                    if test "$target_regions" = 2 && [[ "$overlap" = N || "$overlap" = T ]]; then expect=reject; fi
                    if test "$query_regions" = 2 && [[ "$overlap" = N || "$overlap" = Q ]]; then expect=reject; fi
                fi
                tr=1-2; qr=1-2
                if test "$target_regions" = 2; then tr+=,5-6; fi
                if test "$query_regions" = 2; then qr+=,5-6; fi
                check_regions "$expect" -t CCAACC -q GGAAGG --tRegion="$tr" --qRegion="$qr" \
                    --outOverlap="$overlap" --outPerRegion="$per_region"
            done
        done
        # Decomposition really produces multiple regions on the selected RNA.
        for side in target query; do
            expect=ok
            if test "$side" = target; then
                args=(-t UAUCGGCC -q GG --tRegionLenMax=4)
                if test "$per_region" = false && [[ "$overlap" = N || "$overlap" = T ]]; then expect=reject; fi
            else
                args=(-t GG -q UAUCGGCC --qRegionLenMax=4)
                if test "$per_region" = false && [[ "$overlap" = N || "$overlap" = Q ]]; then expect=reject; fi
            fi
            check_regions "$expect" "${args[@]}" --outOverlap="$overlap" --outPerRegion="$per_region"
        done
    done
    # An automatic-region option alone is fine if no split is needed.
    check_regions ok -t CC -q GG --tRegionLenMax=4 --qRegionLenMax=4 --outOverlap="$overlap"
done

# Shifted indices still obey the same manual-region rule.
check_regions reject -t CCAACC -q GG --tIdxPos0=10 --qIdxPos0=20 --tRegion=10-11,14-15 --outOverlap=T
check_regions ok -t CCAACC -q GG --tIdxPos0=10 --qIdxPos0=20 --tRegion=10-11,14-15 --outOverlap=Q
printf 'start1;end1;start2;end2;E\n10;11;20;21;-2\n14;15;20;21;-2\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/result"
echo 'Regional output overlap checks passed'
