#!/usr/bin/env bash
# Regional overlap guards, including automatic decomposition.
set -euo pipefail
bin="$INTARNABINPATH/src/bin/IntaRNA"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
region_acc=C
common=(--energy=B --seedBP=2 --threads=1 --outMode=C
        --outCsvCols=start1,end1,start2,end2,E --outNumber=10 --default-log-file="$tmp/info.log")

check_regions() {
    local expect=$1
    shift
    : > "$tmp/info.log"
    if "$bin" "${common[@]}" --acc="$region_acc" "$@" > "$tmp/result" 2> "$tmp/error"; then
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

# outMinPu also decomposes regions; blocked positions supply deterministic gaps.
region_acc=N
for overlap in B N T Q; do
    for per_region in false true; do
        for side in target query; do
            expect=ok
            if test "$side" = target; then
                args=(-t CCAACC -q GG --tAccConstr=b:3-4)
                if test "$per_region" = false && [[ "$overlap" = N || "$overlap" = T ]]; then expect=reject; fi
            else
                args=(-t CC -q GGAAGG --qAccConstr=b:3-4)
                if test "$per_region" = false && [[ "$overlap" = N || "$overlap" = Q ]]; then expect=reject; fi
            fi
            check_regions "$expect" "${args[@]}" --outMinPu=0.5 --outOverlap="$overlap" --outPerRegion="$per_region"
        done
    done
done

region_acc=C

# Shifted indices still obey the same manual-region rule.
check_regions reject -t CCAACC -q GG --tIdxPos0=10 --qIdxPos0=20 --tRegion=10-11,14-15 --outOverlap=T
check_regions ok -t CCAACC -q GG --tIdxPos0=10 --qIdxPos0=20 --tRegion=10-11,14-15 --outOverlap=Q
printf 'start1;end1;start2;end2;E\n10;11;20;21;-2\n14;15;20;21;-2\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/result"
echo 'Regional output overlap checks passed'

# Global energy windows must be independent of which region supplies the MFE.
# In this base-pair energy model the optimum is exactly -2; a one-pair site is -1.
for overlap in B T Q; do
    for best_first in true false; do
        if test "$best_first" = true; then seq=CCAAC; regions=1-2,5-5;
        else seq=CAACC; regions=1-1,4-5; fi
        if test "$overlap" = T; then
            seq=${seq//C/G}
            args=(-t CC -q "$seq" --qRegion="$regions")
        else
            args=(-t "$seq" -q GG --tRegion="$regions")
        fi
        for per_region in false true; do
            for delta in 0 0.99 1; do
                "$bin" --energy=B --acc=N --noSeed --model=S --mode=M --threads=1 \
                    --outMode=C --outCsvCols=E --outNumber=10 --outOverlap="$overlap" \
                    --outPerRegion="$per_region" --outDeltaE="$delta" "${args[@]}" \
                    --default-log-file="$tmp/info.log" > "$tmp/energies"
                awk -v delta="$delta" -v per_region="$per_region" '
                    NR == 2 { if ($0 != -2) exit 1 }
                    NR > 1 { if ($0 != -2 && $0 != -1) exit 1; if ($0 == -1) weak++ }
                    END {
                        if (NR < 2) exit 1
                        if (per_region == "false" && delta < 1 && weak) exit 1
                        if ((per_region == "true" || delta == 1) && !weak) exit 1
                    }' "$tmp/energies"
            done
        done
    done
done

# With independent regions all four modes retain the weaker region's optimum.
for overlap in B N T Q; do
    "$bin" --energy=B --acc=N --noSeed --model=S --threads=1 -t CCAAC -q GG \
        --tRegion=1-2,5-5 --outMode=C --outCsvCols=E --outNumber=10 \
        --outOverlap="$overlap" --outPerRegion --outDeltaE=0 \
        --default-log-file="$tmp/info.log" > "$tmp/energies"
    grep -q '^-1$' "$tmp/energies"
done

# An empty merged result has no minimum to dereference.
"$bin" --energy=B --acc=N --noSeed --threads=1 -t AAAAA -q AAAAA --tRegion=1-2,4-5 \
    --outMode=C --outCsvCols=E --outDeltaE=0 --default-log-file="$tmp/info.log" > "$tmp/energies"
printf 'E\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/energies"
echo 'Global and per-region energy window checks passed'
