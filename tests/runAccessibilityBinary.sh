#!/usr/bin/env bash
# Exercise actual filename dispatch, gzip finalization and prediction reuse.
set -euo pipefail
bin="$INTARNABINPATH/src/bin/IntaRNA"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
query=GGGGAAAACCCCUAGC
target=GCUAGGGGUUUUCCCC
common=(--query="$query" --target="$target" --qIntLenMax=5 --tIntLenMax=5
        --seedBP=3 --outMode=C --outCsvCols=id1,start1,end1,id2,start2,end2,E
        --threads=1 --default-log-file=/dev/null)
"$bin" "${common[@]}" --out="qAcc:$tmp/query.agz" --out="tPu:$tmp/target.AGZ" > "$tmp/computed"
test "$(wc -l < "$tmp/computed")" -gt 1
gzip -t "$tmp/query.agz" "$tmp/target.AGZ"
for mode in E P; do
    "$bin" "${common[@]}" --qAcc="$mode" --tAcc="$mode" \
        --qAccFile="$tmp/query.agz" --tAccFile="$tmp/target.AGZ" > "$tmp/reused"
    cmp "$tmp/computed" "$tmp/reused"
done
# A smaller requested band can reuse a wider stored matrix.
smaller=("${common[@]}")
smaller[2]=--qIntLenMax=3
smaller[3]=--tIntLenMax=3
"$bin" "${smaller[@]}" > "$tmp/smaller"
"$bin" "${smaller[@]}" --qAcc=E --tAcc=E \
    --qAccFile="$tmp/query.agz" --tAccFile="$tmp/target.AGZ" > "$tmp/smaller-reused"
cmp "$tmp/smaller" "$tmp/smaller-reused"
# Preserve the existing compressed text format and dispatch.
"$bin" "${common[@]}" --out="qAcc:$tmp/query.txt.gz" --out="tPu:$tmp/target.txt.gz" > "$tmp/text-computed"
gzip -cd "$tmp/query.txt.gz" > "$tmp/query.txt"
grep -q '^#ensemble' "$tmp/query.txt"
"$bin" "${common[@]}" --qAcc=E --tAcc=P --qAccFile="$tmp/query.txt.gz" \
    --tAccFile="$tmp/target.txt.gz" > "$tmp/text-reused"
cmp "$tmp/text-computed" "$tmp/text-reused"
# Multi-sequence suffixes must keep .agz as the extension.
printf '>q1\n%s\n>q2\n%s\n' "$query" "$target" > "$tmp/queries.fa"
multi=("${common[@]}")
multi[0]="--query=$tmp/queries.fa"
"$bin" "${multi[@]}" --out="qPu:$tmp/multi.agz" > "$tmp/multi-computed"
test -s "$tmp/multi-s1.agz"
test -s "$tmp/multi-s2.agz"
"$bin" "${multi[@]}" --qAcc=P --qAccFile="$tmp/multi.agz" > "$tmp/multi-reused"
cmp "$tmp/multi-computed" "$tmp/multi-reused"
# Damaged trailers must fail even if every ED value has been decompressed.
dd if="$tmp/query.agz" of="$tmp/truncated.agz" bs=1 count="$(( $(wc -c < "$tmp/query.agz") - 4 ))" 2>/dev/null
cp "$tmp/query.agz" "$tmp/checksum.agz"
printf '\000\000\000\000' | dd of="$tmp/checksum.agz" bs=1 \
    seek="$(( $(wc -c < "$tmp/query.agz") - 8 ))" conv=notrunc 2>/dev/null
for damaged in truncated checksum; do
    status=0
    "$bin" "${common[@]}" --qAcc=E --qAccFile="$tmp/$damaged.agz" \
        > "$tmp/bad.out" 2> "$tmp/bad.err" || status=$?
    # A controlled error, not success or termination by a signal.
    test "$status" -eq 255
done
echo 'Binary accessibility CLI checks passed'
