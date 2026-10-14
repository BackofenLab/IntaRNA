#!/usr/bin/env bash
# Compare complete accessibility output and predictions across worker counts.
set -euo pipefail
bin="$INTARNABINPATH/src/bin/IntaRNA"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
export OMP_NUM_THREADS=2 OMP_THREAD_LIMIT=2
"$bin" --fullhelp > "$tmp/help"
if ! grep -q -- '--threads' "$tmp/help"; then
    exit 77
fi
cat > "$tmp/sequences.fa" <<'EOF'
>short
GGGGAAAACCCCUAGC
>medium
GCUAGGGGUUUUCCCCGCUA
>long
GGGGAAAACCCCUAGCGGGGAAAA
EOF
for i in $(seq 1 16); do
    printf '%s 0.5\n' "$i"
done > "$tmp/reactivities.shape"
common=(--query="$tmp/sequences.fa" --target="$tmp/sequences.fa"
        --qIntLenMax=8 --tIntLenMax=8 --seedBP=3 --qAccW=14 --tAccW=18
        --qAccL=10 --tAccL=12
        --outMode=C
        --outMaxE=100 --default-log-file=/dev/null)
# D converts reactivities; W reads precomputed pseudoenergies. See the audit
# for pre-existing serial failures of Z and W combined with global Eall output.
for method in none D W; do
    shape=()
    columns=id1,start1,end1,id2,start2,end2,E,ED1,ED2
    if test "$method" != W; then
        columns+=,Eall1,Eall2
    fi
    if test "$method" != none; then
        shape=(--qShape="$tmp/reactivities.shape" --tShape="$tmp/reactivities.shape"
               --qShapeMethod="$method" --tShapeMethod="$method")
    fi
    for threads in 1 2; do
        mkdir "$tmp/$method-$threads"
        "$bin" "${common[@]}" "${shape[@]}" --threads="$threads" --outCsvCols="$columns" \
            --out="qAcc:$tmp/$method-$threads/query.txt" \
            --out="tAcc:$tmp/$method-$threads/target.txt" \
            > "$tmp/$method-$threads/predictions.unsorted" || {
                cat "$tmp/$method-$threads/predictions.unsorted"
                exit 1
            }
        # Workers may publish sequence pairs in a different order.
        LC_ALL=C sort "$tmp/$method-$threads/predictions.unsorted" \
            > "$tmp/$method-$threads/predictions"
        test "$(wc -l < "$tmp/$method-$threads/predictions")" -gt 1
    done
    cmp "$tmp/$method-1/predictions" "$tmp/$method-2/predictions"
    for strand in query target; do
        for sequence in 1 2 3; do
            cmp "$tmp/$method-1/$strand-s$sequence.txt" \
                "$tmp/$method-2/$strand-s$sequence.txt"
        done
    done
done
echo 'Parallel accessibility, SHAPE and prediction checks passed'
