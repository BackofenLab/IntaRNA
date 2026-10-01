#!/usr/bin/env bash
# Exercise evaluation dispatch, validation, energy round trips and personalities.
set -euo pipefail
bin="$INTARNABINPATH/src/bin/IntaRNA"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
common=(--target=GGGGGG --query=UUUUUUU --energy=B --acc=N --threads=1 --default-log-file="$tmp/info.log")
csv=(--outMode=C --outCsvCols=hybridDB,E)
rri='2||.||&2|..|||:1|&1|'
"$bin" "${common[@]}" "${csv[@]}" --rri="$rri" > "$tmp/reference"
printf 'hybridDB;E\n2||.||&2|..|||;-4\n1|&1|;-1\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/reference"

# Distinct inner pairs must survive even with the same energy and boundaries.
"$bin" --target=GGGG --query=CCCC --energy=B --acc=N "${csv[@]}" \
    --rri='1||.|&1||.|:1|.||&1|.||' --default-log-file="$tmp/info.log" > "$tmp/tied"
printf 'hybridDB;E\n1||.|&1||.|;-3\n1|.||&1|.||;-3\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/tied"

# Prediction filters, invalid regions/seeds, length limits and windows are ignored.
"$bin" "${common[@]}" "${csv[@]}" --rri="$rri" --model=B --mode=S --seedBP=20 --seedTQ=invalid \
    --seedQRange=999-1000 --seedNoGU --helixMinBP=4 --helixMaxBP=2 \
    --qRegion=999-1000 --tRegion=999-1000 --qRegionLenMax=1 --tRegionLenMax=1 \
    --intLenMax=1 --qIntLenMax=2 --tIntLenMax=3 --intLoopMax=0 --qIntLoopMax=1 \
    --windowWidth=5 --windowOverlap=100 --outNumber=0 --outOverlap=N --outDeltaE=0 \
    --outMaxE=-999 --outMinPu=1 --outNoLP --outNoGUend --outBestSeedOnly --outPerRegion \
    > "$tmp/ignored" 2> "$tmp/ignored.log"
cmp "$tmp/reference" "$tmp/ignored"
grep -q 'ignoring --outNumber' "$tmp/info.log"
grep -q 'ignoring --seedTQ' "$tmp/info.log"

# Both personality dispatch paths, plus overriding another personality with --rri.
ln -s "$bin" "$tmp/IntaRNAeval"
"$tmp/IntaRNAeval" "${common[@]}" "${csv[@]}" --rri="$rri" > "$tmp/personality"
cmp "$tmp/reference" "$tmp/personality"
for personality in IntaRNAeval IntaRNAseed IntaRNAhelix IntaRNAsTar; do
    "$bin" "${common[@]}" "${csv[@]}" --personality="$personality" --rri="$rri" > "$tmp/personality"
    cmp "$tmp/reference" "$tmp/personality"
done
printf 'rri=%s\nseedTQ=invalid\noutNumber=0\n' "$rri" > "$tmp/parameters"
"$bin" "${common[@]}" "${csv[@]}" --parameterFile="$tmp/parameters" > "$tmp/configured"
cmp "$tmp/reference" "$tmp/configured"

# Energy shifts and positive energies remain visible, even with default outMaxE.
"$bin" "${common[@]}" "${csv[@]}" --rri="$rri" --energyAdd=10 > "$tmp/positive"
printf 'hybridDB;E\n2||.||&2|..|||;6\n1|&1|;9\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/positive"
"$bin" "${common[@]}" "${csv[@]}" --tIdxPos0=-2 --qIdxPos0=8 \
    --rri='-1||.||&9|..|||:-2|&8|' > "$tmp/shifted"
printf 'hybridDB;E\n-1||.||&9|..|||;-4\n-2|&8|;-1\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/shifted"

# Repeated structures do not inflate Zall or the number of results.
for structures in "$rri" "$rri:$rri"; do
    "$bin" "${common[@]}" --outMode=C --rri="$structures" --outCsvCols=hybridDB,E,Zall,Eall,P_E \
        > "$tmp/ensemble-current"
    if test -f "$tmp/ensemble"; then cmp "$tmp/ensemble" "$tmp/ensemble-current"; fi
    cp "$tmp/ensemble-current" "$tmp/ensemble"
done
"$bin" "${common[@]}" --rri="$rri" --outMode=E > "$tmp/ensemble-mode"
test -s "$tmp/ensemble-mode"

# Loops beyond the normal prediction limits can be evaluated.
printf -v dots '%40s' ''
dots=${dots// /.}
printf -v gs '%42s' ''
gs=${gs// /G}
cs=${gs//G/C}
"$bin" --target="$gs" --query="$cs" --acc=N --energy=B --intLoopMax=0 \
    --rri="1|${dots}|&1|${dots}|" --outMode=C --outCsvCols=E --default-log-file="$tmp/info.log" > "$tmp/long-loop"
printf 'E\n-2\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/long-loop"

# Re-evaluate ViennaRNA predictions with the same energy/accessibility settings.
thermo=(--target=AGCGACGCA --query=UGCGUCGCU --accW=0 --accL=0 --temperature=25
        --energyAdd=1.2 --outMode=C --outCsvCols=hybridDB,E,ED1,ED2,E_init,E_loops,E_dangleL,E_dangleR,E_endL,E_endR,E_add
        --threads=1 --default-log-file=/dev/null)
for dangles in '' --energyNoDangles; do
    extra=()
    if test -n "$dangles"; then extra+=("$dangles"); fi
    "$bin" "${thermo[@]}" "${extra[@]}" --noSeed --model=S --mode=M -n 5 > "$tmp/predicted"
    test "$(wc -l < "$tmp/predicted")" -gt 1
    structures=$(awk -F';' 'NR>1 {printf "%s%s", sep, $1; sep=":"}' "$tmp/predicted")
    "$bin" "${thermo[@]}" "${extra[@]}" --rri="$structures" > "$tmp/evaluated"
    cmp "$tmp/predicted" "$tmp/evaluated"
    # All ordinary report formats must accept the evaluated structures.
    for mode in N D; do
        "$bin" --target=AGCGACGCA --query=UGCGUCGCU --rri="$structures" --outMode="$mode" > "$tmp/text"
        test -s "$tmp/text"
    done
done

expect_error() {
    local status=0
    "$bin" "$@" > "$tmp/bad.out" 2> "$tmp/bad.err" || status=$?
    test "$status" -eq 1 || test "$status" -eq 255
}
expect_error "${common[@]}" "${csv[@]}" --personality=IntaRNAeval
for bad in '' '1|&1||' '1...&1...' '0|&1|' '6||&1||' '1|&1|:' '1|&1|::1|&1|' '9223372036854775808|&1|' '1((&1))'; do
    expect_error "${common[@]}" "${csv[@]}" --rri "$bad"
done
expect_error --target=AAAA --query=CCCC --rri='1|&1|'
expect_error --target=GGGG --query=NCCC --rri='1|&1|'
expect_error --target=GGGG --query=CCCC --tIdxPos0=-2 --rri='0|&1|'
printf '>one\nGGGGGG\n>two\nGGGGGG\n' > "$tmp/multiple.fa"
expect_error "${common[@]:1}" --target="$tmp/multiple.fa" --rri="$rri"
expect_error --target=CCCCCC --query="$tmp/multiple.fa" --rri="$rri"
# Retain the selected accessibility model: unavailable intervals fail clearly.
expect_error --target=GGGGGGGG --query=CCCCCCCC --accW=3 --accL=3 --rri='1||||||||&1||||||||'
grep -q 'accessibility' "$tmp/bad.err" "$tmp/bad.out"
echo 'IntaRNAeval CLI checks passed'
