#!/usr/bin/env bash
# Exercise kinetic dispatch, greedy paths, output constraints and energy accounting.
set -euo pipefail
bin="$INTARNABINPATH/src/bin/IntaRNA"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
common=(--target=GGGGGG --query=CCCCCC --energy=B --acc=N
        --threads=1 --default-log-file=/dev/null)
kinetic=(--model=X --mode=K '--seedTQ=3||&3||')
csv=(--outMode=C --outCsvCols=hybridDB,E)

# Six base pairs have independently known energy -6 in the base-pair model.
printf 'hybridDB;E\n1||||||&1||||||;-6\n' > "$tmp/expected"
"$bin" "${common[@]}" "${kinetic[@]}" "${csv[@]}" > "$tmp/default"
cmp "$tmp/expected" "$tmp/default"
for score in A B C; do
    for noLP in false true; do
        "$bin" "${common[@]}" "${kinetic[@]}" "${csv[@]}" \
            --kineticScore="$score" --outNoLP="$noLP" --outNoGUend > "$tmp/score"
        cmp "$tmp/expected" "$tmp/score"
    done
done

# Ties choose the left extension, including the antiparallel query coordinates.
"$bin" "${common[@]}" "${kinetic[@]}" "${csv[@]}" --intLenMax=3 > "$tmp/left"
printf 'hybridDB;E\n2|||&3|||;-3\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/left"

# Boundary-only output and CLI range/index conversions use the same trajectory.
"$bin" "${common[@]}" "${kinetic[@]}" --outMode=C \
    --outCsvCols=start1,end1,start2,end2,E > "$tmp/boundaries"
printf 'start1;end1;start2;end2;E\n1;6;1;6;-6\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/boundaries"
"$bin" "${common[@]}" "${kinetic[@]}" "${csv[@]}" \
    --tRegion=2-5 --qRegion=2-5 > "$tmp/ranges"
printf 'hybridDB;E\n2||||&2||||;-4\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/ranges"
# Signed IntaRNA coordinates skip zero: target position three is labeled 1.
"$bin" "${common[@]}" --model=X --mode=K '--seedTQ=1||&7||' "${csv[@]}" \
    --tIdxPos0=-2 --qIdxPos0=5 > "$tmp/shifted"
printf 'hybridDB;E\n-2||||||&5||||||;-6\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/shifted"

# Every committed prefix is available for suboptimal output, in energy order.
"$bin" "${common[@]}" "${kinetic[@]}" "${csv[@]}" --outNumber=10 > "$tmp/prefixes"
printf 'hybridDB;E\n1||||||&1||||||;-6\n1||||&3||||;-4\n3||&3||;-2\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/prefixes"

# A zero report count still permits energy tracking and never needs traceback.
"$bin" "${common[@]}" "${kinetic[@]}" "${csv[@]}" --outNumber=0 \
    --out="tMinE:$tmp/min-energy" > "$tmp/zero"
printf 'hybridDB;E\n' > "$tmp/empty"
cmp "$tmp/empty" "$tmp/zero"
test -s "$tmp/min-energy"
grep -q -- '-6' "$tmp/min-energy"

# Handler-provided explicit seeds can contain lonely pairs.
"$bin" "${common[@]}" --model=X --mode=K '--seedTQ=3|&3|' \
    "${csv[@]}" --outNoLP > "$tmp/lonely"
printf 'hybridDB;E\n1|||||&1|||||;-5\n' > "$tmp/expected"
cmp "$tmp/expected" "$tmp/lonely"
"$bin" --target=GG --query=UU --energy=B --acc=N --model=X --mode=K \
    '--seedTQ=1||&1||' "${csv[@]}" --outNoGUend --default-log-file=/dev/null > "$tmp/gu"
cmp "$tmp/empty" "$tmp/gu"

# Config-file selection must behave just like command-line selection.
printf 'model=X\nmode=K\nkineticScore=C\nseedTQ=3||&3||\n' > "$tmp/parameters"
"$bin" "${common[@]}" "${csv[@]}" --parameterFile="$tmp/parameters" > "$tmp/configured"
cmp "$tmp/default" "$tmp/configured"

# Re-evaluate all reported structures independently through IntaRNAeval.
# This checks the complete energy, including accessibility, terminal and dangle terms.
thermo=(--target=AGCGACGCA --query=UGCGUCGCU --accW=0 --accL=0 --temperature=25
        --energyAdd=1.2 --outMode=C --outMaxE=100 --outNumber=5
        --outCsvCols=hybridDB,E,ED1,ED2,E_init,E_loops,E_dangleL,E_dangleR,E_endL,E_endR,E_add
        --threads=1 --default-log-file=/dev/null)
for score in A B C; do
    for dangles in false true; do
        "$bin" "${thermo[@]}" --model=X --mode=K '--seedTQ=3|||&5|||' \
            --kineticScore="$score" --energyNoDangles="$dangles" --outNoLP --outNoGUend \
            > "$tmp/predicted"
        test "$(wc -l < "$tmp/predicted")" -gt 1
        structures=$(awk -F';' 'NR>1 {printf "%s%s", sep, $1; sep=":"}' "$tmp/predicted")
        "$bin" "${thermo[@]}" --rri="$structures" --energyNoDangles="$dangles" > "$tmp/evaluated"
        cmp "$tmp/predicted" "$tmp/evaluated"
    done
done

expect_error() {
    local expected="$1"
    shift
    local status=0
    "$bin" "$@" > "$tmp/bad.out" 2> "$tmp/bad.err" || status=$?
    test "$status" -eq 1 || test "$status" -eq 255
    grep -q -- "$expected" "$tmp/bad.out" "$tmp/bad.err"
}
for model in S P B; do
    expect_error 'only with --model=X' "${common[@]}" --model="$model" --mode=K
done
expect_error 'incompatible with --noSeed' "${common[@]}" "${kinetic[@]}" --noSeed
# Explicitly supplying the default A must be rejected outside K as well.
for score in A B C; do
    expect_error 'kineticScore requires' "${common[@]}" --mode=H --kineticScore="$score"
done
expect_error 'kineticScore' "${common[@]}" "${kinetic[@]}" --kineticScore=D
expect_error 'equilibrium ensemble' "${common[@]}" "${kinetic[@]}" --outMode=E
for column in Zall Eall EallTotal P_E; do
    expect_error 'equilibrium ensemble' "${common[@]}" "${kinetic[@]}" \
        --outMode=C --outCsvCols="E,$column"
done
for output in spotProb qSpotProb tSpotProb; do
    expect_error 'equilibrium ensemble' "${common[@]}" "${kinetic[@]}" --out="$output:$tmp/probability"
done
expect_error 'equilibrium ensemble' "${common[@]}" "${kinetic[@]}" --out="spotProb:1&1:$tmp/probability"

# Evaluation ignores prediction controls, including invalid kinetic scores.
"$bin" "${common[@]}" "${csv[@]}" '--rri=1||||||&1||||||' \
    --model=S --mode=K --kineticScore=D > "$tmp/evaluation"
cmp "$tmp/default" "$tmp/evaluation"
# A missing or false flag is promoted once, with a visible INFO message.
logging=("${common[@]}")
logging[${#logging[@]}-1]="--default-log-file=$tmp/info.log"
for setting in absent false true; do
    : > "$tmp/info.log"
    flags=()
    if [ "$setting" != absent ]; then flags=(--outNoLP="$setting"); fi
    "$bin" "${logging[@]}" "${kinetic[@]}" "${csv[@]}" "${flags[@]}" > "$tmp/info.out" 2> "$tmp/info.err"
    if [ "$setting" = true ]; then
        ! grep -q 'setting --outNoLP=true' "$tmp/info.log"
    else
        grep -q 'INFO.*setting --outNoLP=true' "$tmp/info.log"
    fi
done
# Both personality entry points select mode K and noLP by default.
ln -s "$bin" "$tmp/IntaRNAkix"
for invocation in binary option; do
    args=()
    executable="$tmp/IntaRNAkix"
    if [ "$invocation" = option ]; then
        executable="$bin"
        args=(--personality=IntaRNAkix)
    fi
    : > "$tmp/info.log"
    "$executable" "${args[@]}" "${logging[@]}" '--seedTQ=3||&3||' "${csv[@]}" \
        --outNumber=10 > "$tmp/kix-prefixes"
    cmp "$tmp/prefixes" "$tmp/kix-prefixes"
    ! grep -q 'setting --outNoLP=true' "$tmp/info.log"
    "$executable" "${args[@]}" "${common[@]}" '--seedTQ=3||&3||' "${csv[@]}" \
        --mode=S > "$tmp/seed-only"
    printf 'hybridDB;E\n3||&3||;-2\n' > "$tmp/expected"
    cmp "$tmp/expected" "$tmp/seed-only"
done
# Explicitly disabling noLP cannot disable the mode K extension invariant.
"$bin" --personality=IntaRNAkix "${logging[@]}" '--seedTQ=3||&3||' "${csv[@]}" \
    --outNoLP=false > "$tmp/kix-noLP"
cmp "$tmp/default" "$tmp/kix-noLP"
grep -q 'setting --outNoLP=true' "$tmp/info.log"
# Evaluation remains available under the personality and ignores its defaults.
"$tmp/IntaRNAkix" "${common[@]}" "${csv[@]}" '--rri=1||||||&1||||||' > "$tmp/kix-eval"
cmp "$tmp/default" "$tmp/kix-eval"
echo 'Kinetic seed-extension and IntaRNAkix CLI checks passed'
