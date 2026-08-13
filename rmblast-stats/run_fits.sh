#!/bin/sh
# Run ALP Gumbel fits for the RepeatMasker matrix set at canonical gap costs.
# Settings reproduce the original rmblast_*_values fits exactly:
# ALP defaults eps_lambda=0.01 eps_K=0.05, seed=1, iad=false, max_time=10.
cd "$(dirname "$0")" || exit 1
MATDIR=/usr/local/RepeatMasker/Matrices/ncbi/nt
CMPDIR=/usr/local/RepeatModeler/Matrices/ncbi/nt
FIT=./alp-fit/alp_fit
OUT=./fits
mkdir -p "$OUT"

fit_one() {
    name=$1; path=$2; go=$3; ge=$4
    echo "=== $name @ $go/$ge"
    $FIT -matrix "$path" -gapopen "$go" -gapextend "$ge" -max_time 10 \
         -out "$OUT/$name.par" > "$OUT/$name.fit" 2> "$OUT/$name.err" \
        || { echo "FAILED: $name"; return 1; }
    grep '^ncbi_row' "$OUT/$name.fit"
}

for gc in 35 37 39 41 43 45 47 49 51 53; do
    fit_one "14p${gc}g" "$MATDIR/14p${gc}g.matrix" 29 6
    fit_one "18p${gc}g" "$MATDIR/18p${gc}g.matrix" 28 5
    fit_one "20p${gc}g" "$MATDIR/20p${gc}g.matrix" 25 5
    fit_one "25p${gc}g" "$MATDIR/25p${gc}g.matrix" 22 5
done
fit_one "30p53g" "$MATDIR/30p53g.matrix" 22 5
fit_one "comparison" "$CMPDIR/comparison.matrix" 20 5
echo DONE
