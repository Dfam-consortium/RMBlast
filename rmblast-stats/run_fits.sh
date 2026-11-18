#!/bin/sh
# Run ALP Gumbel fits for the RepeatMasker matrix set at canonical gap costs.
# Settings reproduce the original rmblast_*_values fits exactly:
# ALP defaults eps_lambda=0.01 eps_K=0.05, seed=1, iad=false, max_time=10.
#
# MATDIR (RepeatMasker's ncbi/nt matrices), CMPDIR (the directory holding
# comparison.matrix) and OUT (where the .fit/.par/.err files go, read by
# gen_tables.py) can be set in the environment.  Exits non-zero, naming
# them, if any fit fails.
cd "$(dirname "$0")" || exit 1
MATDIR=${MATDIR:-/usr/local/RepeatMasker/Matrices/ncbi/nt}
CMPDIR=${CMPDIR:-/usr/local/RepeatModeler/Matrices/ncbi/nt}
FIT=./alp-fit/alp_fit
OUT=${OUT:-./fits}
mkdir -p "$OUT" || exit 1
[ -x "$FIT" ] || { echo "run_fits.sh: $FIT not built" >&2; exit 1; }

failed=""
fit_one() {
    name=$1; path=$2; go=$3; ge=$4
    echo "=== $name @ $go/$ge"
    if [ ! -f "$path" ]; then
        echo "FAILED: $name ($path not found)"; failed="$failed $name"; return
    fi
    if ! $FIT -matrix "$path" -gapopen "$go" -gapextend "$ge" -max_time 10 \
            -out "$OUT/$name.par" > "$OUT/$name.fit" 2> "$OUT/$name.err"; then
        echo "FAILED: $name (see $OUT/$name.err)"; failed="$failed $name"; return
    fi
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

if [ -n "$failed" ]; then
    echo "run_fits.sh: failed:$failed" >&2
    exit 1
fi
echo DONE
