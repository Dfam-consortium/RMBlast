#!/bin/bash
# Tests for the wrapper scripts in wrappers/ (and the verify_wrappers/ rmblastn).
#
#   bash tests/wrappers.sh            # uses ../target/release binaries
#   RMBLAST_BIN_DIR=<dir> bash tests/wrappers.sh
#
# Where NCBI RMBlast 2.17.1 is installed (NCBI_DIR, default
# /usr/local/rmblast-2.17.1/bin) the blastdbcmd cases are also compared with
# NCBI's blastdbcmd, and the verify wrapper is exercised; otherwise those
# checks are skipped and say so.  The UCSC tools (faToTwoBit, twoBitToFa,
# twoBitInfo) are needed for the 2bit cases; they are looked up the way the
# wrappers do (UCSCTOOLS_DIR, then PATH).
#
# Exit status: 0 if every check passed, 1 otherwise.

set -u
HERE="$(cd "$(dirname "$0")/.." && pwd)"
W="$HERE/wrappers"
VW="$HERE/verify_wrappers"
DATA="$HERE/rmblastn/tests/data"
NCBI_DIR="${NCBI_DIR:-/usr/local/rmblast-2.17.1/bin}"
T="$(mktemp -d "${TMPDIR:-/tmp}/rmblast-wrapper-tests.XXXXXX")"
trap 'rm -rf "$T"' EXIT

pass=0; fail=0; skip=0
ok()   { pass=$((pass+1)); echo "ok   $1"; }
bad()  { fail=$((fail+1)); echo "FAIL $1"; [ -n "${2:-}" ] && printf '     %s\n' "$2"; }
skp()  { skip=$((skip+1)); echo "skip $1 ($2)"; }

# run <name> <expected rc> <stderr substring or ""> <command...>
# Leaves stdout in $T/out and stderr in $T/err.
run() {
    local name="$1" want_rc="$2" want_err="$3"; shift 3
    "$@" > "$T/out" 2> "$T/err"; local rc=$?
    if [ "$rc" != "$want_rc" ]; then
        bad "$name" "exit $rc, expected $want_rc; stderr: $(head -c 300 "$T/err")"
    elif [ -n "$want_err" ] && ! grep -qF -- "$want_err" "$T/err"; then
        bad "$name" "stderr lacks '$want_err': $(head -c 300 "$T/err")"
    else
        ok "$name"
    fi
}

have_ucsc=0
for t in faToTwoBit twoBitToFa twoBitInfo; do
    [ -x "${UCSCTOOLS_DIR:-/usr/local/ucscTools}/$t" ] || command -v "$t" > /dev/null || { have_ucsc=0; break; }
    have_ucsc=1
done
have_ncbi=0; [ -x "$NCBI_DIR/blastdbcmd" ] && [ -x "$NCBI_DIR/makeblastdb" ] && have_ncbi=1

# A database with ids fam1..fam17 (the CLI fixture records, renamed).
awk '/^>/ { n++; sub(/^>/, ">fam" n " ") } { print }' "$DATA/db.fa" > "$T/lib.fa"

# ── makeblastdb ─────────────────────────────────────────────────────────────
mkdir "$T/mk"
run "makeblastdb: RepeatMasker-style in-place build" 0 "Done" \
    "$W/makeblastdb" -dbtype nucl -parse_seqids -blastdb_version 4 -out "$T/lib.fa" -in "$T/lib.fa"
for ext in nhr nin nsq njs; do
    [ -s "$T/lib.fa.$ext" ] || bad "makeblastdb: dummy .$ext written"
done
if [ $have_ucsc = 1 ]; then
    [ -s "$T/lib.fa.2bit" ] && ok "makeblastdb: .2bit written" || bad "makeblastdb: .2bit written"
fi
run "makeblastdb: -out elsewhere copies the FASTA" 0 "copy of" \
    "$W/makeblastdb" -out "$T/mk/db" -in "$T/lib.fa" -dbtype nucl
cmp -s "$T/lib.fa" "$T/mk/db" && ok "makeblastdb: copy is identical" || bad "makeblastdb: copy is identical"
ls "$T/mk" | grep -q '\.tmp\.' && bad "makeblastdb: no temp files left" || ok "makeblastdb: no temp files left"
run "makeblastdb: unknown option refused" 1 "not supported by this wrapper: -bogus" \
    "$W/makeblastdb" -in "$T/lib.fa" -out "$T/mk/x" -bogus
run "makeblastdb: -dbtype prot refused" 1 "nucleotide only" \
    "$W/makeblastdb" -in "$T/lib.fa" -out "$T/mk/p" -dbtype prot
[ -e "$T/mk/p.nhr" ] || [ -e "$T/mk/p.phr" ] && bad "makeblastdb: prot writes nothing" || ok "makeblastdb: prot writes nothing"
run "makeblastdb: -mask_data refused" 1 "-mask_data" \
    "$W/makeblastdb" -in "$T/lib.fa" -out "$T/mk/m" -mask_data x.asnb -dbtype nucl
run "makeblastdb: -version" 0 "" "$W/makeblastdb" -version
grep -q "makeblastdb" "$T/out" && ok "makeblastdb: -version names itself" || bad "makeblastdb: -version names itself"
# Without faToTwoBit anywhere: still a usable database, with a warning.
mkdir "$T/nopath"
for t in bash cp mv rm cat sed awk dirname basename; do ln -s "$(command -v $t)" "$T/nopath/$t"; done
if [ -x /usr/local/ucscTools/faToTwoBit ]; then
    skp "makeblastdb: works without faToTwoBit" "/usr/local/ucscTools/faToTwoBit exists (the last fallback)"
else
    run "makeblastdb: works without faToTwoBit" 0 "faToTwoBit not found" \
        env -u UCSCTOOLS_DIR PATH="$T/nopath" "$W/makeblastdb" -in "$T/lib.fa" -out "$T/mk/no2bit" -dbtype nucl
    [ -e "$T/mk/no2bit.2bit" ] && bad "makeblastdb: no .2bit without faToTwoBit" || ok "makeblastdb: no .2bit without faToTwoBit"
fi

# ── blastdbcmd (FASTA backend, and 2bit when the UCSC tools are there) ─────
backends="fasta"; [ $have_ucsc = 1 ] && backends="fasta 2bit"
for be in $backends; do
    if [ "$be" = fasta ]; then db="$T/fa/lib"; mkdir -p "$T/fa"; cp "$T/lib.fa" "$db"
    else db="$T/tb/lib"; mkdir -p "$T/tb"; "$W/makeblastdb" -in "$T/lib.fa" -out "$T/tb/lib.fa" -dbtype nucl 2>/dev/null
         mv "$T/tb/lib.fa.2bit" "$T/tb/lib.2bit"; rm -f "$T/tb/lib.fa"*; fi
    run "blastdbcmd[$be]: -entry fam3" 0 "" "$W/blastdbcmd" -db "$db" -entry fam3
    grep -q '^>fam3' "$T/out" && ok "blastdbcmd[$be]: fam3 header" || bad "blastdbcmd[$be]: fam3 header"
    cp "$T/out" "$T/fam3.$be"
    run "blastdbcmd[$be]: -entry nosuch is an error" 1 "Entry not found: nosuch" \
        "$W/blastdbcmd" -db "$db" -entry nosuch
    run "blastdbcmd[$be]: flags take no value" 0 "" \
        "$W/blastdbcmd" -db "$db" -target_only -entry all -outfmt "%i %l"
    [ "$(wc -l < "$T/out")" = 17 ] && ok "blastdbcmd[$be]: 17 ids listed" || bad "blastdbcmd[$be]: 17 ids listed" "$(wc -l < "$T/out") lines"
    printf 'fam3 1-50\nfam3 900000-900100\nnosuch 1-10\nfam2 10-20\nfam3 100-90000\n' > "$T/batch"
    run "blastdbcmd[$be]: batch with bad entries exits 1" 1 "Skipped nosuch" \
        "$W/blastdbcmd" -db "$db" -entry_batch "$T/batch"
    grep -qF "Skipped fam3: interval 900000-900100 is outside" "$T/err" \
        && ok "blastdbcmd[$be]: start past the end skipped" || bad "blastdbcmd[$be]: start past the end skipped" "$(cat "$T/err")"
    [ "$(grep -c '^>' "$T/out")" = 3 ] && ok "blastdbcmd[$be]: the 3 good intervals printed" \
        || bad "blastdbcmd[$be]: the 3 good intervals printed" "$(grep '^>' "$T/out")"
    cp "$T/out" "$T/batch.$be"
    echo "keep me" > "$T/outfile"
    run "blastdbcmd[$be]: bad -outfmt refused" 1 "not supported" \
        "$W/blastdbcmd" -db "$db" -entry all -outfmt "%s" -out "$T/outfile"
    [ "$(cat "$T/outfile")" = "keep me" ] && ok "blastdbcmd[$be]: -out untouched on error" || bad "blastdbcmd[$be]: -out untouched on error"
    run "blastdbcmd[$be]: -info" 0 "" "$W/blastdbcmd" -db "$db" -info
    grep -q "^	17 sequences; " "$T/out" && ok "blastdbcmd[$be]: -info counts" || bad "blastdbcmd[$be]: -info counts" "$(head -3 "$T/out")"
done
# FASTA-only edge cases
printf '\n>  spaced desc\nACGTACGTAC\n>b\nAAAA\n' > "$T/sp.fa"
run "blastdbcmd: leading blank line, '> name' defline" 0 "" "$W/blastdbcmd" -db "$T/sp.fa" -entry spaced
grep -q '^ACGTACGTAC$' "$T/out" && ok "blastdbcmd: '> name' record found" || bad "blastdbcmd: '> name' record found" "$(cat "$T/out")"
: > "$T/empty.fa"; echo ">x" > /dev/null
run "blastdbcmd: -info on an empty FASTA" 0 "" "$W/blastdbcmd" -db "$T/empty.fa" -info
grep -q "0 sequences; 0 total bases" "$T/out" && ok "blastdbcmd: empty -info prints 0s" || bad "blastdbcmd: empty -info prints 0s" "$(head -3 "$T/out")"

# Against NCBI: the same bases for -entry and for the batch intervals.
if [ $have_ncbi = 1 ]; then
    mkdir "$T/ncbi"; cp "$T/lib.fa" "$T/ncbi/lib.fa"
    "$NCBI_DIR/makeblastdb" -in "$T/ncbi/lib.fa" -dbtype nucl -parse_seqids -out "$T/ncbi/lib" > /dev/null
    "$NCBI_DIR/blastdbcmd" -db "$T/ncbi/lib" -entry fam3 > "$T/fam3.ncbi"
    printf 'fam3 1-50\nfam2 10-20\nfam3 100-90000\n' > "$T/goodbatch"
    "$NCBI_DIR/blastdbcmd" -db "$T/ncbi/lib" -entry_batch "$T/goodbatch" > "$T/batch.ncbi"
    bases() { grep -v '^>' "$1" | tr -d '\n'; }
    for be in $backends; do
        [ "$(bases "$T/fam3.$be")" = "$(bases "$T/fam3.ncbi")" ] && ok "blastdbcmd[$be]: -entry bases = NCBI" || bad "blastdbcmd[$be]: -entry bases = NCBI"
        [ "$(bases "$T/batch.$be")" = "$(bases "$T/batch.ncbi")" ] && ok "blastdbcmd[$be]: batch bases = NCBI" || bad "blastdbcmd[$be]: batch bases = NCBI"
    done
    "$NCBI_DIR/blastdbcmd" -db "$T/ncbi/lib" -entry nosuch > /dev/null 2>&1; nrc=$?
    [ $nrc = 1 ] && ok "blastdbcmd: NCBI also exits 1 for a missing entry" || bad "blastdbcmd: NCBI exit for a missing entry" "$nrc"
else
    skp "blastdbcmd: comparison with NCBI" "no NCBI blastdbcmd at $NCBI_DIR"
fi

# ── rmblastn wrapper ────────────────────────────────────────────────────────
RM_ARGS=(-num_alignments 9999999 -gapopen 20 -gapextend 5 -mask_level 80 -complexity_adjust
         -word_size 14 -xdrop_ungap 400 -xdrop_gap_final 200 -xdrop_gap 100 -min_raw_gapped_score 200
         -dust yes -outfmt "6 score qseqid qstart qend sstrand sseqid sstart send" -num_threads 2)
BIN="${RMBLAST_BIN_DIR:-$HERE/target/release}/rmblastn"
run "rmblastn wrapper: RepeatMasker-style call" 0 "" \
    env BLASTMAT="$DATA" "$W/rmblastn" -db "$T/lib.fa" -query "$DATA/query.fa" -matrix comparison.matrix "${RM_ARGS[@]}"
cp "$T/out" "$T/rm.wrapper"
"$BIN" --db "$T/lib.fa" --query "$DATA/query.fa" --matrix "$DATA/comparison.matrix" \
    "${RM_ARGS[@]:2}" > "$T/rm.direct" 2>/dev/null
cmp -s "$T/rm.wrapper" "$T/rm.direct" && ok "rmblastn wrapper: output = binary's" || bad "rmblastn wrapper: output = binary's"
run "rmblastn wrapper: .2bit db name maps to the FASTA" 0 "" \
    env BLASTMAT="$DATA" "$W/rmblastn" -db "$T/lib.fa.2bit" -query "$DATA/query.fa" -matrix comparison.matrix "${RM_ARGS[@]}"
cmp -s "$T/out" "$T/rm.direct" && ok "rmblastn wrapper: .2bit name, same output" || bad "rmblastn wrapper: .2bit name, same output"
# Negative values pass through (3.0.12: "invalid value '--4.9' for '--matrix-beta'")
run "rmblastn wrapper: negative values pass through" 0 "" \
    env BLASTMAT="$DATA" "$W/rmblastn" -db "$T/lib.fa" -query "$DATA/query.fa" -matrix comparison.matrix \
        -gapopen 20 -gapextend 5 -word_size 14 -matrix_beta -4.9 -mask_level -1 -outfmt "6 score"
run "rmblastn wrapper: -version in NCBI's shape" 0 "" "$W/rmblastn" -version
grep -Eq '^rmblastn: [0-9]+\.[0-9]+\.[0-9]+\+$' <(head -1 "$T/out") && ok "rmblastn wrapper: 'rmblastn: X.Y.Z+'" \
    || bad "rmblastn wrapper: 'rmblastn: X.Y.Z+'" "$(head -1 "$T/out")"
# The wrapper's list of boolean flags must match the binary's.
# (only the "Options:" section, which ends at the first unindented line)
"$BIN" --help | awk '/^Options:/ { on = 1; next } on && /^[^ ]/ { on = 0 }
        on && /^ +(-[a-zA-Z], )?--[a-z]/ && !/</ { for (i = 1; i <= NF; i++) if ($i ~ /^--/) { sub(/^--/, "", $i); print $i } }' \
    | grep -vx 'help\|version' | sort > "$T/flags.bin"
sed -n 's/^RMBLASTN_FLAGS=" \(.*\) "$/\1/p' "$W/rmblastn" | tr ' ' '\n' | sort > "$T/flags.wrapper"
cmp -s "$T/flags.bin" "$T/flags.wrapper" && ok "rmblastn wrapper: flag list = binary's" \
    || bad "rmblastn wrapper: flag list = binary's" "binary: $(tr '\n' ' ' < "$T/flags.bin") wrapper: $(tr '\n' ' ' < "$T/flags.wrapper")"
grep -q "^RMBLASTN_FLAGS=\"$(sed -n 's/^RMBLASTN_FLAGS="\(.*\)"$/\1/p' "$W/rmblastn")\"" "$VW/rmblastn" \
    && ok "verify rmblastn: same flag list" || bad "verify rmblastn: same flag list"

# ── verify wrapper (needs NCBI) ─────────────────────────────────────────────
if [ -x "$NCBI_DIR/rmblastn" ]; then
    cp "$T/lib.fa" "$T/v.fa"
    NCBI_DIR="$NCBI_DIR" "$VW/makeblastdb" -in "$T/v.fa" -out "$T/v.fa" -dbtype nucl -parse_seqids > /dev/null 2>&1
    mkdir "$T/tmp"
    run "verify rmblastn: RepeatMasker-style call" 0 "Rust MATCHES NCBI" \
        env BLASTMAT="$DATA" NCBI_DIR="$NCBI_DIR" RMBLAST_VERIFY_DIR="$T/vd" TMPDIR="$T/tmp" \
        "$VW/rmblastn" -db "$T/v.fa" -query "$DATA/query.fa" -matrix comparison.matrix "${RM_ARGS[@]}" -matrix_beta -4.9
    run "verify rmblastn: query on stdin" 0 "Rust MATCHES NCBI" \
        env BLASTMAT="$DATA" NCBI_DIR="$NCBI_DIR" RMBLAST_VERIFY_DIR="$T/vd" TMPDIR="$T/tmp" \
        bash -c '"$0" -db "$1" -query - -matrix comparison.matrix "${@:3}" < "$2"' \
        "$VW/rmblastn" "$T/v.fa" "$DATA/query.fa" "${RM_ARGS[@]}"
    [ -s "$T/out" ] && ok "verify rmblastn: stdin query gave hits" || bad "verify rmblastn: stdin query gave hits"
    ls "$T/tmp"/verify_rmblastn.* > /dev/null 2>&1 && bad "verify rmblastn: no temp dirs left" "$(ls -d "$T/tmp"/verify_rmblastn.* | head -3)" \
        || ok "verify rmblastn: no temp dirs left"
    # RepeatMasker 4.2.6+ drives a 3.x rmblastn with GNU options and a
    # matrix path (RMBlastSearchEngine.pm); the verify wrapper must hand NCBI
    # its own spelling.
    run "verify rmblastn: RepeatMasker 4.2.6 GNU-style call" 0 "Rust MATCHES NCBI" \
        env -u BLASTMAT NCBI_DIR="$NCBI_DIR" RMBLAST_VERIFY_DIR="$T/vd" TMPDIR="$T/tmp" \
        "$VW/rmblastn" --db "$T/v.fa" --query "$DATA/query.fa" --gapopen 20 --gapextend 5 \
            --mask-level 80 --complexity-adjust --word-size 14 --xdrop-ungap 400 --xdrop-gap-final 200 \
            --xdrop-gap 100 --min-raw-gapped-score 200 --dust no --num-threads 2 \
            --outfmt="6 score qseqid qstart qend sstrand sseqid sstart send" --matrix "$DATA/comparison.matrix"
else
    skp "verify rmblastn" "no NCBI rmblastn at $NCBI_DIR"
fi

echo
echo "$pass passed, $fail failed, $skip skipped"
[ $fail = 0 ]
