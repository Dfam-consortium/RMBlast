#!/bin/bash
# Print a version's CHANGELOG.md section as GitHub release notes.
# Usage: release_notes.sh <version> [changelog]
#
# The section runs from "## <version> ..." to the next "## " heading, without the
# heading itself.  Hard-wrapped lines are joined, because GitHub shows each line
# break in a release description as a break on the page.  Exits non-zero if the
# section is missing or empty.
set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "usage: $0 <version> [changelog]" >&2
    exit 2
fi
version="$1"
changelog="${2:-CHANGELOG.md}"

notes="$(awk -v version="$version" '
    function flush() { if (buf != "") print buf; buf = "" }
    function gap() { if (blank) { print ""; blank = 0 } printed = 1 }
    /^## / && !in_fence { if (in_section) exit; in_section = ($2 == version); next }
    !in_section { next }
    # fenced code: copied verbatim, never joined
    /^[[:space:]]*(```|~~~)/ { flush(); gap(); print; in_fence = !in_fence; next }
    in_fence { print; next }
    /^[[:space:]]*$/ { flush(); if (printed) blank = 1; next }
    # headings and table rows stand alone: nothing joins onto them
    /^#/ || /^[[:space:]]*\|/ { flush(); gap(); print; next }
    {
        starts_block = ($0 ~ /^[[:space:]]*([-*+]|[0-9]+\.)[[:space:]]/)
        if (buf != "" && !starts_block) {
            line = $0
            sub(/^[[:space:]]+/, "", line)
            buf = buf " " line
            next
        }
        flush()
        gap()
        buf = $0
    }
    END { flush() }
' "$changelog")"

if [ -z "$notes" ]; then
    echo "no \"## $version\" section with any content in $changelog" >&2
    exit 1
fi
printf '%s\n' "$notes"
