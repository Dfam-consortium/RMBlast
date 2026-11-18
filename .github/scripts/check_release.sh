#!/usr/bin/env bash
# Pre-release check, run by `make check-release [TAG=X.Y.Z]` after a release
# build.  Fails unless
#   - the tracked tree has no uncommitted changes,
#   - Cargo.toml's version equals TAG (if given; a leading v is allowed),
#   - an existing tag TAG points at HEAD,
#   - CHANGELOG.md has a section for the version, and
#   - each binary in target/release reports exactly "<name> <version> (<HEAD>)".
# The last check is what would have caught 3.0.11's target/release, which
# reported 3.0.11 while holding older code.
#
# Usage: check_release.sh <version> <tag or ""> <binary>...
set -euo pipefail

version=$1 tag=$2
shift 2

fail() { echo "check-release: $*" >&2; exit 1; }

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    git status --short --untracked-files=no >&2
    fail "uncommitted changes (above)"
fi

if [ -n "$tag" ]; then
    [ "${tag#v}" = "$version" ] || fail "TAG $tag but Cargo.toml says $version"
    if git rev-parse -q --verify "refs/tags/$tag" >/dev/null; then
        [ "$(git rev-parse "$tag^{commit}")" = "$(git rev-parse HEAD)" ] \
            || fail "tag $tag does not point at HEAD"
    fi
fi

bash "$(dirname "$0")/release_notes.sh" "$version" >/dev/null \
    || fail "no CHANGELOG.md section for $version"

commit=$(git rev-parse --short=7 HEAD)
for b in "$@"; do
    got=$("target/release/$b" --version)
    [ "$got" = "$b $version ($commit)" ] \
        || fail "$b reports '$got', expected '$b $version ($commit)'"
done

echo "check-release: $version at $commit is ready to tag"
