#!/usr/bin/env bash
# Build the academic-use Phobius bundle offered on the download page.
#
# The bundle is the standalone command-line predictor: phobius.pl plus the
# decoder, model and decoder options. It is licensed, so it is NOT part of the
# container image -- place the result on the mounted storage, in the download/
# sub-directory of a searched mount (e.g. /home/data/download).
#
#   ./scripts/build_phobius_bundle.sh ~/bin /home/data/download
#
# Do not distribute the old web-server tree as "Phobius": it is 90 MB of
# obsolete JREs and vendored BioPerl, and it contains the legacy predict.pl CGI
# whose FASTA-header handling allowed command injection.
set -euo pipefail

SRC="${1:-$HOME/bin}"
DEST="${2:-/home/data/download}"
VERSION="${VERSION:-1.01}"

need() {
    [ -f "$SRC/$1" ] || { echo "missing: $SRC/$1" >&2; exit 1; }
}
need phobius.pl
need phobius.model
need phobius.options

# phobius.pl invokes "$PHOBIUS_DIR/decodeanhmm", so the binary must carry that
# exact name inside the bundle whatever it is called on disk.
BINARY=""
for candidate in decodeanhmm decodeanhmm.64bit; do
    [ -f "$SRC/$candidate" ] && { BINARY="$candidate"; break; }
done
[ -n "$BINARY" ] || { echo "missing: $SRC/decodeanhmm" >&2; exit 1; }

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
stage="$work/phobius"
mkdir -p "$stage"

install -m 755 "$SRC/phobius.pl"       "$stage/phobius.pl"
install -m 755 "$SRC/$BINARY"          "$stage/decodeanhmm"
install -m 644 "$SRC/phobius.model"    "$stage/phobius.model"
install -m 644 "$SRC/phobius.options"  "$stage/phobius.options"

cat > "$stage/README" <<'README'
Phobius - a combined transmembrane topology and signal peptide predictor

Usage:
    ./phobius.pl [-short] [-plp] sequences.fasta

Requires Perl 5. The decodeanhmm binary is a 64-bit Linux build; rebuild it
from the anhmm source for other platforms.

Please cite:
    Kall L, Krogh A, Sonnhammer ELL. A combined transmembrane topology and
    signal peptide prediction method. J Mol Biol 338:1027-1036 (2004).

LICENCE
    Phobius is free for academic use. It may not be redistributed or used for
    commercial purposes without permission from the authors. The model file
    (phobius.model) and the decoder (decodeanhmm) are covered by this licence.
README

mkdir -p "$DEST"
tarball="$DEST/phobius-${VERSION}.tar.gz"
tar czf "$tarball" -C "$work" phobius

echo "Wrote $tarball ($(du -h "$tarball" | cut -f1))"
tar tzf "$tarball" | sed 's/^/  /'
