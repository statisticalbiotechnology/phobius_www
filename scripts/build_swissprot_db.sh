#!/usr/bin/env bash
# Build the DIAMOND database used by the PolyPhobius homology search.
#
# Run this once against the mounted data volume; re-run to refresh. Swiss-Prot
# is used rather than TrEMBL because the whole index must fit in the platform's
# 5 GB volume (the resulting .dmnd is roughly 250 MB).
set -euo pipefail

DEST="${1:-/mnt/data}"
DIAMOND="${DIAMOND:-diamond}"
URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

mkdir -p "$DEST"
echo "Downloading Swiss-Prot to $DEST ..."
curl -fSL --retry 3 -o "$DEST/uniprot_sprot.fasta.gz" "$URL"

echo "Building DIAMOND database ..."
"$DIAMOND" makedb --in "$DEST/uniprot_sprot.fasta.gz" --db "$DEST/swissprot" --quiet

rm -f "$DEST/uniprot_sprot.fasta.gz"
echo "Done: $DEST/swissprot.dmnd ($(du -h "$DEST/swissprot.dmnd" | cut -f1))"
echo "Set PHOBIUS_DIAMOND_DB=$DEST/swissprot.dmnd"
