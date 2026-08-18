#!/usr/bin/env bash
# Build the BLAST database used by the PolyPhobius homology search.
#
# Run this once against the mounted data volume; re-run to refresh. Swiss-Prot is
# used rather than TrEMBL because the whole database must fit in the platform's
# 5 GB volume (the result is roughly 340 MB).
#
# -parse_seqids is required: without it blastdbcmd cannot retrieve the full
# subject sequences that the alignment step needs.
set -euo pipefail

DEST="${1:-/home/data}"
MAKEBLASTDB="${MAKEBLASTDB:-makeblastdb}"
URL="https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

mkdir -p "$DEST"
echo "Downloading Swiss-Prot to $DEST ..."
curl -fSL --retry 3 -o "$DEST/uniprot_sprot.fasta.gz" "$URL"

echo "Building BLAST database ..."
zcat "$DEST/uniprot_sprot.fasta.gz" \
    | "$MAKEBLASTDB" -in - -dbtype prot -title swissprot \
        -out "$DEST/swissprot" -parse_seqids

rm -f "$DEST/uniprot_sprot.fasta.gz"
echo "Done: $DEST/swissprot.* ($(du -ch "$DEST"/swissprot.* | tail -1 | cut -f1))"
echo "No configuration needed if $DEST is a searched mount path; otherwise set"
echo "PHOBIUS_BLAST_DB=$DEST/swissprot"
