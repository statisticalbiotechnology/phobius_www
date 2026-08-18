"""Homology search for PolyPhobius.

The legacy server ran ``blastget``: legacy NCBI ``blastall`` against a local
UniProt/TrEMBL database, with hit sequences pulled back out of a BioPerl
``Bio::Index::Fasta`` index. That is not deployable on a platform with a 5 GB
volume cap, and it dragged in 721 vendored BioPerl files.

Modern BLAST+ replaces it: ``blastp`` for the search and ``blastdbcmd`` to pull
full subject sequences straight out of the same database, so no separate index
is needed.

DIAMOND was tried first and rejected on measurement. For a *single* query its
runtime is dominated by building a seed index over the whole database: against
Swiss-Prot it took 42.8 s where blastp takes 1.2 s, and it returned a strict
subset of blastp's hits. DIAMOND wins on bulk searches; this is the opposite
case.

Predictions from this path are **not** identical to the legacy ones: the
database is Swiss-Prot rather than TrEMBL. The "supply your own alignment" path
remains the reproducible one.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from .config import Settings, settings as default_settings
from .fasta import Record, format_fasta, parse

#: Legacy blastget thresholds (blastget: frac_aligned_hit / frac_aligned_query).
MIN_COVERAGE = 0.75
MAX_HITS = 50
EVALUE = "1e-5"


class HomologyError(RuntimeError):
    """The homology search failed or found too little to work with."""


def _run(cmd: list[str], timeout: int, what: str, cwd: Path | None = None) -> str:
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, check=False, cwd=cwd
        )
    except subprocess.TimeoutExpired as exc:
        raise HomologyError(f"{what} did not finish within {timeout}s.") from exc
    except FileNotFoundError as exc:
        raise HomologyError(f"{what}: executable not found ({cmd[0]}).") from exc
    if proc.returncode != 0:
        detail = (proc.stderr or "").strip().splitlines()
        raise HomologyError(f"{what} failed: {detail[-1] if detail else proc.returncode}")
    return proc.stdout


def search(query: Record, cfg: Settings | None = None) -> list[Record]:
    """Return the query followed by its homologues, unaligned."""
    cfg = cfg or default_settings
    if not any(Path(f"{cfg.blast_db}{ext}").is_file() for ext in (".pin", ".pal")):
        raise HomologyError(
            "No homology database is configured on this server. Submit an "
            "alignment instead, or use the standalone package."
        )

    with tempfile.TemporaryDirectory(prefix="phobius-search-") as tmp:
        query_path = Path(tmp) / "query.fa"
        query_path.write_text(format_fasta([query]))

        # qlen/slen/length let us reproduce blastget's frac_aligned_query and
        # frac_aligned_hit; BLAST has no subject-coverage output specifier.
        table = _run(
            [
                cfg.blastp,
                "-db", str(cfg.blast_db),
                "-query", str(query_path),
                "-evalue", EVALUE,
                "-max_target_seqs", str(MAX_HITS),
                "-num_threads", "2",
                "-outfmt", "6 sseqid qlen slen length",
            ],
            cfg.homology_timeout,
            "BLAST",
            cwd=Path(tmp),
        )

        ids: list[str] = []
        seen: set[str] = set()
        for line in table.splitlines():
            fields = line.split("\t")
            if len(fields) != 4:
                continue
            sseqid, qlen, slen, length = fields
            try:
                if int(length) / int(qlen) <= MIN_COVERAGE:
                    continue
                if int(length) / int(slen) <= MIN_COVERAGE:
                    continue
            except (ValueError, ZeroDivisionError):
                continue
            if sseqid in seen:
                continue
            seen.add(sseqid)
            ids.append(sseqid)

        if not ids:
            raise HomologyError(
                "No homologues passed the coverage threshold, so a "
                "homology-supported prediction would be no different from the "
                "plain one. Try the normal prediction instead."
            )

        # Retrieve the *full* subject sequences. The alignment must cover whole
        # proteins, not just the aligned segments BLAST reports.
        id_file = Path(tmp) / "ids.txt"
        id_file.write_text("\n".join(ids) + "\n")
        fasta = _run(
            [
                cfg.blastdbcmd,
                "-db", str(cfg.blast_db),
                "-entry_batch", str(id_file),
                "-outfmt", "%f",
            ],
            cfg.homology_timeout,
            "blastdbcmd",
            cwd=Path(tmp),
        )

    homologues = [r for r in parse(fasta) if r.name != query.name]
    if not homologues:
        raise HomologyError(
            "The homologues found could not be retrieved from the database. It "
            "may have been built without -parse_seqids."
        )
    return [query, *homologues]


def align(records: list[Record], cfg: Settings | None = None) -> list[Record]:
    """Align sequences with Kalign, keeping the query first.

    PolyPhobius predicts for whichever sequence comes first in the alignment, so
    the query is moved back to the front if the aligner reorders it.
    """
    cfg = cfg or default_settings
    if len(records) < 2:
        raise HomologyError("At least two sequences are needed to build an alignment.")

    with tempfile.TemporaryDirectory(prefix="phobius-align-") as tmp:
        infile = Path(tmp) / "in.fa"
        outfile = Path(tmp) / "out.fa"
        infile.write_text(format_fasta(records))
        _run(
            [cfg.kalign, "-i", str(infile), "-o", str(outfile), "-f", "fasta"],
            cfg.homology_timeout,
            "Kalign",
            cwd=Path(tmp),
        )
        aligned = parse(outfile.read_text())

    if not aligned:
        raise HomologyError("Kalign produced an empty alignment.")

    query_name = records[0].name
    ordered = [r for r in aligned if r.name == query_name]
    ordered += [r for r in aligned if r.name != query_name]
    if not ordered or ordered[0].name != query_name:
        raise HomologyError("The query sequence was lost during alignment.")

    lengths = {len(r) for r in ordered}
    if len(lengths) != 1:
        raise HomologyError("Kalign returned rows of differing length.")
    return ordered


def search_and_align(query: Record, cfg: Settings | None = None) -> list[Record]:
    """Full pipeline: find homologues, then align them with the query."""
    return align(search(query, cfg), cfg)
