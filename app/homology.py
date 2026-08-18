"""Homology search for PolyPhobius.

The legacy server ran ``blastget``: legacy NCBI ``blastall`` against a local
UniProt/TrEMBL database, with hit sequences pulled back out of a BioPerl
``Bio::Index::Fasta`` index. That is not deployable on a platform with a 5 GB
volume cap, and it dragged in 721 vendored BioPerl files.

DIAMOND replaces all of it. ``--outfmt 6 ... full_sseq`` returns the subject
sequence inline, so the separate index and fetch step disappear entirely.

Predictions from this path are **not** identical to the legacy ones: the
database is Swiss-Prot rather than TrEMBL and the search sensitivity model is
different. The "supply your own alignment" path remains the reproducible one.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from .config import Settings, settings as default_settings
from .fasta import Record, format_fasta, parse

#: Legacy blastget thresholds (blastget:  frac_aligned_hit / frac_aligned_query).
MIN_COVERAGE = 75.0
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
    if not cfg.diamond_db.exists():
        raise HomologyError(
            "No homology database is configured on this server. Submit an "
            "alignment instead, or use the standalone package."
        )

    with tempfile.TemporaryDirectory(prefix="phobius-search-") as tmp:
        query_path = Path(tmp) / "query.fa"
        query_path.write_text(format_fasta([query]))

        table = _run(
            [
                cfg.diamond, "blastp",
                "--db", str(cfg.diamond_db),
                "--query", str(query_path),
                "--very-sensitive",
                "--evalue", EVALUE,
                "--max-target-seqs", str(MAX_HITS),
                "--outfmt", "6", "sseqid", "qcovhsp", "scovhsp", "full_sseq",
                "--threads", "1",
                "--quiet",
                # DIAMOND writes scratch files to the working directory by
                # default, which is not writable when running as a non-root
                # user in the container. Keep them with the query instead.
                "--tmpdir", tmp,
            ],
            cfg.homology_timeout,
            "DIAMOND",
            cwd=Path(tmp),
        )

    homologues: list[Record] = []
    seen: set[str] = set()
    for line in table.splitlines():
        fields = line.split("\t")
        if len(fields) != 4:
            continue
        sseqid, qcov, scov, sequence = fields
        try:
            if float(qcov) <= MIN_COVERAGE or float(scov) <= MIN_COVERAGE:
                continue
        except ValueError:
            continue
        if sseqid in seen or sseqid == query.name:
            continue
        seen.add(sseqid)
        homologues.append(Record(sseqid, sequence.replace("-", "").upper()))

    if not homologues:
        raise HomologyError(
            "No homologues passed the coverage threshold, so a homology-supported "
            "prediction would be no different from the plain one. Try the normal "
            "prediction instead."
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
