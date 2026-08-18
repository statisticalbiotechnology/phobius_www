"""Homology search pipeline (BLAST+ and Kalign).

These tests build a throwaway BLAST database rather than depending on Swiss-Prot
being present, so they run anywhere the tools are installed.
"""

import os
import random
import shutil
import subprocess

import pytest

from app.config import Settings
from app.fasta import Record
from app.homology import HomologyError, align, search

MAKEBLASTDB = shutil.which("makeblastdb")
BLASTP = shutil.which("blastp")
KALIGN = shutil.which("kalign")

pytestmark = pytest.mark.skipif(
    not (MAKEBLASTDB and BLASTP and KALIGN),
    reason="ncbi-blast+ and kalign are not installed",
)

QUERY = (
    "MYGKIIFVLLLSAIVSISASSTTGVAMHTSTSSSVTKSYISSQTNDTHKRDTYAATPRAHEVSEISVRT"
    "VYPPEEETGERVQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKKSPSDVKPLPSPDTDVPLS"
    "SVEIENPETSDQ"
)


@pytest.fixture(scope="module")
def tiny_db(tmp_path_factory):
    """A database of near-identical homologues plus unrelated decoys."""
    directory = tmp_path_factory.mktemp("blastdb")
    fasta = directory / "db.fasta"
    rng = random.Random(3)

    records = []
    for i in range(12):
        mutated = list(QUERY)
        for _ in range(int(len(mutated) * 0.12)):
            mutated[rng.randrange(len(mutated))] = rng.choice("ACDEFGHIKLMNPQRSTVWY")
        records.append(f">sp|TEST{i:02d}|HOMOLOG_{i}\n" + "".join(mutated))
    for i in range(5):
        decoy = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(300))
        records.append(f">sp|DEC{i:02d}|DECOY_{i}\n{decoy}")

    fasta.write_text("\n".join(records) + "\n")
    subprocess.run(
        [MAKEBLASTDB, "-in", str(fasta), "-dbtype", "prot", "-title", "test",
         "-out", str(directory / "swissprot"), "-parse_seqids"],
        check=True, capture_output=True, cwd=directory,
    )
    return directory / "swissprot"


@pytest.fixture
def cfg(tiny_db):
    return Settings(blast_db=tiny_db)


def test_search_finds_homologues_and_drops_decoys(cfg):
    hits = search(Record("query", QUERY), cfg)
    assert hits[0].name == "query"
    assert len(hits) > 1
    assert all("DECOY" not in r.name for r in hits[1:])


def test_search_returns_full_subject_sequences(cfg):
    """blastdbcmd retrieves whole proteins, not just the aligned segments.

    BLAST's own `sseq` output would give only the aligned region, which would
    truncate the homologues before they reach the aligner.
    """
    hits = search(Record("query", QUERY), cfg)
    assert all(len(r.sequence) == len(QUERY) for r in hits[1:])
    assert all("-" not in r.sequence for r in hits[1:])


def test_alignment_is_rectangular_with_query_first(cfg):
    aligned = align(search(Record("query", QUERY), cfg), cfg)
    assert aligned[0].name == "query"
    assert len({len(r) for r in aligned}) == 1


def test_pipeline_does_not_write_to_the_working_directory(cfg, tmp_path, monkeypatch):
    """Regression: a search tool defaulting to scratch files in the working
    directory fails when the container runs as UID 1000 against a read-only /app."""
    workdir = tmp_path / "readonly"
    workdir.mkdir()
    monkeypatch.chdir(workdir)
    os.chmod(workdir, 0o555)
    try:
        aligned = align(search(Record("query", QUERY), cfg), cfg)
        assert aligned
        assert list(workdir.iterdir()) == []
    finally:
        os.chmod(workdir, 0o755)


def test_missing_database_gives_a_useful_message(tmp_path):
    cfg = Settings(blast_db=tmp_path / "nope")
    with pytest.raises(HomologyError, match="No homology database"):
        search(Record("query", QUERY), cfg)


def test_no_homologues_is_reported_clearly(cfg):
    unrelated = Record("q", "WWWWWWWWWWCCCCCCCCCCWWWWWWWWWWCCCCCCCCCCWWWWWWWWWW")
    with pytest.raises(HomologyError, match="No homologues|coverage"):
        search(unrelated, cfg)


def test_database_without_parse_seqids_is_reported(tmp_path):
    """Without -parse_seqids the search succeeds but retrieval returns nothing."""
    fasta = tmp_path / "db.fasta"
    fasta.write_text(f">sp|X|HOMOLOG\n{QUERY}\n")
    subprocess.run(
        [MAKEBLASTDB, "-in", str(fasta), "-dbtype", "prot", "-out", str(tmp_path / "swissprot")],
        check=True, capture_output=True, cwd=tmp_path,
    )
    cfg = Settings(blast_db=tmp_path / "swissprot")
    with pytest.raises(HomologyError):
        search(Record("query", QUERY), cfg)
