"""Homology search pipeline (DIAMOND + Kalign).

These tests build a throwaway DIAMOND database rather than depending on
Swiss-Prot being present, so they run anywhere the two tools are installed.
"""

import os
import random
import shutil
import subprocess

import pytest

from app.config import Settings
from app.fasta import Record
from app.homology import HomologyError, align, search

DIAMOND = shutil.which("diamond")
KALIGN = shutil.which("kalign")

pytestmark = pytest.mark.skipif(
    not (DIAMOND and KALIGN), reason="diamond and kalign are not installed"
)

QUERY = (
    "MYGKIIFVLLLSAIVSISASSTTGVAMHTSTSSSVTKSYISSQTNDTHKRDTYAATPRAHEVSEISVRT"
    "VYPPEEETGERVQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKKSPSDVKPLPSPDTDVPLS"
    "SVEIENPETSDQ"
)


@pytest.fixture(scope="module")
def tiny_db(tmp_path_factory):
    """A database of near-identical homologues plus unrelated decoys."""
    directory = tmp_path_factory.mktemp("diamonddb")
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
        [DIAMOND, "makedb", "--in", str(fasta), "--db", str(directory / "test"), "--quiet"],
        check=True, cwd=directory,
    )
    return directory / "test.dmnd"


@pytest.fixture
def cfg(tiny_db):
    return Settings(diamond_db=tiny_db)


def test_search_finds_homologues_and_drops_decoys(cfg):
    hits = search(Record("query", QUERY), cfg)
    assert hits[0].name == "query"
    assert len(hits) > 1
    assert all("DECOY" not in r.name for r in hits)


def test_search_returns_full_subject_sequences(cfg):
    # --outfmt full_sseq is what let us delete the BioPerl fasta index entirely.
    hits = search(Record("query", QUERY), cfg)
    assert all(len(r.sequence) > 100 for r in hits[1:])
    assert all("-" not in r.sequence for r in hits[1:])


def test_alignment_is_rectangular_with_query_first(cfg):
    aligned = align(search(Record("query", QUERY), cfg), cfg)
    assert aligned[0].name == "query"
    assert len({len(r) for r in aligned}) == 1


def test_pipeline_does_not_write_to_the_working_directory(cfg, tmp_path, monkeypatch):
    """Regression: DIAMOND defaults to scratch files in the working directory,
    which fails when the container runs as UID 1000 against a read-only /app."""
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
    cfg = Settings(diamond_db=tmp_path / "nope.dmnd")
    with pytest.raises(HomologyError, match="No homology database"):
        search(Record("query", QUERY), cfg)


def test_no_homologues_is_reported_clearly(cfg):
    unrelated = Record("q", "WWWWWWWWWWCCCCCCCCCCWWWWWWWWWWCCCCCCCCCCWWWWWWWWWW")
    with pytest.raises(HomologyError, match="No homologues|coverage"):
        search(unrelated, cfg)
