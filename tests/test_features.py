"""Feature rendering, checked against the format documented for the old server."""

from app.features import (
    SHORT_HEADER,
    count_transmembrane,
    has_signal_peptide,
    label_runs,
    long_format,
    short_format,
    topology_string,
)

# The worked example from the legacy instructions page.
MTH_DROME = (
    "nnn" + "h" * 16 + "ccccc" + "C"
    + "o" * 194 + "M" * 20 + "i" * 11 + "M" * 20 + "o" * 11
)


def test_label_runs_are_one_based_and_inclusive():
    runs = label_runs("nnhhMMM")
    assert [(r.label, r.start, r.stop) for r in runs] == [
        ("n", 1, 2), ("h", 3, 4), ("M", 5, 7)
    ]
    assert len(runs[2]) == 3


def test_long_format_matches_documented_example():
    lines = long_format("MTH_DROMEa signal peptide", label_runs(MTH_DROME)).splitlines()
    assert lines[0] == "ID   MTH_DROMEa signal peptide"
    # SIGNAL is emitted before the regions it summarises, as in the original.
    assert lines[1] == "FT   SIGNAL        1     24       "
    assert lines[2] == "FT   REGION        1      3       N-REGION."
    assert lines[3] == "FT   REGION        4     19       H-REGION."
    assert lines[4] == "FT   REGION       20     24       C-REGION."
    assert lines[5] == "FT   TOPO_DOM     25    219       NON CYTOPLASMIC."
    assert lines[6] == "FT   TRANSMEM    220    239       "
    assert lines[-1] == "//"


def test_topology_string_matches_documented_example():
    assert topology_string(label_runs(MTH_DROME)).startswith("n4-19c24/25o220-239i")


def test_short_format_columns():
    line = short_format("MTH_DROME", label_runs(MTH_DROME))
    assert line.startswith("MTH_DROME" + " " * 21)
    assert line[30:36] == "  2  Y"
    assert SHORT_HEADER.startswith("SEQENCE ID")  # legacy spelling, kept on purpose


def test_soluble_protein_reported_as_cytoplasmic():
    # predict.pl:397 rewrote a wholly non-cytoplasmic topology to "i".
    assert topology_string(label_runs("o" * 50)) == "i"


def test_two_region_prediction_keeps_legacy_cytoplasmic_quirk():
    # predict.pl:343: with exactly two regions an 'o' is labelled CYTOPLASMIC.
    two = long_format("x", label_runs("M" * 10 + "o" * 10))
    assert "CYTOPLASMIC." in two and "NON CYTOPLASMIC." not in two
    three = long_format("x", label_runs("i" * 5 + "M" * 10 + "o" * 10))
    assert "NON CYTOPLASMIC." in three


def test_counts():
    runs = label_runs(MTH_DROME)
    assert count_transmembrane(runs) == 2
    assert has_signal_peptide(runs) is True
    assert has_signal_peptide(label_runs("i" * 20)) is False
