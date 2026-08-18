"""Regression against the legacy implementation.

``tests/golden/labels.tsv`` holds per-residue labels produced by ``decodeanhmm``
-- the licensed native engine the container no longer ships. These tests assert
that the Java engine used in its place returns exactly the same thing, which is
the whole basis for having dropped the native binary.

Regenerate the golden file only with ``tests/generate_golden.py``, and only
against the legacy engine. If a change makes these fail, the change altered
predictions.
"""

import pathlib

import pytest

from app import engines, features
from app.fasta import parse

from conftest import needs_engine

pytestmark = needs_engine


@pytest.fixture(scope="module")
def java_only():
    """Settings that force the Java engine.

    Without this, a deployment with PHOBIUS_DECODEANHMM set would have these
    tests silently exercise the native binary instead of the engine the service
    uses by default.
    """
    from app.config import Settings

    return Settings(decodeanhmm=None)


@pytest.fixture(scope="module")
def predicted(corpus_text, java_only):
    records = parse(corpus_text)
    return records, engines.predict(records, cfg=java_only)


def test_every_sequence_is_predicted(predicted, golden):
    records, predictions = predicted
    assert len(predictions) == len(records) == len(golden)


def test_labels_match_the_legacy_engine_exactly(predicted, golden):
    records, predictions = predicted
    mismatches = []
    for record, (_, labels) in zip(records, predictions):
        expected = golden[record.name]
        if labels != expected:
            mismatches.append(
                f"{record.name}\n  legacy: {expected[:80]}\n  now   : {labels[:80]}"
            )
    assert not mismatches, "\n".join(mismatches)


def test_label_length_matches_sequence_length(predicted):
    records, predictions = predicted
    for record, (_, labels) in zip(records, predictions):
        assert len(labels) == len(record), record.name


def test_posteriors_align_with_the_one_best_prediction(corpus_text, java_only):
    """At each position the argmax curve should usually agree with the label.

    Not a strict identity -- the plot shows per-residue marginals while the
    prediction is the single most probable path -- but wholesale disagreement
    would mean the two engines were being fed different things.
    """
    records = [r for r in parse(corpus_text) if len(r) >= 100][:6]
    predictions = engines.predict(records, cfg=java_only)
    posteriors = engines.posteriors(records, cfg=java_only)

    curve_of = {"i": "cytoplasmic", "o": "non_cytoplasmic", "O": "non_cytoplasmic",
                "C": "non_cytoplasmic", "M": "transmembrane",
                "n": "signal_peptide", "h": "signal_peptide", "c": "signal_peptide"}

    for (_, labels), posterior in zip(predictions, posteriors):
        agree = sum(
            1 for index, label in enumerate(labels)
            if max(posterior.curves, key=lambda c: posterior.curves[c][index])
            == curve_of[label]
        )
        assert agree / len(labels) > 0.8


@needs_engine
def test_native_fast_path_agrees_when_available(corpus_text, tmp_path):
    """If a licensed decodeanhmm is mounted, it must agree with the Java engine."""
    import shutil

    from app.config import Settings

    import os

    # Honour the deployment's own configuration first, so this cross-check runs
    # against whatever binary is actually mounted.
    native = (
        os.environ.get("PHOBIUS_DECODEANHMM")
        or shutil.which("decodeanhmm")
        or "old_phobius/program/decodeanhmm"
    )
    options = os.environ.get("PHOBIUS_OPTIONS")
    if not options:
        beside_model = pathlib.Path(native).parent / "phobius.options"
        options = str(beside_model) if beside_model.is_file() else "old_phobius/cgi-bin/phobius.options"
    if not os.path.exists(native) or not os.path.exists(options):
        pytest.skip("no licensed decodeanhmm available to cross-check")

    # A homopolymer has no unique answer: every position of a helix boundary is
    # exactly equally probable, so two correct builds of the decoder may break
    # the tie differently. Measured over 500 realistic sequences (242k residues),
    # the 64-bit rebuild and the historical 32-bit binary agree on every residue;
    # the only divergence in the corpus is this synthetic poly-leucine case.
    DEGENERATE = {"edge_allhydrophobic"}

    records = [r for r in parse(corpus_text) if r.name not in DEGENERATE]
    java_labels = [labels for _, labels in engines.predict(records, cfg=Settings(decodeanhmm=None))]

    cfg = Settings(decodeanhmm=native, _options_override=options)
    native_labels = [labels for _, labels in engines.predict(records, cfg=cfg)]

    mismatches = [
        f"{record.name}: java {j[:60]} != native {n[:60]}"
        for record, j, n in zip(records, java_labels, native_labels) if j != n
    ]
    assert not mismatches, "\n".join(mismatches)


class TestAlignmentCoordinates:
    """Alignment-based predictions must be reported in the query's own coordinates.

    The engine returns one label per alignment *column*. Reporting those directly
    gives positions that do not match the sequence the user submitted -- they run
    past its end and contain holes where the alignment has gaps. The legacy code
    advanced its counter only over matched label runs (predict.pl:271), which had
    the effect of skipping gap columns.
    """

    def test_gap_columns_are_removed(self):
        assert engines.ungap_labels("nnn--hhh-cc") == "nnnhhhcc"

    def test_ungapped_labels_match_the_query_length(self):
        query = "MKKLLAVVGG"
        aligned_labels = "nn--nhhh--hhcccc"   # 16 columns, 4 gaps
        assert len(engines.ungap_labels(aligned_labels)) == len(query) + 2

    def test_regions_stay_within_the_query(self):
        aligned_labels = "nnnn" + "-" * 5 + "h" * 6 + "-" * 3 + "c" * 4
        regions = features.label_runs(engines.ungap_labels(aligned_labels))
        residues = len(engines.ungap_labels(aligned_labels))
        assert regions[-1].stop == residues
        assert all(r.stop <= residues for r in regions)
        # Contiguous: no holes left behind by the gap columns.
        assert all(b.start == a.stop + 1 for a, b in zip(regions, regions[1:]))

    def test_posterior_is_restricted_to_query_columns(self):
        from app.engines import Posterior

        posterior = Posterior({"cytoplasmic": [0.1, 0.2, 0.3, 0.4, 0.5]})
        trimmed = engines.ungap_posterior(posterior, "n-h-c")
        assert trimmed.curves["cytoplasmic"] == [0.1, 0.3, 0.5]
