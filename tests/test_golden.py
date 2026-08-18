"""Regression against the legacy implementation.

``tests/golden/labels.tsv`` holds per-residue labels produced by ``decodeanhmm``
-- the licensed native engine the container no longer ships. These tests assert
that the Java engine used in its place returns exactly the same thing, which is
the whole basis for having dropped the native binary.

Regenerate the golden file only with ``tests/generate_golden.py``, and only
against the legacy engine. If a change makes these fail, the change altered
predictions.
"""

import pytest

from app import engines, features
from app.fasta import parse

from conftest import needs_engine

pytestmark = needs_engine


@pytest.fixture(scope="module")
def predicted(corpus_text):
    records = parse(corpus_text)
    return records, engines.predict(records)


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


def test_posteriors_align_with_the_one_best_prediction(corpus_text):
    """At each position the argmax curve should usually agree with the label.

    Not a strict identity -- the plot shows per-residue marginals while the
    prediction is the single most probable path -- but wholesale disagreement
    would mean the two engines were being fed different things.
    """
    records = [r for r in parse(corpus_text) if len(r) >= 100][:6]
    predictions = engines.predict(records)
    posteriors = engines.posteriors(records)

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

    native = shutil.which("decodeanhmm") or "old_phobius/program/decodeanhmm"
    options = os.environ.get("PHOBIUS_OPTIONS", "old_phobius/cgi-bin/phobius.options")
    if not os.path.exists(native) or not os.path.exists(options):
        pytest.skip("no licensed decodeanhmm available to cross-check")

    records = parse(corpus_text)
    java_labels = [labels for _, labels in engines.predict(records)]

    cfg = Settings(decodeanhmm=native, _options_override=options)
    native_labels = [labels for _, labels in engines.predict(records, cfg=cfg)]

    assert java_labels == native_labels
