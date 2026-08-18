"""Input normalisation. These rules decide what the engine actually sees."""

import pytest

from app.fasta import FastaError, Record, normalise, parse, parse_alignment


def test_lowercase_is_uppercased():
    # Not cosmetic: the native engine rejects lowercase input outright, so this
    # is what made the legacy server work at all.
    assert parse(">a\nmykliv\n")[0].sequence == "MYKLIV"


def test_digits_and_whitespace_are_stripped():
    assert normalise("  MY 12 KL\tIV ") == "MYKLIV"


def test_gap_characters_survive_but_dots_do_not():
    # Documented legacy quirk: '.' is deleted rather than becoming a gap,
    # because the deletion rule runs first. See app/fasta.py.
    assert normalise("AC-DE") == "AC-DE"
    assert normalise("AC.DE") == "ACDE"


def test_comment_lines_are_ignored():
    assert parse(">a\n# comment\n%another\n;third\nMKKL\n")[0].sequence == "MKKL"


def test_bare_sequence_gets_a_name():
    record = parse("MKKL\n")[0]
    assert record.name == "UNNAMED"
    assert record.sequence == "MKKL"


def test_crlf_input():
    assert parse(">a\r\nMKKL\r\nMKKL\r\n")[0].sequence == "MKKLMKKL"


def test_empty_records_are_dropped():
    assert parse(">a\n\n>b\nMKKL\n") == [Record("b", "MKKL")]


def test_name_is_first_token():
    assert parse(">sp|P1|X some description\nMK\n")[0].name == "sp|P1|X"


def test_alignment_requires_equal_lengths():
    with pytest.raises(FastaError, match="same length"):
        parse_alignment(">a\nMKKL\n>b\nMKK\n")


def test_alignment_requires_two_sequences():
    with pytest.raises(FastaError, match="at least two"):
        parse_alignment(">a\nMKKL\n")


def test_alignment_accepts_gaps():
    records = parse_alignment(">a\nMK-KL\n>b\nMKQKL\n")
    assert [len(r) for r in records] == [5, 5]
