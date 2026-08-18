"""FASTA parsing and residue normalisation.

This module reproduces the normalisation that ``read_fasta()`` in the legacy
``predict.pl`` performed.  That behaviour is load-bearing, not cosmetic: the
prediction engines reject lowercase input outright, so the uppercasing below is
what made the legacy server work at all.  Any change here changes predictions.

Legacy equivalent (predict.pl:200-202)::

    s/[^a-zA-Z-]//g;       # everything but letters and '-' is deleted
    s/[.]/-/g;             # dead code: '.' was already deleted by the line above
    tr/a-z/A-Z/;           # uppercase

Note the second rule never fires, in the original or here: '.' does not match
``[a-zA-Z-]`` so it is deleted before the replacement runs. This matters for
alignments that use '.' as a gap character -- those residues silently vanish
and the alignment loses its column structure. We keep the legacy order so that
predictions stay byte-identical, and rely on the equal-length check in
:func:`parse_alignment` to turn that silent corruption into a loud error.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

#: Characters deleted from sequence lines. Mirrors the legacy ``$ignore`` regex.
_NON_RESIDUE = re.compile(r"[^a-zA-Z-]")

#: Lines starting with these are ignored entirely (legacy ``$comments``).
_COMMENT = re.compile(r"^[#;%]")


class FastaError(ValueError):
    """Raised when input cannot be read as FASTA."""


@dataclass(frozen=True)
class Record:
    """A single normalised FASTA record."""

    header: str
    """The description line with the leading '>' removed, verbatim."""

    sequence: str
    """Residues after normalisation. Uppercase, may contain '-' for alignments."""

    @property
    def name(self) -> str:
        """First whitespace-delimited token of the header.

        This is the identifier the legacy code used (``entry_name``,
        predict.pl:234-239). It is untrusted user input and must never be
        interpolated into markup, shell commands or plot scripts without
        escaping -- that was the source of the injection bugs in predict.pl.
        """
        parts = self.header.split()
        return parts[0] if parts else "UNNAMED"

    def __len__(self) -> int:
        return len(self.sequence)


def normalise(line: str) -> str:
    """Normalise one sequence line exactly as the legacy reader did."""
    # Order is deliberate and matches the legacy reader: the substitution
    # below deletes '.' before the (consequently dead) '.' -> '-' mapping.
    # See the module docstring before changing this.
    return _NON_RESIDUE.sub("", line).replace(".", "-").upper()


def parse(text: str) -> list[Record]:
    """Parse FASTA text into normalised records.

    Bare sequence input (no '>' line at all) is accepted and given the header
    ``UNNAMED``, matching the legacy behaviour at predict.pl:678.
    """
    text = text.replace("\r\n", "\n").replace("\r", "\n")

    if text.strip() and ">" not in text:
        text = ">UNNAMED\n" + text

    records: list[Record] = []
    header: str | None = None
    chunks: list[str] = []

    for line in text.split("\n"):
        if _COMMENT.match(line):
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(Record(header, "".join(chunks)))
            header = line[1:].strip()
            chunks = []
        elif header is not None:
            chunks.append(normalise(line))

    if header is not None:
        records.append(Record(header, "".join(chunks)))

    return [r for r in records if r.sequence]


def format_fasta(records: list[Record], width: int = 60) -> str:
    """Render records back to FASTA, for feeding to a subprocess."""
    out: list[str] = []
    for record in records:
        out.append(f">{record.header}")
        seq = record.sequence
        out.extend(seq[i : i + width] for i in range(0, len(seq), width))
    return "\n".join(out) + "\n"


def parse_alignment(text: str) -> list[Record]:
    """Parse an aligned FASTA block, requiring equal lengths.

    PolyPhobius predicts for the *first* sequence in the alignment; the rest
    supply the homology signal.
    """
    records = parse(text)
    if not records:
        raise FastaError("No sequences found in the alignment.")
    if len(records) < 2:
        raise FastaError(
            "An alignment needs at least two sequences: the query plus one or "
            "more homologues."
        )

    lengths = {len(r) for r in records}
    if len(lengths) != 1:
        first = records[0]
        odd = next(r for r in records if len(r) != len(first))
        raise FastaError(
            f"Aligned sequences must all be the same length, but "
            f"'{first.name}' is {len(first)} and '{odd.name}' is {len(odd)}."
        )
    return records
