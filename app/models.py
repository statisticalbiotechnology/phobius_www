"""Validated request models.

Everything the user can influence is parsed and bounded here, before it reaches
an engine. In particular :func:`build_constraints` replaces the legacy
``getConstraint`` (predict.pl:561-607), which validated positions but then built
a command fragment by string concatenation.
"""

from __future__ import annotations

import re
from enum import Enum

from pydantic import BaseModel, Field, field_validator

from .config import settings
from .fasta import Record

#: A single constraint token: a position, or an inclusive range.
_TOKEN = re.compile(r"^(\d+)(?:-(\d+))?$")

#: Location codes accepted from the constrained-prediction form.
LOCATION_NAMES: dict[str, str] = {
    "M": "Membrane",
    "i": "Cytoplasmic loop",
    "o": "Non-cytoplasmic loop",
    "n": "Signal peptide",
}


class OutputFormat(str, Enum):
    """Output selector. Values match the legacy form for URL compatibility."""

    SHORT = "short"
    LONG = "nog"
    LONG_WITH_PLOT = "plp"
    PLOT_ONLY = "aplp"


class InputFormat(str, Enum):
    SEARCH = "blast"
    ALIGNMENT = "align"


class ConstraintError(ValueError):
    """A constraint could not be understood or does not fit the sequence."""


class Constraint(BaseModel):
    """One user-specified constraint on a residue or range."""

    location: str = Field(pattern="^[Mion]$")
    start: int = Field(ge=1)
    stop: int = Field(ge=1)

    @property
    def description(self) -> str:
        where = f"{self.start}-{self.stop}" if self.stop != self.start else str(self.start)
        return f"Amino acid {where} in {LOCATION_NAMES[self.location]}"

    def to_argument(self) -> list[str]:
        """Engine tokens for this constraint.

        A non-cytoplasmic constraint emits both 'o' and 'O' because the model
        has two non-cytoplasmic states; the legacy code did the same at
        predict.pl:589.
        """
        span = f"{self.start}-{self.stop}" if self.stop != self.start else str(self.start)
        tokens = [f"{self.location}_{span}"]
        if self.location == "o":
            tokens.append(f"O_{span}")
        return tokens


def parse_constraint_field(location: str, text: str, sequence_length: int) -> list[Constraint]:
    """Parse one constrained-prediction input box.

    Accepts space-separated positions and ranges. 'C' stands for the C-terminus
    and is replaced by the sequence length, as documented on the form.
    """
    if not text or not text.strip():
        return []

    normalised = re.sub(r"[Cc]", str(sequence_length), text.strip())
    constraints: list[Constraint] = []

    for token in normalised.split():
        match = _TOKEN.match(token)
        if not match:
            raise ConstraintError(
                f"Could not read '{token}' as a position or range. Use numbers "
                f"like '33' or '33-40', separated by spaces."
            )
        start = int(match.group(1))
        stop = int(match.group(2)) if match.group(2) else start
        if start > stop:
            start, stop = stop, start
        for value in (start, stop):
            if value < 1:
                raise ConstraintError(f"Position {value} is less than 1.")
            if value > sequence_length:
                raise ConstraintError(
                    f"Position {value} is beyond the end of the sequence "
                    f"({sequence_length} residues)."
                )
        constraints.append(Constraint(location=location, start=start, stop=stop))

    return constraints


def build_constraints(
    sequence_length: int,
    membrane: str = "",
    cytoplasmic: str = "",
    non_cytoplasmic: str = "",
    signal_peptide: bool = False,
) -> list[Constraint]:
    """Collect all constraints from the constrained-prediction form."""
    constraints: list[Constraint] = []
    constraints += parse_constraint_field("M", membrane, sequence_length)
    constraints += parse_constraint_field("o", non_cytoplasmic, sequence_length)
    constraints += parse_constraint_field("i", cytoplasmic, sequence_length)
    if signal_peptide:
        constraints.append(Constraint(location="n", start=1, stop=1))
    return constraints


def constraint_tokens(constraints: list[Constraint]) -> list[str]:
    """Render constraints as separate argv tokens.

    Each token must be its own argument. The legacy code got that for free from
    Perl's ``split(/ +/, ...)`` (predict.pl:255); passing them joined makes the
    engine treat the whole string as an input filename and crash.

    Every token is machine-generated from validated integers and a location code
    matching ``^[Mion]$``, so nothing user-supplied reaches argv verbatim.
    """
    tokens: list[str] = []
    for constraint in constraints:
        tokens.extend(constraint.to_argument())
    return tokens


class SubmissionError(ValueError):
    """The submission is empty or exceeds a configured limit."""


def check_limits(records: list[Record]) -> None:
    """Reject submissions that would tie up a worker for too long.

    The legacy server capped only the pasted textarea (10 000 characters) and
    left uploads unbounded, while serialising every request behind one global
    lock -- so a single large upload stalled the whole site.
    """
    if not records:
        raise SubmissionError(
            "No sequences found. Paste a protein sequence in FASTA format, or "
            "choose a file to upload."
        )
    if len(records) > settings.max_sequences:
        raise SubmissionError(
            f"{len(records)} sequences submitted, but this server accepts at most "
            f"{settings.max_sequences} per request. For larger jobs please use the "
            f"standalone package."
        )

    total = sum(len(r) for r in records)
    if total > settings.max_residues_total:
        raise SubmissionError(
            f"{total:,} residues submitted, but this server accepts at most "
            f"{settings.max_residues_total:,} per request. For larger jobs please "
            f"use the standalone package."
        )

    for record in records:
        if len(record) > settings.max_residues_single:
            raise SubmissionError(
                f"Sequence '{record.name}' is {len(record):,} residues; the limit "
                f"is {settings.max_residues_single:,}."
            )


class ApiRequest(BaseModel):
    """JSON request body for ``POST /api/predict``."""

    sequence: str = Field(description="One or more protein sequences in FASTA format.")
    plot: bool = Field(default=False, description="Include an SVG posterior plot per sequence.")

    @field_validator("sequence")
    @classmethod
    def _not_empty(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("sequence must not be empty")
        return value
