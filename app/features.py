"""Turn an engine label string into Phobius features.

The engines emit one label per residue, e.g. ``nnnnhhhhhcccCoooMMMMiii``:

===== =========================================================
``n`` signal peptide n-region
``h`` signal peptide h-region
``c`` signal peptide c-region
``C`` the first residue of the mature protein (cleavage marker)
``i`` cytoplasmic loop
``o`` non-cytoplasmic loop
``O`` non-cytoplasmic loop (alternative state used by some modes)
``M`` transmembrane helix
===== =========================================================

The rendering below is a faithful port of ``printLongN1`` (predict.pl:308-363)
and ``printShortN1`` (predict.pl:365-399). Quirks of the original are preserved
deliberately and flagged in comments -- the regression tests in
``tests/test_golden.py`` compare this output byte-for-byte against the legacy
implementation, so "fixing" one of them will fail the suite on purpose.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

#: Labels that mark part of a signal peptide.
SIGNAL_LABELS = frozenset("nhc")

#: Labels treated as non-cytoplasmic.
NON_CYTOPLASMIC = frozenset("oO")

_RUN = re.compile(r"(\w)\1*")


@dataclass(frozen=True)
class Region:
    """A run of identical labels, with 1-based inclusive coordinates."""

    label: str
    start: int
    stop: int

    def __len__(self) -> int:
        return self.stop - self.start + 1


def label_runs(labels: str) -> list[Region]:
    """Collapse a per-residue label string into runs.

    >>> label_runs("nnhhMMMii")
    [Region(label='n', start=1, stop=2), Region(label='h', start=3, stop=4), \
Region(label='M', start=5, stop=7), Region(label='i', start=8, stop=9)]
    """
    return [
        Region(m.group(1), m.start() + 1, m.end())
        for m in _RUN.finditer(labels.strip())
    ]


def count_transmembrane(regions: list[Region]) -> int:
    """Number of predicted transmembrane helices."""
    return sum(1 for r in regions if r.label == "M")


def has_signal_peptide(regions: list[Region]) -> bool:
    """Whether a signal peptide was predicted.

    The legacy code set its ``$sp`` flag on the h-region specifically
    (predict.pl:328), not on any signal label, so we do the same.
    """
    return any(r.label == "h" for r in regions)


def cleavage_site(regions: list[Region]) -> int | None:
    """Position of the last residue of the signal peptide, if any."""
    for r in regions:
        if r.label == "c":
            return r.stop
    return None


# The UniProt-style feature-table line used by the long output format.
# Trailing whitespace on description-less lines is part of the legacy format.
_FT = "FT   {key:<8} {start:>6} {stop:>6}       {description}"


def long_format(name: str, regions: list[Region]) -> str:
    """Render the long (UniProt feature table) output.

    Port of ``printLongN1``. The SIGNAL line is *prepended* when the c-region is
    reached, which is why it appears above the N/H/C-REGION lines it summarises.
    """
    lines: list[str] = []
    hanger = 0

    for region in regions:
        label, start, stop = region.label, region.start, region.stop

        if label == "n":
            lines.append(_FT.format(key="REGION", start=start, stop=stop,
                                    description="N-REGION."))
        elif label == "h":
            lines.append(_FT.format(key="REGION", start=start, stop=stop,
                                    description="H-REGION."))
        elif label == "c":
            lines.append(_FT.format(key="REGION", start=start, stop=stop,
                                    description="C-REGION."))
            # Prepended, exactly as predict.pl:332 does.
            lines.insert(0, _FT.format(key="SIGNAL", start=1, stop=stop,
                                       description=""))
        elif label == "C":
            # The cleavage marker is folded into the following loop region.
            hanger = start
        elif label in NON_CYTOPLASMIC:
            if hanger > 0:
                start, hanger = hanger, 0
            # Legacy quirk (predict.pl:343): when the whole prediction consists
            # of exactly two regions, a non-cytoplasmic loop is reported as
            # CYTOPLASMIC. Preserved for byte-compatibility.
            description = "CYTOPLASMIC." if len(regions) == 2 else "NON CYTOPLASMIC."
            lines.append(_FT.format(key="TOPO_DOM", start=start, stop=stop,
                                    description=description))
        elif label == "M":
            lines.append(_FT.format(key="TRANSMEM", start=start, stop=stop,
                                    description=""))
        elif label == "i":
            lines.append(_FT.format(key="TOPO_DOM", start=start, stop=stop,
                                    description="CYTOPLASMIC."))

    body = "".join(line + "\n" for line in lines)
    return f"ID   {name}\n{body}//\n"


#: Header for the short output. The misspelling is in the legacy output and is
#: kept so that existing parsers of this format keep working.
SHORT_HEADER = f"{'SEQENCE ID':<30} {'TM':>2} {'SP':>2} PREDICTION"


def topology_string(regions: list[Region]) -> str:
    """The compact topology string, e.g. ``n4-19c24/25o219-238i``."""
    parts: list[str] = []
    for region in regions:
        label, start, stop = region.label, region.start, region.stop
        if label in ("n", "o", "i"):
            parts.append(label)
        elif label == "h":
            parts.append(f"{start}-{stop}")
        elif label == "c":
            parts.append(f"c{stop}/")
        elif label == "C":
            parts.append(str(start))
        elif label == "O":
            parts.append("o")
        elif label == "M":
            parts.append(f"{start}-{stop}")

    topology = "".join(parts)
    # A wholly non-cytoplasmic prediction is reported as cytoplasmic
    # (predict.pl:397); see the note in instructions about not over-reading
    # location for sequences with no helices.
    return "i" if topology == "o" else topology


def short_format(name: str, regions: list[Region]) -> str:
    """Render one line of the short output. Port of ``printShortN1``."""
    tm = count_transmembrane(regions)
    sp = "Y" if has_signal_peptide(regions) else "0"
    return f"{name:<30} {tm:>2} {sp:>2} {topology_string(regions)}"


def to_dict(name: str, sequence: str, regions: list[Region]) -> dict:
    """Structured representation for the JSON API."""
    return {
        "id": name,
        "length": len(sequence),
        "transmembrane_count": count_transmembrane(regions),
        "signal_peptide": has_signal_peptide(regions),
        "cleavage_site": cleavage_site(regions),
        "topology": topology_string(regions),
        "features": [
            {"label": r.label, "start": r.start, "stop": r.stop}
            for r in regions
        ],
    }
