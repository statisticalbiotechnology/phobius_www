"""Render posterior label probabilities as a self-contained SVG.

This replaces the legacy gnuplot pipeline, which wrote a generated script and a
data file to a predictable path under the web root and then shelled out. That
design leaked every user's submitted sequences to anyone who guessed a filename,
and interpolated the unescaped FASTA header into the script -- where gnuplot's
backquote substitution turns it into a shell command.

Here the plot is built as a string and returned; nothing touches the filesystem,
and every piece of user-controlled text goes through :func:`xml.sax.saxutils.escape`.

SVG rather than a raster plot keeps the image sharp at any zoom, keeps the
container free of a plotting stack and its font configuration, and makes the
output diffable in tests.
"""

from __future__ import annotations

from xml.sax.saxutils import escape, quoteattr

from .engines import Posterior
from .features import Region

# Okabe-Ito palette: distinguishable under the common forms of colour vision
# deficiency. The mapping of curve to role matches the legacy gnuplot line
# types (transmembrane/cytoplasmic/non-cytoplasmic/signal peptide).
COLOURS: dict[str, str] = {
    "transmembrane": "#D55E00",
    "cytoplasmic": "#009E73",
    "non_cytoplasmic": "#0072B2",
    "signal_peptide": "#CC79A7",
}

LABELS: dict[str, str] = {
    "transmembrane": "transmembrane",
    "cytoplasmic": "cytoplasmic",
    "non_cytoplasmic": "non-cytoplasmic",
    "signal_peptide": "signal peptide",
}

#: Colour of each region type in the prediction track under the axis.
TRACK_COLOURS: dict[str, str] = {
    "M": COLOURS["transmembrane"],
    "i": COLOURS["cytoplasmic"],
    "o": COLOURS["non_cytoplasmic"],
    "O": COLOURS["non_cytoplasmic"],
    "C": COLOURS["non_cytoplasmic"],
    "n": COLOURS["signal_peptide"],
    "h": COLOURS["signal_peptide"],
    "c": COLOURS["signal_peptide"],
}

_W, _H = 900, 380
_LEFT, _RIGHT, _TOP, _BOTTOM = 62, 18, 40, 74
_TRACK_H = 12
_TRACK_GAP = 8


def _nice_ticks(limit: int, target: int = 8) -> list[int]:
    """Round tick positions covering 1..limit."""
    if limit <= 1:
        return [1]
    raw = limit / target
    magnitude = 10 ** max(0, len(str(int(raw))) - 1)
    for factor in (1, 2, 2.5, 5, 10):
        step = int(magnitude * factor)
        if step and limit / step <= target:
            break
    else:  # pragma: no cover - defensive
        step = max(1, limit // target)
    ticks = list(range(step, limit + 1, step))
    return [1] + [t for t in ticks if t > step // 2]


def render(
    posterior: Posterior,
    regions: list[Region],
    title: str,
    *,
    subtitle: str = "",
) -> str:
    """Return an SVG document for one sequence's posterior probabilities."""
    length = len(posterior)
    if length == 0:
        return ""

    plot_w = _W - _LEFT - _RIGHT
    plot_h = _H - _TOP - _BOTTOM

    def x(pos: float) -> float:
        """Map residue position (1-based) to an x coordinate."""
        if length == 1:
            return _LEFT + plot_w / 2
        return _LEFT + (pos - 1) / (length - 1) * plot_w

    def y(prob: float) -> float:
        return _TOP + (1 - prob) * plot_h

    out: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {_W} {_H}" '
        f'width="100%" role="img" aria-label={quoteattr("Posterior label probabilities for " + title)} '
        f'font-family="system-ui, -apple-system, Segoe UI, Helvetica, Arial, sans-serif">',
        f'<rect width="{_W}" height="{_H}" fill="#ffffff"/>',
    ]

    # --- title -------------------------------------------------------------
    out.append(
        f'<text x="{_W / 2}" y="20" text-anchor="middle" font-size="14" '
        f'font-weight="600" fill="#111">{escape(title)}</text>'
    )
    if subtitle:
        out.append(
            f'<text x="{_W / 2}" y="34" text-anchor="middle" font-size="11" '
            f'fill="#555">{escape(subtitle)}</text>'
        )

    # --- grid and axes -----------------------------------------------------
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        gy = y(frac)
        out.append(
            f'<line x1="{_LEFT}" y1="{gy:.1f}" x2="{_LEFT + plot_w}" y2="{gy:.1f}" '
            f'stroke="#e6e6e6" stroke-width="1"/>'
        )
        out.append(
            f'<text x="{_LEFT - 8}" y="{gy + 4:.1f}" text-anchor="end" font-size="11" '
            f'fill="#555">{frac:g}</text>'
        )

    for tick in _nice_ticks(length):
        tx = x(tick)
        out.append(
            f'<line x1="{tx:.1f}" y1="{_TOP + plot_h}" x2="{tx:.1f}" '
            f'y2="{_TOP + plot_h + 4}" stroke="#999" stroke-width="1"/>'
        )
        out.append(
            f'<text x="{tx:.1f}" y="{_TOP + plot_h + 16}" text-anchor="middle" '
            f'font-size="11" fill="#555">{tick}</text>'
        )

    out.append(
        f'<text transform="translate(16,{_TOP + plot_h / 2}) rotate(-90)" '
        f'text-anchor="middle" font-size="11" fill="#333">'
        f'Posterior label probability</text>'
    )

    # --- curves ------------------------------------------------------------
    # Transmembrane is drawn filled, matching the legacy plot's use of impulses
    # for that series; it reads as a solid block where a helix is confident.
    tm = posterior.curves.get("transmembrane", [])
    if tm:
        pts = " ".join(f"{x(i + 1):.1f},{y(v):.1f}" for i, v in enumerate(tm))
        out.append(
            f'<polygon points="{x(1):.1f},{y(0):.1f} {pts} {x(length):.1f},{y(0):.1f}" '
            f'fill="{COLOURS["transmembrane"]}" fill-opacity="0.28"/>'
        )

    for curve in ("transmembrane", "cytoplasmic", "non_cytoplasmic", "signal_peptide"):
        values = posterior.curves.get(curve)
        if not values:
            continue
        pts = " ".join(f"{x(i + 1):.1f},{y(v):.1f}" for i, v in enumerate(values))
        out.append(
            f'<polyline points="{pts}" fill="none" stroke="{COLOURS[curve]}" '
            f'stroke-width="1.8" stroke-linejoin="round"/>'
        )

    out.append(
        f'<rect x="{_LEFT}" y="{_TOP}" width="{plot_w}" height="{plot_h}" '
        f'fill="none" stroke="#999" stroke-width="1"/>'
    )

    # --- prediction track --------------------------------------------------
    # The legacy plot drew the one-best prediction as coloured bars just below
    # the axis (documented in instructions.html as "between -0.04 and 0").
    track_y = _TOP + plot_h + _TRACK_GAP + 14
    for index, region in enumerate(regions):
        colour = TRACK_COLOURS.get(region.label)
        if colour is None:
            continue
        # Each block runs to where the next one starts, so adjacent regions abut
        # exactly instead of leaving sub-pixel gaps along the track.
        x0 = x(region.start)
        x1 = x(regions[index + 1].start) if index + 1 < len(regions) else x(length)
        out.append(
            f'<rect x="{x0:.1f}" y="{track_y}" width="{max(x1 - x0, 1):.1f}" '
            f'height="{_TRACK_H}" fill="{colour}">'
            f'<title>{escape(region.label)} {region.start}-{region.stop}</title></rect>'
        )
    out.append(
        f'<text x="{_LEFT - 8}" y="{track_y + _TRACK_H - 2}" text-anchor="end" '
        f'font-size="10" fill="#555">prediction</text>'
    )

    # --- legend ------------------------------------------------------------
    legend_y = _H - 14
    span = plot_w / len(COLOURS)
    for index, curve in enumerate(
        ("transmembrane", "cytoplasmic", "non_cytoplasmic", "signal_peptide")
    ):
        lx = _LEFT + index * span
        out.append(
            f'<line x1="{lx:.1f}" y1="{legend_y - 4}" x2="{lx + 18:.1f}" '
            f'y2="{legend_y - 4}" stroke="{COLOURS[curve]}" stroke-width="3"/>'
        )
        out.append(
            f'<text x="{lx + 24:.1f}" y="{legend_y}" font-size="11" fill="#333">'
            f'{escape(LABELS[curve])}</text>'
        )

    out.append("</svg>")
    return "".join(out)
