"""Wrappers around the Phobius prediction engines.

Two things differ structurally from the legacy CGI:

1. **Batching.** ``predict.pl`` spawned one engine process *per sequence*. Here
   a whole request goes through the engine in a single invocation, so a 100
   sequence submission costs two JVM starts rather than two hundred.
2. **No shell.** Every command is built as an argument list and run without a
   shell, and no user-supplied text is ever interpolated into a command,
   filename or plot script.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import Settings, native_engine_status, settings as default_settings
from .fasta import Record, format_fasta

#: Aggregation of per-label posteriors into the four plotted curves.
#:
#: 'C' (first residue of the mature protein) is grouped with the non-cytoplasmic
#: curve. The legacy code was inconsistent here: its PolyPhobius path did the
#: same, but its plain path counted 'C' toward the signal-peptide curve
#: (predict.pl:424 vs :434). We use the PolyPhobius grouping everywhere, which
#: affects a single residue of the plotted curve and nothing in the prediction.
CURVES: dict[str, tuple[str, ...]] = {
    "cytoplasmic": ("i",),
    "non_cytoplasmic": ("o", "O", "C"),
    "transmembrane": ("M",),
    "signal_peptide": ("n", "h", "c"),
}


#: Character the engines emit for an alignment column where the query has a gap.
GAP = "-"


def ungap_labels(labels: str) -> str:
    """Drop alignment gap columns, giving labels in query-sequence coordinates.

    An alignment-based prediction returns one label per alignment *column*, so
    positions derived from it directly would not match the sequence the user
    submitted. The legacy code advanced its position counter only across matched
    label runs (predict.pl:271), which had the same effect.
    """
    return labels.replace(GAP, "")


def ungap_posterior(posterior: "Posterior", labels: str) -> "Posterior":
    """Restrict posterior curves to the columns where the query has a residue."""
    keep = [i for i, label in enumerate(labels) if label != GAP]
    return Posterior({
        curve: [values[i] for i in keep if i < len(values)]
        for curve, values in posterior.curves.items()
    })


class EngineError(RuntimeError):
    """An engine failed, timed out, or produced output we cannot parse."""


@dataclass
class Posterior:
    """Posterior label probabilities for one sequence."""

    curves: dict[str, list[float]]

    def __len__(self) -> int:
        return len(next(iter(self.curves.values()))) if self.curves else 0


def _run(cmd: list[str], stdin: str, timeout: int, what: str) -> str:
    """Run a command with no shell, returning stdout.

    The engines print banners and harmless numerical diagnostics to stderr, so
    stderr is only surfaced when the command actually fails.
    """
    try:
        proc = subprocess.run(
            cmd,
            input=stdin,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise EngineError(
            f"{what} did not finish within {timeout}s. Try submitting fewer or "
            f"shorter sequences."
        ) from exc
    except FileNotFoundError as exc:
        raise EngineError(f"{what}: executable not found ({cmd[0]}).") from exc

    if proc.returncode != 0:
        detail = (proc.stderr or "").strip().splitlines()
        tail = detail[-1] if detail else f"exit status {proc.returncode}"
        raise EngineError(f"{what} failed: {tail}")
    return proc.stdout


# --------------------------------------------------------------------------
# output parsing
# --------------------------------------------------------------------------

def parse_labeled_fasta(stdout: str) -> list[tuple[str, str]]:
    """Parse ``>name`` / ``?X labels`` engine output into (name, labels) pairs.

    Label lines are recognised by a '?' followed by a single non-space tag and a
    space, as in the legacy reader (predict.pl:267).
    """
    results: list[tuple[str, str]] = []
    name: str | None = None
    labels: list[str] = []

    for line in stdout.splitlines():
        if line.startswith(">"):
            if name is not None:
                results.append((name, "".join(labels)))
            name = line[1:].strip()
            labels = []
        elif len(line) > 2 and line[0] == "?" and not line[1].isspace() and line[2] == " ":
            labels.append(line[2:].replace(" ", ""))

    if name is not None:
        results.append((name, "".join(labels)))
    return results


def parse_posteriors(stdout: str) -> list[Posterior]:
    """Parse ``Run -plp`` output into one :class:`Posterior` per sequence.

    A '#' line introduces each sequence and names the columns. The column order
    is *not* stable between invocations -- it genuinely varies -- so the header
    must be parsed rather than assumed. The legacy code had the same discovery
    baked into ``setPlpPos`` (predict.pl:438).

    The blocks carry no sequence identifier, so callers must match them to the
    input by position and verify the lengths.
    """
    blocks: list[Posterior] = []
    order: list[str] = []
    current: dict[str, list[float]] | None = None

    def flush() -> None:
        if current is not None:
            blocks.append(Posterior({
                curve: [
                    sum(current[label][i] for label in labels if label in current)
                    for i in range(len(next(iter(current.values()))))
                ]
                for curve, labels in CURVES.items()
            } if current and next(iter(current.values()), None) is not None else {}))

    for line in stdout.splitlines():
        if line.startswith("#"):
            flush()
            order = line[1:].split()
            current = {label: [] for label in order}
            continue
        if current is None or not line.strip():
            continue
        fields = line.split()
        # First field is the row number; the rest are probabilities in header order.
        if len(fields) != len(order) + 1:
            continue
        for label, value in zip(order, fields[1:]):
            current[label].append(float(value))

    flush()
    return blocks


# --------------------------------------------------------------------------
# engine invocations
# --------------------------------------------------------------------------

def _java(cls: str, args: list[str], cfg: Settings) -> list[str]:
    return [cfg.java, *cfg.java_opts, "-cp", cfg.classpath, cls, *args]


def predict(
    records: list[Record],
    constraints: list[str] | None = None,
    cfg: Settings | None = None,
) -> list[tuple[str, str]]:
    """One-best prediction for each record. Returns (name, labels) pairs.

    ``constraints`` is a list of already-validated tokens such as
    ``["M_40-60", "o_150", "O_150"]``; see :func:`app.models.constraint_tokens`.
    """
    cfg = cfg or default_settings

    if not constraints and native_engine_status(cfg)[0]:
        return _predict_decodeanhmm(records, cfg)

    args = ["-raw"]
    if constraints:
        args += ["-c", *constraints]
    stdout = _run(
        _java("se.ki.cgb.hmmdecode.Phobius", args, cfg),
        format_fasta(records),
        cfg.engine_timeout,
        "Phobius",
    )
    return _validate(parse_labeled_fasta(stdout), records)


def _predict_decodeanhmm(records: list[Record], cfg: Settings) -> list[tuple[str, str]]:
    """Optional fast path using the licensed native binary, if one is mounted.

    Verified to produce byte-identical labels to the Java engine; see
    ``tests/test_golden.py``.
    """
    cmd = [str(cfg.decodeanhmm)]
    if cfg.phobius_options.is_file():
        cmd += ["-f", str(cfg.phobius_options)]
    cmd.append(str(cfg.model))
    stdout = _run(cmd, format_fasta(records), cfg.engine_timeout, "decodeanhmm")
    return _validate(parse_labeled_fasta(stdout), records)


def predict_alignment(
    records: list[Record],
    cfg: Settings | None = None,
) -> tuple[str, str]:
    """PolyPhobius prediction from an alignment. Predicts for the first sequence."""
    cfg = cfg or default_settings
    with _alignment_file(records, cfg) as path:
        stdout = _run(
            _java("se.ki.cgb.hmmdecode.Run", ["-a", str(cfg.model), str(path)], cfg),
            "",
            cfg.engine_timeout,
            "PolyPhobius",
        )
    parsed = parse_labeled_fasta(stdout)
    if not parsed:
        raise EngineError("PolyPhobius returned no prediction.")
    name, labels = parsed[0]
    expected = len(records[0])
    if len(labels) != expected:
        raise EngineError(
            f"PolyPhobius returned {len(labels)} labels for a {expected}-column "
            f"alignment."
        )
    return name, labels


def posteriors(
    records: list[Record],
    constraints: list[str] | None = None,
    cfg: Settings | None = None,
) -> list[Posterior]:
    """Posterior label probabilities for each record."""
    cfg = cfg or default_settings
    args = ["-plp"]
    if constraints:
        args += ["-c", *constraints]
    args.append(str(cfg.model))
    stdout = _run(
        _java("se.ki.cgb.hmmdecode.Run", args, cfg),
        format_fasta(records),
        cfg.engine_timeout,
        "Phobius (posterior probabilities)",
    )
    blocks = parse_posteriors(stdout)
    if len(blocks) != len(records):
        raise EngineError(
            f"Expected posteriors for {len(records)} sequences, got {len(blocks)}."
        )
    for block, record in zip(blocks, records):
        if len(block) != len(record):
            raise EngineError(
                f"Posterior length {len(block)} does not match sequence "
                f"'{record.name}' of length {len(record)}."
            )
    return blocks


def posteriors_alignment(
    records: list[Record],
    cfg: Settings | None = None,
) -> Posterior:
    """Posterior label probabilities from an alignment."""
    cfg = cfg or default_settings
    with _alignment_file(records, cfg) as path:
        stdout = _run(
            _java("se.ki.cgb.hmmdecode.Run", ["-a", "-plp", str(cfg.model), str(path)], cfg),
            "",
            cfg.engine_timeout,
            "PolyPhobius (posterior probabilities)",
        )
    blocks = parse_posteriors(stdout)
    if not blocks:
        raise EngineError("PolyPhobius returned no posterior probabilities.")
    return blocks[0]


class _alignment_file:
    """Write an alignment to a private temporary file for the engine to read.

    ``Run -a`` only accepts a path, not stdin. The file is created with a random
    name in a private directory and removed afterwards -- unlike the legacy
    server, which wrote alignments to predictable names under the web root where
    anyone could fetch other users' sequences.
    """

    def __init__(self, records: list[Record], cfg: Settings) -> None:
        self._records = records
        self._cfg = cfg

    def __enter__(self) -> Path:
        import tempfile

        self._dir = tempfile.TemporaryDirectory(prefix="phobius-")
        path = Path(self._dir.name) / "alignment.fa"
        path.write_text(format_fasta(self._records))
        return path

    def __exit__(self, *exc: object) -> None:
        self._dir.cleanup()


def _validate(parsed: list[tuple[str, str]], records: list[Record]) -> list[tuple[str, str]]:
    if len(parsed) != len(records):
        raise EngineError(
            f"Expected predictions for {len(records)} sequences, got {len(parsed)}."
        )
    for (name, labels), record in zip(parsed, records):
        if len(labels) != len(record):
            raise EngineError(
                f"Prediction for '{record.name}' has {len(labels)} labels but the "
                f"sequence has {len(record)} residues."
            )
    return parsed
