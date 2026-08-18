"""Capture reference output from the *legacy* implementation.

Run this once against a checkout that still has the old engine, then commit the
result. ``test_golden.py`` asserts the new service reproduces it byte-for-byte,
which is what makes it safe to have replaced predict.pl at all.

    python tests/generate_golden.py --decodeanhmm old_phobius/program/decodeanhmm \
        --options old_phobius/cgi-bin/phobius.options

The labels are produced by ``decodeanhmm``, the licensed native engine that the
container no longer ships. Recording them here means the substitution stays
checkable without redistributing the binary.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from app.engines import parse_labeled_fasta  # noqa: E402
from app.fasta import format_fasta, parse  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--decodeanhmm", required=True)
    ap.add_argument("--model", default=str(HERE.parent / "engine" / "phobius.model"))
    ap.add_argument("--options", default=None)
    ap.add_argument("--corpus", default=str(HERE / "data" / "regression.fa"))
    args = ap.parse_args()

    records = parse(pathlib.Path(args.corpus).read_text())
    cmd = [args.decodeanhmm]
    if args.options:
        cmd += ["-f", args.options]
    cmd.append(args.model)

    proc = subprocess.run(cmd, input=format_fasta(records), capture_output=True,
                          text=True, check=True)
    predictions = parse_labeled_fasta(proc.stdout)
    if len(predictions) != len(records):
        raise SystemExit(
            f"legacy engine returned {len(predictions)} predictions for "
            f"{len(records)} sequences"
        )

    dest = HERE / "golden" / "labels.tsv"
    dest.write_text(
        "".join(f"{r.name}\t{labels}\n" for r, (_, labels) in zip(records, predictions))
    )
    print(f"wrote {dest} ({len(predictions)} sequences)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
