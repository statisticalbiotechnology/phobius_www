"""Runtime configuration, all overridable by environment variable.

The licensed model file is deliberately *not* baked into the container image; it
is mounted from private storage at runtime, and the service refuses to start
without it rather than failing on the first request.

Where that mount lands differs by host, and SciLifeLab Serve lets the operator
choose the mount path while offering no way to set an environment variable for a
custom app. So rather than requiring configuration, the usual locations are
probed for ``phobius.model``. Anything explicit still wins:

1. ``PHOBIUS_MODEL`` -- a full path to the file
2. ``PHOBIUS_DATA_DIR`` -- a directory searched before the defaults
3. the directories in :data:`DATA_DIRS`
4. ``engine/phobius.model`` next to the source, for local development

Both variables can be set without platform support by dropping a ``phobius.env``
file on the same mounted storage; ``start-script.sh`` loads it before startup.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent

#: Directories probed for mounted data files, in order. These are the paths
#: hosting platforms commonly mount project storage at.
DATA_DIRS: tuple[str, ...] = ("/mnt/data", "/home/data", "/data", "/app/data")

#: Sub-directories of each mount that are searched as well as the mount itself.
#: A BLAST database is nine files sharing a prefix, so keeping it in a ``db/``
#: sub-directory next to the model is markedly tidier than having it loose.
DATA_SUBDIRS: tuple[str, ...] = ("", "db")


def data_dirs() -> list[Path]:
    """Directories searched for mounted data files, most specific first.

    Each mount point is searched both directly and in its ``db/``
    sub-directory, so files may be organised either way.
    """
    mounts: list[str] = []
    explicit = os.environ.get("PHOBIUS_DATA_DIR")
    if explicit:
        mounts.append(explicit)
    mounts.extend(DATA_DIRS)

    searched: list[Path] = []
    # dict.fromkeys keeps order while removing duplicates.
    for mount in dict.fromkeys(mounts):
        for sub in DATA_SUBDIRS:
            searched.append(Path(mount) / sub if sub else Path(mount))
    return searched


def discover(filename: str, env_var: str, fallback: Path) -> Path:
    """Locate a mounted data file.

    An explicit environment variable always wins; otherwise the candidate
    directories are searched, falling back to a development-local path.
    """
    explicit = os.environ.get(env_var)
    if explicit:
        return Path(explicit)
    for directory in data_dirs():
        candidate = directory / filename
        if candidate.is_file():
            return candidate
    return fallback


def discover_blast_db() -> Path:
    """Locate the BLAST database, returning its *prefix*.

    A BLAST database is a set of files sharing a prefix, so discovery keys off
    an index file: ``.pin`` for a plain database, ``.pal`` for a split one.
    """
    explicit = os.environ.get("PHOBIUS_BLAST_DB")
    if explicit:
        return Path(explicit)
    for directory in data_dirs():
        for marker in ("swissprot.pin", "swissprot.pal"):
            if (directory / marker).is_file():
                return directory / "swissprot"
    return _ROOT / "data" / "swissprot"


def discover_downloads_dir() -> Path:
    """Locate the directory holding downloadable bundles.

    Searched like every other mounted resource, so ``/home/data/download`` works
    with no configuration.
    """
    explicit = os.environ.get("PHOBIUS_DOWNLOADS")
    if explicit:
        return Path(explicit)
    for directory in data_dirs():
        candidate = directory / "download"
        if candidate.is_dir():
            return candidate
    return _ROOT / "download"


def _path(env: str, default: Path) -> Path:
    return Path(os.environ.get(env, str(default)))


def _int(env: str, default: int) -> int:
    return int(os.environ.get(env, default))


@dataclass(frozen=True)
class Settings:
    # --- engine assets -----------------------------------------------------
    model: Path = field(default_factory=lambda: discover(
        "phobius.model", "PHOBIUS_MODEL", _ROOT / "engine" / "phobius.model"))
    engine_dir: Path = field(default_factory=lambda: _path("PHOBIUS_ENGINE_DIR", _ROOT / "engine"))
    java: str = field(default_factory=lambda: os.environ.get("PHOBIUS_JAVA", "java"))

    # --- homology search (PolyPhobius "search for homologues" mode) --------
    # BLAST rather than DIAMOND: for a single query DIAMOND spends its time
    # building a seed index over the whole database, which dwarfs the search.
    # Measured against Swiss-Prot, blastp takes 1.2 s where DIAMOND took 42.8 s,
    # and returned a superset of its hits.
    blastp: str = field(default_factory=lambda: os.environ.get("PHOBIUS_BLASTP", "blastp"))
    blastdbcmd: str = field(default_factory=lambda: os.environ.get("PHOBIUS_BLASTDBCMD", "blastdbcmd"))
    kalign: str = field(default_factory=lambda: os.environ.get("PHOBIUS_KALIGN", "kalign"))
    blast_db: Path = field(default_factory=lambda: discover_blast_db())

    # --- limits ------------------------------------------------------------
    # Measured throughput of the Java engine is ~2.6k residues/s, so the
    # default residue cap corresponds to roughly 20 seconds of compute.
    max_sequences: int = field(default_factory=lambda: _int("PHOBIUS_MAX_SEQUENCES", 100))
    max_residues_total: int = field(default_factory=lambda: _int("PHOBIUS_MAX_RESIDUES", 50_000))
    max_residues_single: int = field(default_factory=lambda: _int("PHOBIUS_MAX_RESIDUES_SINGLE", 10_000))
    max_upload_bytes: int = field(default_factory=lambda: _int("PHOBIUS_MAX_UPLOAD_BYTES", 4_000_000))
    engine_timeout: int = field(default_factory=lambda: _int("PHOBIUS_ENGINE_TIMEOUT", 120))
    homology_timeout: int = field(default_factory=lambda: _int("PHOBIUS_HOMOLOGY_TIMEOUT", 120))
    max_concurrency: int = field(default_factory=lambda: _int("PHOBIUS_MAX_CONCURRENCY", 2))

    # --- downloads ---------------------------------------------------------
    # The standalone Phobius bundle is licensed, so like the model it is served
    # from mounted storage rather than shipped in the image.
    downloads_dir: Path = field(default_factory=lambda: discover_downloads_dir())

    # --- optional fast path ------------------------------------------------
    # decodeanhmm is licensed and must never ship in a public image. If a copy
    # is present on a mounted volume it is ~8x faster for plain predictions;
    # the golden tests assert both engines agree byte-for-byte.
    decodeanhmm: str | None = field(default_factory=lambda: os.environ.get("PHOBIUS_DECODEANHMM"))
    _options_override: str | None = field(default_factory=lambda: os.environ.get("PHOBIUS_OPTIONS"))

    @property
    def phobius_options(self) -> Path:
        """Decoder options for the native engine.

        Defaults to sitting beside the model, because both are part of the
        licensed distribution and are mounted together. Without it decodeanhmm
        omits the label lines entirely and its output cannot be parsed.
        """
        if self._options_override:
            return Path(self._options_override)
        return self.model.parent / "phobius.options"

    @property
    def classpath(self) -> str:
        """Classpath for the Java engines.

        ``se.ki.cgb.hmmdecode.Phobius`` loads ``phobius.model`` as a *classpath
        resource*, not as an argument, so the directory holding the model has to
        be on the classpath. ``Run`` takes the model as an argument instead.
        """
        return os.pathsep.join([
            str(self.engine_dir / "homologhmm.jar"),
            str(self.engine_dir / "biojava.jar"),
            str(self.model.parent),
        ])

    @property
    def java_opts(self) -> list[str]:
        """Flags that cut JVM startup, which dominates single-sequence latency."""
        return ["-XX:TieredStopAtLevel=1", "-XX:+UseSerialGC", "-Xmx512m"]

    def homology_search_available(self) -> bool:
        """Whether the PolyPhobius homology search can actually be offered."""
        return homology_search_status(self)[0]

    def check(self) -> list[str]:
        """Return a list of fatal configuration problems (empty if healthy)."""
        problems: list[str] = []
        if not self.model.is_file():
            searched = "\n      ".join(str(d / "phobius.model") for d in data_dirs())
            problems.append(
                "phobius.model not found. It is licensed and is deliberately not "
                "shipped in this image.\n"
                "    Mount it under any of these paths:\n"
                f"      {searched}\n"
                "    or set PHOBIUS_MODEL / PHOBIUS_DATA_DIR -- which can be done "
                "without platform\n    support by putting a phobius.env file next "
                "to the model on the same storage."
            )
        if self.model.name != "phobius.model":
            problems.append(
                f"The model must be named exactly 'phobius.model' (found "
                f"'{self.model.name}') because the Java engine loads it by that "
                f"resource name from the classpath."
            )
        for jar in ("homologhmm.jar", "biojava.jar"):
            if not (self.engine_dir / jar).is_file():
                problems.append(f"Missing engine jar: {self.engine_dir / jar}")
        if not shutil.which(self.java):
            problems.append(f"Java runtime '{self.java}' not found on PATH.")
        # The native engine is deliberately not checked here. It is an optional
        # accelerator producing identical output, so a broken one degrades to the
        # Java engine rather than taking the service down --
        # see :func:`native_engine_status`.
        return problems


#: A short sequence run through the native engine to prove it works end to end.
_PROBE_SEQUENCE = "MKKLLAVVGGSILAWQPTRDELAAAWWWLLLIIIVVVGGGAAAFFFMMM"

#: Cache for :func:`native_engine_status`, keyed by binary and options path.
_native_status: dict[tuple[str, str], tuple[bool, str]] = {}


def native_engine_status(cfg: "Settings") -> tuple[bool, str]:
    """Whether the optional native engine is configured *and actually runs*.

    Existence is not enough. The historical decodeanhmm is a 32-bit binary and
    the container is 64-bit with no i386 libraries, so it fails at exec time
    with a bare "not found". Checking only for the file would leave the service
    reporting healthy and then failing every prediction, which is worse than not
    using it at all. So it is executed once and, if it does not work, the reason
    is logged and the Java engine is used instead.
    """
    if not cfg.decodeanhmm:
        return False, "not configured"

    key = (str(cfg.decodeanhmm), str(cfg.phobius_options))
    if key in _native_status:
        return _native_status[key]

    binary = shutil.which(cfg.decodeanhmm) or cfg.decodeanhmm
    status: tuple[bool, str]
    if not Path(binary).is_file():
        status = (False, f"'{cfg.decodeanhmm}' does not exist")
    elif not os.access(binary, os.X_OK):
        status = (False, f"'{binary}' is not executable (chmod +x it on the volume)")
    elif not cfg.phobius_options.is_file():
        status = (
            False,
            f"its options file is missing at {cfg.phobius_options}; without it "
            f"decodeanhmm emits output this service cannot parse",
        )
    else:
        # Probe with a real prediction rather than `-h`, which exits 1 even on a
        # healthy binary. This also proves the options file and model are usable
        # and that the output can still be parsed.
        try:
            proc = subprocess.run(
                [binary, "-f", str(cfg.phobius_options), str(cfg.model)],
                input=f">probe\n{_PROBE_SEQUENCE}\n",
                capture_output=True, text=True, timeout=60,
            )
        except OSError as exc:
            status = (
                False,
                f"'{binary}' could not be executed ({exc.strerror}). The historical "
                f"build is 32-bit and cannot run in this 64-bit image; rebuild it "
                f"for 64-bit from the anhmm source.",
            )
        except subprocess.TimeoutExpired:
            status = (False, f"'{binary}' did not respond to a probe prediction")
        else:
            labels = "".join(
                line[2:].replace(" ", "")
                for line in proc.stdout.splitlines()
                if len(line) > 2 and line[0] == "?" and line[2] == " "
            )
            if len(labels) == len(_PROBE_SEQUENCE):
                status = (True, "ok")
            else:
                detail = (proc.stderr or proc.stdout).strip().splitlines()
                status = (
                    False,
                    f"a probe prediction did not return usable output "
                    f"({detail[-1] if detail else f'exit {proc.returncode}'})",
                )

    _native_status[key] = status
    return status


#: Cache for :func:`homology_search_status`, keyed by database path.
_homology_status: dict[str, tuple[bool, str]] = {}


def homology_search_status(cfg: "Settings") -> tuple[bool, str]:
    """Whether the BLAST homology search is usable, and why not if it isn't.

    The database files existing is not enough: it must also be readable, and
    built with ``-parse_seqids`` so full subject sequences can be retrieved.
    Checking only for a file would advertise the "search for homologues" option
    in the form and then fail every search, so the database is opened first.
    """
    for tool in (cfg.blastp, cfg.blastdbcmd, cfg.kalign):
        if not shutil.which(tool):
            return False, f"'{tool}' is not installed"

    if not any(Path(f"{cfg.blast_db}{ext}").is_file() for ext in (".pin", ".pal")):
        searched = ", ".join(str(d / "swissprot") for d in data_dirs())
        return False, f"no BLAST database; looked for {searched}"

    key = str(cfg.blast_db)
    if key in _homology_status:
        return _homology_status[key]

    try:
        proc = subprocess.run(
            [cfg.blastdbcmd, "-db", str(cfg.blast_db), "-info"],
            capture_output=True, text=True, timeout=60,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        status = (False, f"could not read {cfg.blast_db}: {exc}")
    else:
        if proc.returncode == 0:
            status = (True, "ok")
        else:
            detail = (proc.stderr or proc.stdout).strip().splitlines()
            status = (
                False,
                f"{cfg.blast_db} could not be opened ("
                f"{detail[-1] if detail else 'unknown error'}). Rebuild it with "
                f"scripts/build_swissprot_db.sh, which passes -parse_seqids.",
            )

    _homology_status[key] = status
    return status


#: Extensions offered for download. Anything else in the directory is ignored,
#: so a stray file cannot be published by dropping it on the volume.
DOWNLOAD_SUFFIXES = (".tar.gz", ".tgz", ".zip", ".pdf")


def available_downloads(cfg: "Settings") -> list[Path]:
    """Files offered on the download page, newest name order.

    Returning concrete paths discovered by scanning means no request value is
    ever joined onto a filesystem path, so the download route cannot be walked
    out of this directory.
    """
    if not cfg.downloads_dir.is_dir():
        return []
    return sorted(
        f for f in cfg.downloads_dir.iterdir()
        if f.is_file() and f.name.endswith(DOWNLOAD_SUFFIXES)
    )


settings = Settings()
