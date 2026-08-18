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
from dataclasses import dataclass, field
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent

#: Directories probed for mounted data files, in order. These are the paths
#: hosting platforms commonly mount project storage at.
DATA_DIRS: tuple[str, ...] = ("/mnt/data", "/home/data", "/data", "/app/data")


def data_dirs() -> list[Path]:
    """Directories searched for mounted data files, most specific first."""
    candidates: list[str] = []
    explicit = os.environ.get("PHOBIUS_DATA_DIR")
    if explicit:
        candidates.append(explicit)
    candidates.extend(DATA_DIRS)
    # dict.fromkeys keeps order while removing duplicates.
    return [Path(d) for d in dict.fromkeys(candidates)]


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
    diamond: str = field(default_factory=lambda: os.environ.get("PHOBIUS_DIAMOND", "diamond"))
    kalign: str = field(default_factory=lambda: os.environ.get("PHOBIUS_KALIGN", "kalign"))
    diamond_db: Path = field(default_factory=lambda: discover(
        "swissprot.dmnd", "PHOBIUS_DIAMOND_DB", _ROOT / "data" / "swissprot.dmnd"))

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
        return bool(shutil.which(self.diamond)) and bool(shutil.which(self.kalign)) and self.diamond_db.exists()

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
        if self.decodeanhmm:
            if not shutil.which(self.decodeanhmm) and not Path(self.decodeanhmm).is_file():
                problems.append(f"PHOBIUS_DECODEANHMM is set but '{self.decodeanhmm}' is not executable.")
            elif not self.phobius_options.is_file():
                problems.append(
                    f"PHOBIUS_DECODEANHMM is set but its options file is missing at "
                    f"{self.phobius_options}; decodeanhmm would emit unparseable output."
                )
        return problems


settings = Settings()
