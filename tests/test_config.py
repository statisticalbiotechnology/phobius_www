"""Locating the mounted model.

Hosting platforms mount project storage wherever the operator configures it, and
SciLifeLab Serve provides no way to set an environment variable for a custom
app. The service therefore has to find the licensed model without being told
where it is.
"""

import pytest

from app.config import DATA_DIRS, Settings, data_dirs, discover


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    # Clear every variable these tests reason about, so a deployment-shaped
    # environment (e.g. PHOBIUS_DECODEANHMM set on the host) cannot leak in.
    for var in ("PHOBIUS_MODEL", "PHOBIUS_DATA_DIR",
                "PHOBIUS_DECODEANHMM", "PHOBIUS_OPTIONS", "PHOBIUS_BLAST_DB"):
        monkeypatch.delenv(var, raising=False)


@pytest.fixture
def no_mounts(monkeypatch, tmp_path):
    """Point the default search path at directories that do not exist.

    Without this, tests asserting "nothing is mounted" pass on a laptop and fail
    inside a container that really does have /mnt/data or /home/data.
    """
    import app.config as config

    empty = tuple(str(tmp_path / "absent" / name) for name in ("a", "b"))
    monkeypatch.setattr(config, "DATA_DIRS", empty)
    return empty


def test_explicit_path_wins(monkeypatch, tmp_path):
    explicit = tmp_path / "elsewhere.model"
    monkeypatch.setenv("PHOBIUS_MODEL", str(explicit))
    assert discover("phobius.model", "PHOBIUS_MODEL", tmp_path / "fallback") == explicit


def test_data_dir_is_searched_first(monkeypatch, tmp_path):
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(tmp_path))
    assert data_dirs()[0] == tmp_path
    assert data_dirs()[1] == tmp_path / "db"


def test_image_does_not_pin_the_data_paths():
    """The Dockerfile must not set PHOBIUS_MODEL or PHOBIUS_BLAST_DB.

    An explicit value beats discovery, so pinning them in the image breaks every
    deployment whose storage is mounted anywhere other than the pinned path.
    """
    import pathlib

    dockerfile = pathlib.Path(__file__).resolve().parent.parent / "Dockerfile"
    body = "\n".join(
        line for line in dockerfile.read_text().splitlines()
        if not line.lstrip().startswith("#")
    )
    assert "PHOBIUS_MODEL=" not in body
    assert "PHOBIUS_BLAST_DB=" not in body


def test_model_found_on_a_mounted_directory(monkeypatch, tmp_path):
    mount = tmp_path / "home" / "data"
    mount.mkdir(parents=True)
    (mount / "phobius.model").write_bytes(b"x")
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(mount))

    found = discover("phobius.model", "PHOBIUS_MODEL", tmp_path / "fallback")
    assert found == mount / "phobius.model"


def test_falls_back_when_nothing_is_mounted(tmp_path, no_mounts):
    fallback = tmp_path / "engine" / "phobius.model"
    assert discover("phobius.model", "PHOBIUS_MODEL", fallback) == fallback


def test_db_subdirectory_is_searched(monkeypatch, tmp_path):
    """A BLAST database is nine files; keeping it in db/ keeps the mount tidy."""
    mount = tmp_path / "data"
    (mount / "db").mkdir(parents=True)
    (mount / "db" / "swissprot.pin").write_bytes(b"x")
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(mount))
    assert Settings().blast_db == mount / "db" / "swissprot"


def test_mount_is_searched_before_its_db_subdirectory(monkeypatch, tmp_path):
    mount = tmp_path / "data"
    (mount / "db").mkdir(parents=True)
    (mount / "swissprot.pin").write_bytes(b"x")
    (mount / "db" / "swissprot.pin").write_bytes(b"x")
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(mount))
    assert Settings().blast_db == mount / "swissprot"


def test_model_may_also_live_in_the_db_subdirectory(monkeypatch, tmp_path):
    """The same rule applies to every mounted file, not just the database."""
    mount = tmp_path / "data"
    (mount / "db").mkdir(parents=True)
    (mount / "db" / "phobius.model").write_bytes(b"x")
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(mount))
    settings = Settings()
    assert settings.model == mount / "db" / "phobius.model"
    # The Java engine loads the model as a classpath resource, so wherever it is
    # found, that directory has to reach the classpath.
    assert str(mount / "db") in settings.classpath


def test_every_mount_gets_a_db_subdirectory_entry():
    from app.config import DATA_SUBDIRS

    searched = [str(d) for d in data_dirs()]
    assert DATA_SUBDIRS == ("", "db")
    for mount in DATA_DIRS:
        assert mount in searched
        assert f"{mount}/db" in searched


def test_search_order_has_no_duplicates(monkeypatch):
    # A PHOBIUS_DATA_DIR that repeats a default must not be probed twice.
    monkeypatch.setenv("PHOBIUS_DATA_DIR", DATA_DIRS[0])
    assert len(data_dirs()) == len(set(data_dirs()))


def test_serve_mount_path_is_covered():
    # SciLifeLab Serve deployments mount project storage here.
    assert "/home/data" in DATA_DIRS


def test_missing_model_error_lists_every_path_searched(monkeypatch, tmp_path, no_mounts):
    # Simulate the container: nothing mounted, and no development copy either.
    import app.config as config

    monkeypatch.setattr(config, "_ROOT", tmp_path / "no-such-checkout")
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(tmp_path))
    problems = "\n".join(Settings().check())
    assert "phobius.model not found" in problems
    for directory in (str(tmp_path), *no_mounts):
        assert directory in problems
    assert "phobius.env" in problems


def test_classpath_follows_the_discovered_model(monkeypatch, tmp_path):
    """The Java engine loads the model as a classpath resource, so the directory
    it was discovered in has to end up on the classpath."""
    mount = tmp_path / "data"
    mount.mkdir()
    (mount / "phobius.model").write_bytes(b"x")
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(mount))
    assert str(mount) in Settings().classpath


def test_blast_database_is_discovered_alongside_the_model(monkeypatch, tmp_path):
    """A BLAST database is a set of files sharing a prefix, so discovery keys off
    the index file and returns the prefix, not a filename."""
    mount = tmp_path / "data"
    mount.mkdir()
    (mount / "swissprot.pin").write_bytes(b"x")
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(mount))
    assert Settings().blast_db == mount / "swissprot"


def test_split_blast_database_is_discovered(monkeypatch, tmp_path):
    # A database large enough to be split has a .pal alias instead of a .pin.
    mount = tmp_path / "data"
    mount.mkdir()
    (mount / "swissprot.pal").write_bytes(b"x")
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(mount))
    assert Settings().blast_db == mount / "swissprot"


class TestNativeEngine:
    """The optional native engine must never be able to break the service.

    It is an accelerator producing byte-identical output, so anything wrong with
    it degrades to the Java engine instead of failing startup or predictions.
    """

    def test_unconfigured_is_not_an_error(self):
        from app.config import native_engine_status

        usable, reason = native_engine_status(Settings())
        assert usable is False
        assert reason == "not configured"

    def test_missing_binary_is_reported_not_raised(self, tmp_path):
        from app.config import native_engine_status

        usable, reason = native_engine_status(Settings(decodeanhmm=str(tmp_path / "absent")))
        assert usable is False
        assert "does not exist" in reason

    def test_missing_options_file_is_reported(self, tmp_path):
        from app.config import native_engine_status

        binary = tmp_path / "decodeanhmm"
        binary.write_text("#!/bin/sh\nexit 0\n")
        binary.chmod(0o755)
        usable, reason = native_engine_status(
            Settings(decodeanhmm=str(binary), _options_override=str(tmp_path / "absent.options"))
        )
        assert usable is False
        assert "options file is missing" in reason

    def test_a_binary_that_cannot_exec_is_reported(self, tmp_path):
        """Reproduces the 32-bit-binary-in-a-64-bit-image case."""
        from app.config import native_engine_status

        binary = tmp_path / "decodeanhmm"
        binary.write_bytes(b"\x7fELF\x01\x01\x01\x00" + b"\x00" * 56)  # unloadable ELF
        binary.chmod(0o755)
        options = tmp_path / "phobius.options"
        options.write_text("N 1\n")

        usable, reason = native_engine_status(
            Settings(decodeanhmm=str(binary), _options_override=str(options))
        )
        assert usable is False
        assert "could not be executed" in reason
        assert "32-bit" in reason

    def test_a_non_executable_file_is_reported(self, tmp_path):
        """Uploading to a volume commonly loses the execute bit."""
        from app.config import native_engine_status

        binary = tmp_path / "decodeanhmm"
        binary.write_text("#!/bin/sh\nexit 0\n")
        binary.chmod(0o644)
        usable, reason = native_engine_status(Settings(decodeanhmm=str(binary)))
        assert usable is False
        assert "not executable" in reason

    def test_check_does_not_make_a_broken_native_engine_fatal(self, tmp_path):
        problems = Settings(decodeanhmm=str(tmp_path / "absent")).check()
        assert not any("DECODEANHMM" in p for p in problems)
