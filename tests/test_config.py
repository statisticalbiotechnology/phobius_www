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
    for var in ("PHOBIUS_MODEL", "PHOBIUS_DATA_DIR", "PHOBIUS_DIAMOND_DB"):
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
    assert [str(d) for d in data_dirs()[1:]] == list(DATA_DIRS)


def test_image_does_not_pin_the_data_paths():
    """The Dockerfile must not set PHOBIUS_MODEL or PHOBIUS_DIAMOND_DB.

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
    assert "PHOBIUS_DIAMOND_DB=" not in body


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


def test_diamond_database_is_discovered_alongside_the_model(monkeypatch, tmp_path):
    mount = tmp_path / "data"
    mount.mkdir()
    (mount / "swissprot.dmnd").write_bytes(b"x")
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(mount))
    assert Settings().diamond_db == mount / "swissprot.dmnd"
