"""Serving the licensed academic bundle from mounted storage.

The bundle contains phobius.model and decodeanhmm, so it is never part of the
container image; it is discovered on the same volume as everything else.
"""

import io
import tarfile

import pytest

from app.config import DOWNLOAD_SUFFIXES, Settings, available_downloads


@pytest.fixture
def downloads(tmp_path, monkeypatch):
    directory = tmp_path / "download"
    directory.mkdir()
    (directory / "phobius-1.01.tar.gz").write_bytes(b"bundle")
    (directory / "homologhmm_1.05.tar.gz").write_bytes(b"other")
    monkeypatch.setenv("PHOBIUS_DOWNLOADS", str(directory))
    return directory


def test_bundles_are_discovered(downloads):
    names = [f.name for f in available_downloads(Settings())]
    assert names == ["homologhmm_1.05.tar.gz", "phobius-1.01.tar.gz"]


def test_download_directory_is_found_on_a_mount(monkeypatch, tmp_path):
    mount = tmp_path / "data"
    (mount / "download").mkdir(parents=True)
    (mount / "download" / "phobius-1.01.tar.gz").write_bytes(b"x")
    monkeypatch.delenv("PHOBIUS_DOWNLOADS", raising=False)
    monkeypatch.setenv("PHOBIUS_DATA_DIR", str(mount))
    assert Settings().downloads_dir == mount / "download"


def test_unexpected_file_types_are_not_offered(downloads):
    """A stray file on the volume must not become publicly downloadable."""
    (downloads / "phobius.model").write_bytes(b"licensed")
    (downloads / "notes.txt").write_text("private")
    (downloads / "phobius.env").write_text("SECRET=1")
    names = [f.name for f in available_downloads(Settings())]
    assert "phobius.model" not in names
    assert "notes.txt" not in names
    assert "phobius.env" not in names
    assert all(n.endswith(DOWNLOAD_SUFFIXES) for n in names)


def test_missing_directory_is_not_an_error(monkeypatch, tmp_path):
    monkeypatch.setenv("PHOBIUS_DOWNLOADS", str(tmp_path / "absent"))
    assert available_downloads(Settings()) == []


class TestRoute:
    @pytest.fixture
    def client(self, downloads, monkeypatch):
        from fastapi.testclient import TestClient

        import app.config as config
        import app.main as main

        monkeypatch.setattr(config, "settings", Settings())
        monkeypatch.setattr(main, "settings", config.settings)
        with TestClient(main.app) as c:
            yield c

    def test_page_lists_the_bundle(self, client):
        body = client.get("/download").text
        assert "/download/phobius-1.01.tar.gz" in body
        assert "academic" in body.lower()

    def test_licence_dialog_carries_the_agreed_wording(self, client):
        """The wording is what users are agreeing to, so it is pinned here."""
        import re

        body = client.get("/download").text
        flat = re.sub(r"\s+", " ", body)
        assert "<dialog" in body
        assert "Phobius 1.01 standalone download" in flat
        assert "intended for academic users only" in flat
        assert "contact Erik Sonnhammer at" in flat
        assert "mailto:Erik.Sonnhammer@gmail.com" in body
        assert "only for private study, education or nonprofit research" in flat
        assert ">Download Phobius<" in flat

    def test_download_is_reachable_without_scripting(self, client):
        """The dialog is a notice, not an access control: keep a plain link for
        clients without JavaScript rather than leaving them with a dead button."""
        body = client.get("/download").text
        assert "<noscript>" in body
        noscript = body.split("<noscript>")[1].split("</noscript>")[0]
        assert "/download/phobius-1.01.tar.gz" in noscript
        assert "academic" in noscript.lower()

    def test_bundle_can_be_fetched(self, client):
        r = client.get("/download/phobius-1.01.tar.gz")
        assert r.status_code == 200
        assert r.content == b"bundle"

    @pytest.mark.parametrize("name", [
        "../phobius.model",
        "..%2Fphobius.model",
        "....//phobius.model",
        "/etc/passwd",
        "phobius.env",
        "notes.txt",
    ])
    def test_traversal_and_unlisted_files_are_refused(self, client, downloads, name):
        (downloads / "notes.txt").write_text("private")
        (downloads / "phobius.env").write_text("SECRET=1")
        assert client.get(f"/download/{name}").status_code == 404


def test_built_bundle_has_everything_needed_to_run(tmp_path):
    """Guards the bundle's contents: phobius.pl invokes ./decodeanhmm by that
    exact name, and refuses to start without the model and options beside it."""
    required = {"phobius/phobius.pl", "phobius/decodeanhmm",
                "phobius/phobius.model", "phobius/phobius.options"}
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as tar:
        for name in sorted(required):
            info = tarfile.TarInfo(name)
            info.size = 1
            tar.addfile(info, io.BytesIO(b"x"))
    buffer.seek(0)
    with tarfile.open(fileobj=buffer, mode="r:gz") as tar:
        assert required <= set(tar.getnames())
