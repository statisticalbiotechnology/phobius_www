import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from app.config import settings  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent

#: Tests that actually invoke an engine need the licensed model to be present.
needs_engine = pytest.mark.skipif(
    bool(settings.check()),
    reason=f"engine not configured: {'; '.join(settings.check())}",
)


@pytest.fixture(scope="session")
def corpus_text() -> str:
    return (HERE / "data" / "regression.fa").read_text()


@pytest.fixture(scope="session")
def golden() -> dict[str, str]:
    """Reference labels captured from the legacy decodeanhmm engine."""
    path = HERE / "golden" / "labels.tsv"
    out = {}
    for line in path.read_text().splitlines():
        name, labels = line.split("\t")
        out[name] = labels
    return out


@pytest.fixture(scope="session")
def client():
    from fastapi.testclient import TestClient

    from app.main import app

    with TestClient(app) as c:
        yield c
