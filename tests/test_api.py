"""End-to-end tests, including the specific defects found in predict.pl."""

import pytest

from conftest import needs_engine

pytestmark = needs_engine

SEQ = (
    "MYGKIIFVLLLSAIVSISASSTTGVAMHTSTSSSVTKSYISSQTNDTHKRDTYAATPRAHEVSEISVRT"
    "VYPPEEETGERVQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKKSPSDVKPLPSPDTDVPLS"
    "SVEIENPETSDQ"
)


@pytest.mark.parametrize("path", ["/", "/constrained", "/poly", "/instructions",
                                  "/download", "/api", "/healthz"])
def test_pages_render(client, path):
    assert client.get(path).status_code == 200


def test_json_prediction(client):
    r = client.post("/api/predict", json={"sequence": f">glpa\n{SEQ}"})
    assert r.status_code == 200
    prediction = r.json()["predictions"][0]
    assert prediction["id"] == "glpa"
    assert prediction["length"] == 150
    assert prediction["transmembrane_count"] == 1
    assert prediction["signal_peptide"] is True
    assert prediction["cleavage_site"] == 19
    assert prediction["topology"] == "n5-15c19/20o92-114i"


def test_json_prediction_with_plot(client):
    r = client.post("/api/predict", json={"sequence": SEQ, "plot": True})
    svg = r.json()["predictions"][0]["svg"]
    assert svg.startswith("<svg") and svg.endswith("</svg>")


def test_text_api_short_and_long(client):
    short = client.post("/api/predict.txt?format=short", json={"sequence": f">glpa\n{SEQ}"})
    assert short.text.splitlines()[1].startswith("glpa")
    long = client.post("/api/predict.txt", json={"sequence": f">glpa\n{SEQ}"})
    assert long.text.startswith("ID   glpa\nFT   SIGNAL")


def test_legacy_cgi_route_still_works(client):
    """Existing bookmarks and scripts post to /cgi-bin/predict.pl."""
    r = client.post("/cgi-bin/predict.pl",
                    data={"protseq": f">glpa\n{SEQ}", "format": "short"})
    assert r.status_code == 200
    assert "n5-15c19/20o92-114i" in r.text


def test_file_upload_takes_precedence_over_textarea(client):
    r = client.post("/predict",
                    data={"protseq": ">ignored\nMKKL", "format": "short"},
                    files={"protfile": ("q.fa", f">fromfile\n{SEQ}", "text/plain")})
    assert "fromfile" in r.text and "ignored" not in r.text


# --- the vulnerabilities that were live in predict.pl ---------------------

def test_hostile_fasta_header_is_escaped_not_executed(client):
    """predict.pl:540 wrote the FASTA header into HTML unescaped."""
    r = client.post("/predict",
                    data={"protseq": f'>"><script>alert(1)</script>\n{SEQ}',
                          "format": "plp"})
    assert r.status_code == 200
    assert "<script>alert(1)</script>" not in r.text
    assert "&lt;script&gt;" in r.text


def test_backquotes_in_header_reach_no_shell(client, tmp_path):
    """predict.pl:482 interpolated the header into a gnuplot script, where
    backquoted text is substituted with the output of a shell command."""
    marker = tmp_path / "pwned"
    r = client.post("/predict",
                    data={"protseq": f">`touch {marker}`\n{SEQ}", "format": "plp"})
    assert r.status_code == 200
    assert not marker.exists()


def test_constraints_reject_shell_metacharacters(client):
    r = client.post("/predict", data={"constrained": "Y", "protseq": SEQ,
                                      "M": "$(id)", "format": "nog"})
    assert r.status_code == 400
    assert "Could not read" in r.text


def test_constraints_reject_out_of_range_positions(client):
    r = client.post("/predict", data={"constrained": "Y", "protseq": SEQ,
                                      "M": "9999", "format": "nog"})
    assert r.status_code == 400
    assert "beyond the end" in r.text


def test_constraints_are_honoured(client):
    """A constraint says a residue *must* be in that state, not that the region
    must stop there -- the engine may extend a helix past the constrained span."""
    import re

    r = client.post("/predict", data={"constrained": "Y", "protseq": SEQ,
                                      "M": "40-60", "format": "nog"})
    assert r.status_code == 200
    helices = [
        (int(a), int(b))
        for a, b in re.findall(r"FT   TRANSMEM\s+(\d+)\s+(\d+)", r.text)
    ]
    assert any(start <= 40 and stop >= 60 for start, stop in helices), helices


def test_unconstrained_prediction_differs_from_constrained(client):
    """Guards against the constraints being silently dropped on the way to the
    engine -- which is exactly what happened when they were passed joined."""
    plain = client.post("/predict", data={"protseq": SEQ, "format": "nog"}).text
    forced = client.post("/predict", data={"constrained": "Y", "protseq": SEQ,
                                           "M": "40-60", "format": "nog"}).text
    assert "FT   TRANSMEM" in plain and "FT   TRANSMEM" in forced
    assert plain != forced


# --- limits ---------------------------------------------------------------

def test_too_many_residues_is_rejected(client):
    from app.config import settings

    big = ">x\n" + "A" * (settings.max_residues_single + 1)
    r = client.post("/api/predict", json={"sequence": big})
    assert r.status_code == 400
    assert "residues" in r.json()["detail"]


def test_too_many_sequences_is_rejected(client):
    from app.config import settings

    many = "".join(f">s{i}\nMKKLLA\n" for i in range(settings.max_sequences + 1))
    r = client.post("/api/predict", json={"sequence": many})
    assert r.status_code == 400


def test_empty_submission_is_rejected(client):
    r = client.post("/predict", data={"protseq": "", "format": "short"})
    assert r.status_code == 400
    assert "No sequences found" in r.text


def test_unknown_output_format_is_rejected(client):
    r = client.post("/predict", data={"protseq": SEQ, "format": "bogus"})
    assert r.status_code == 400


def test_plot_only_requires_an_alignment(client):
    r = client.post("/predict", data={"protseq": SEQ, "format": "aplp"})
    assert r.status_code == 400


# --- PolyPhobius ----------------------------------------------------------

def test_polyphobius_on_supplied_alignment(client):
    alignment = f">query\n{SEQ}\n>homolog\n{SEQ[:70]}A{SEQ[71:]}\n"
    r = client.post("/predict", data={"poly": "Y", "informat": "align",
                                      "protseq": alignment, "format": "nog"})
    assert r.status_code == 200
    assert "ID   query" in r.text


def test_polyphobius_rejects_ragged_alignment(client):
    r = client.post("/predict", data={"poly": "Y", "informat": "align",
                                      "protseq": ">a\nMKKL\n>b\nMKK\n",
                                      "format": "nog"})
    assert r.status_code == 400
    assert "same length" in r.text


# --- instructions page ----------------------------------------------------

def test_instructions_carry_the_example_sequence(client):
    """Users asked for the example the old server had, so it is pinned here."""
    body = client.get("/instructions").text
    assert "Q8TCT8|PSL2_HUMAN" in body
    assert "MGPQRRLSPAGAALLWGFLLQLTAAQEAILHASGNGTTKDYCMLYNPYWTALPSTLENAT" in body


def test_example_sequence_is_valid_input(client):
    """The documented example must actually run, or it is worse than nothing."""
    import re

    body = client.get("/instructions").text
    block = re.search(r'<pre id="example-fasta">(.*?)</pre>', body, re.S).group(1)
    fasta = block.replace("&gt;", ">").replace("&amp;", "&")

    r = client.post("/api/predict", json={"sequence": fasta})
    assert r.status_code == 200
    prediction = r.json()["predictions"][0]
    assert prediction["id"] == "Q8TCT8|PSL2_HUMAN"
    assert prediction["length"] == 520
    assert prediction["transmembrane_count"] > 0


def test_example_has_a_copy_button(client):
    body = client.get("/instructions").text
    assert 'class="copy-example' in body
    assert 'data-target="example-fasta"' in body
    assert 'id="example-fasta"' in body
