# Phobius web service

Web service for [Phobius](https://doi.org/10.1016/j.jmb.2004.03.016), a combined
transmembrane topology and signal peptide predictor, and for
[PolyPhobius](https://doi.org/10.1093/bioinformatics/bti1014), its
homology-supported variant.

This replaces the Perl CGI server that ran from 2004. The prediction engines are
unchanged — the labels this service returns are byte-for-byte identical to those
of the legacy engine, and `tests/test_golden.py` enforces that — but everything
around them has been rebuilt as a single container.

## Quick start

```bash
docker run -p 8000:8000 -v /path/to/licensed:/mnt/data:ro ghcr.io/statisticalbiotechnology/phobius:latest
```

The mounted directory must contain `phobius.model`. See **Licensed files** below.

## Layout

| Path | What it is |
|---|---|
| `app/fasta.py` | FASTA parsing and residue normalisation |
| `app/features.py` | label strings → feature tables and topology strings |
| `app/engines.py` | subprocess wrappers around the prediction engines |
| `app/homology.py` | DIAMOND + Kalign pipeline for PolyPhobius |
| `app/plot.py` | posterior probability plot, rendered as inline SVG |
| `app/models.py` | request validation and submission limits |
| `app/main.py` | FastAPI routes and templates |
| `engine/` | `homologhmm.jar` (GPL), `biojava.jar` (LGPL) |
| `tests/golden/` | reference output captured from the legacy engine |

## Licensed files

**Nothing licensed is in the container image**, which is what lets the image be
published freely.

| File | Terms | Where it lives |
|---|---|---|
| `phobius.model` | academic licence | mounted at `/mnt/data`, `PHOBIUS_MODEL` |
| `decodeanhmm` | licensed executable | optional, not required — see below |
| `homologhmm.jar` | GPL | in the image |
| `biojava.jar` | LGPL | in the image |

The Java engine covers every prediction mode, so `decodeanhmm` is not needed at
all. It is roughly 8× faster for plain single-sequence predictions (41 ms versus
356 ms, dominated by JVM startup; throughput differs by only ~1.4×). If you have
a licensed copy on the mounted volume you can enable it by setting
`PHOBIUS_DECODEANHMM` — put `phobius.options` beside it. The golden tests assert
both engines produce identical output.

## Configuration

| Variable | Default | Meaning |
|---|---|---|
| `PHOBIUS_MODEL` | `/mnt/data/phobius.model` | licensed model file; must be named `phobius.model` |
| `PHOBIUS_ENGINE_DIR` | `/app/engine` | directory holding the jars |
| `PHOBIUS_DIAMOND_DB` | `/mnt/data/swissprot.dmnd` | homology database; absent disables that mode |
| `PHOBIUS_DECODEANHMM` | unset | optional licensed native engine |
| `PHOBIUS_MAX_SEQUENCES` | `100` | sequences per request |
| `PHOBIUS_MAX_RESIDUES` | `50000` | residues per request (~20 s of compute) |
| `PHOBIUS_MAX_RESIDUES_SINGLE` | `10000` | residues in one sequence |
| `PHOBIUS_ENGINE_TIMEOUT` | `120` | seconds before an engine is killed |
| `PHOBIUS_MAX_CONCURRENCY` | `2` | concurrent predictions; match the CPU allocation |

The model must keep the filename `phobius.model`: the Java engine loads it as a
classpath resource by that exact name, not as an argument.

## Homology search (optional)

PolyPhobius works two ways. Supplying your own aligned FASTA needs no database
and is the reproducible option. The automatic search needs a DIAMOND database:

```bash
./scripts/build_swissprot_db.sh /mnt/data
```

That produces a ~250 MB index. The original server searched UniProt/TrEMBL with
legacy BLAST; Swiss-Prot with DIAMOND is a different search, so predictions from
this path will not reproduce historical ones.

## Publishing the image

`.github/workflows/ci.yml` runs the tests, then builds the image, then verifies
it, and only then pushes to `ghcr.io/<owner>/phobius`. Pushes to `master`
publish `latest` and a short-SHA tag; a `v*` git tag additionally publishes
semver tags:

```bash
git tag v2.0.0 && git push --tags     # -> :2.0.0, :2.0, :latest
```

Three checks run against the built image **before** anything is pushed, because
a licensed file that reached a public registry could not be un-published:

| Check | Fails if |
|---|---|
| licensed files | `phobius.model`, `phobius.options` or `decodeanhmm` is in the image |
| runtime user | the image does not run as UID 1000 |
| smoke test | the container does not boot and serve `/healthz`, `/`, `/docs` |

The smoke test uses a placeholder model, so CI never needs the licensed file.

> **On the first push the GHCR package is private.** Serve can only pull public
> images, so go to *Packages → phobius → Package settings → Change visibility →
> Public* once. Later pushes keep that setting.

The image is amd64 only. To also build arm64, add
`platforms: linux/amd64,linux/arm64` to the push step — expect a much longer
build, since the arm64 layers are emulated.

## Deploying on SciLifeLab Serve

The image is built for the platform's constraints: it runs as **UID 1000**,
listens on **port 8000** (Serve allows 3000–9999), holds no licensed or private
data, and needs no persistent volume unless you enable the homology search.

1. Publish the image and make the GHCR package public, as above.
2. In your project, define a storage mount and place `phobius.model` there —
   Serve's file storage is not public, so the licensed model stays private.
3. Create an app of type *Other*, point it at
   `ghcr.io/<owner>/phobius:latest`, set the port to 8000, and select the
   mount path.
4. `GET /healthz` returns 503 with a specific message if the model is missing.
   Without it the container refuses to start rather than failing per-request.

Default resources (2 vCPU, 4 GB) are ample. If you raise them, raise
`PHOBIUS_MAX_CONCURRENCY` to match.

## API

`POST /api/predict` with `{"sequence": "<FASTA>"}` returns structured
predictions; add `"plot": true` for an SVG per sequence. `POST /api/predict.txt`
returns the feature table as plain text, `?format=short` for one line per
protein. Interactive documentation is at `/docs`.

The legacy `POST /cgi-bin/predict.pl` endpoint still works with its original
field names, so existing bookmarks and scripts keep running.

## Development

```bash
pip install -e ".[dev]"
PHOBIUS_MODEL=engine/phobius.model uvicorn app.main:app --reload
pytest
```

Tests that need an engine skip cleanly when the model is absent. `pytest` also
skips the homology tests unless `diamond` and `kalign` are on `PATH`; the
container has both.

### Changing prediction behaviour

`tests/golden/labels.tsv` holds labels captured from the legacy engine for a
52-sequence corpus covering edge cases (1 residue, 3000 residues, all-hydrophobic,
ambiguous residues, hostile FASTA headers). If a change makes those tests fail,
the change altered predictions. Regenerate them only with
`tests/generate_golden.py`, and only against the legacy engine.

Two deliberate quirks are preserved for byte-compatibility and are covered by
tests: a two-region prediction reports a non-cytoplasmic loop as `CYTOPLASMIC.`,
and the short-format header retains its original spelling. Both are noted in the
code.

## Citing

Käll L, Krogh A, Sonnhammer ELL. *A combined transmembrane topology and signal
peptide prediction method.* J Mol Biol 338:1027–1036 (2004).

Käll L, Krogh A, Sonnhammer ELL. *An HMM posterior decoder for sequence feature
prediction that includes homology information.* Bioinformatics 21(Suppl 1):i251–i257 (2005).
