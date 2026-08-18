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
| `app/homology.py` | BLAST+ and Kalign pipeline for PolyPhobius |
| `app/plot.py` | posterior probability plot, rendered as inline SVG |
| `app/models.py` | request validation and submission limits |
| `app/main.py` | FastAPI routes and templates |
| `start-script.sh` | container entry point; Serve runs this from the WORKDIR |
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

### The optional native engine

The Java engine covers every prediction mode, so `decodeanhmm` is **not needed**.
It is about 8× faster for a single sequence (41 ms versus 356 ms, almost all of
it JVM startup), but only ~1.4× on throughput, so on a web form the difference is
imperceptible. Most deployments should skip it.

If you do want it, put all three files on the same mounted storage:

```
/home/data/phobius.model      # required
/home/data/phobius.options    # required by decodeanhmm; looked for beside the model
/home/data/decodeanhmm        # the licensed binary, chmod +x
/home/data/phobius.env        # PHOBIUS_DECODEANHMM=/home/data/decodeanhmm
```

`phobius.options` is looked for next to the model, not next to the binary;
`PHOBIUS_OPTIONS` overrides that.

**The historical binary will not work.** It is a 32-bit build and the image is
64-bit with no i386 libraries, so it cannot be executed at all. Rebuild it for
64-bit from the anhmm source first.

Nothing here can break the service. At startup the binary is probed with a real
prediction, and if anything is wrong — missing, not executable, wrong
architecture, missing options file — the reason is logged and the Java engine is
used instead:

```
WARNING  phobius: native fast path requested but unusable, falling back to the
Java engine: '/home/data/decodeanhmm' could not be executed (No such file or
directory). The historical build is 32-bit and cannot run in this 64-bit image;
rebuild it for 64-bit from the anhmm source.
```

The golden tests assert both engines produce byte-identical output, so switching
between them cannot change results.

## Configuration

| Variable | Default | Meaning |
|---|---|---|
| `PHOBIUS_MODEL` | *discovered* | full path to the licensed model; must be named `phobius.model` |
| `PHOBIUS_DATA_DIR` | *unset* | directory searched before the defaults |
| `PHOBIUS_ENGINE_DIR` | `/app/engine` | directory holding the jars |
| `PHOBIUS_BLAST_DB` | *discovered* | BLAST database prefix; absent disables that mode |
| `PHOBIUS_DECODEANHMM` | unset | optional licensed native engine |
| `PHOBIUS_MAX_SEQUENCES` | `100` | sequences per request |
| `PHOBIUS_MAX_RESIDUES` | `50000` | residues per request (~20 s of compute) |
| `PHOBIUS_MAX_RESIDUES_SINGLE` | `10000` | residues in one sequence |
| `PHOBIUS_ENGINE_TIMEOUT` | `120` | seconds before an engine is killed |
| `PHOBIUS_MAX_CONCURRENCY` | `2` | concurrent predictions; match the CPU allocation |

The model must keep the filename `phobius.model`: the Java engine loads it as a
classpath resource by that exact name, not as an argument.

### Finding the model

You normally do not need to configure anything. Mount the storage holding
`phobius.model` anywhere, and it is found by searching, in order:

1. `PHOBIUS_MODEL`, if set — a full path, and always wins
2. `PHOBIUS_DATA_DIR`, if set
3. `/mnt/data`, `/home/data`, `/data`, `/app/data`
4. `engine/phobius.model` beside the source, for local development

The same search locates the BLAST database, keying off `swissprot.pin`
(or `swissprot.pal` for a split database) and using its prefix. The image deliberately does **not**
pin `PHOBIUS_MODEL`, because an explicit value beats discovery and would break
every deployment whose storage is mounted somewhere else.

If the model cannot be found, the container exits with a message listing every
path it tried *and what each directory actually contains*, which is usually
enough to spot a wrong mount path or a misnamed file.

### Settings without environment variables

SciLifeLab Serve provides no way to set environment variables for a custom app.
Any setting can instead go in a `phobius.env` file placed on the mounted storage,
alongside the model — `start-script.sh` loads it at boot:

```sh
# phobius.env
PHOBIUS_MODEL=/home/data/licensed/phobius.model
PHOBIUS_MAX_SEQUENCES=25
PHOBIUS_MAX_CONCURRENCY=4
```

Only `KEY=VALUE` lines are honoured; anything else is ignored rather than
executed, so the file cannot be used to run arbitrary commands from the volume.

## Homology search (optional)

PolyPhobius works two ways. Supplying your own aligned FASTA needs no database
and is the reproducible option. The automatic search needs a BLAST database on
the mounted storage:

```bash
./scripts/build_swissprot_db.sh /home/data
```

That downloads Swiss-Prot and produces ~340 MB of `swissprot.*` files, well
inside the 5 GB volume cap. The option appears in the PolyPhobius form by itself
once the database is readable.

`-parse_seqids` is required and the script passes it: without it the search still
runs but `blastdbcmd` cannot retrieve the full subject sequences the aligner
needs. A database missing it is reported at startup rather than failing per
request.

The original server searched UniProt/TrEMBL, so predictions from this path will
not reproduce historical ones.

### Why BLAST rather than DIAMOND

DIAMOND is the faster tool for bulk searches, but its cost for a *single* query
is dominated by building a seed index over the whole database. Measured against
Swiss-Prot on this image, one query:

| | time | homologues found |
|---|---|---|
| DIAMOND `--very-sensitive` | 42.8 s | 7 |
| DIAMOND default, `-c1` | 16.0 s | 5 |
| **blastp** | **1.2 s** | **8** |

blastp is 36× faster and returned a strict superset. End to end the request went
from 78 s to 1.25 s. The database location makes no difference — reading it from
the mounted volume cost 16.0 s versus 17.5 s from container-local disk, so this
is compute, not I/O.

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
| start-up script | `./start-script.sh` is missing or not executable at the WORKDIR |
| smoke test | the container does not boot and serve `/healthz`, `/`, `/docs` |

The smoke test launches the container the way Serve does — overriding the
entrypoint and invoking `./start-script.sh` by path — rather than relying on the
`ENTRYPOINT`, because those two paths can differ and only the first is what
actually runs in production.

The smoke test uses a placeholder model, so CI never needs the licensed file.

> **On the first push the GHCR package is private.** Serve can only pull public
> images, so go to *Packages → phobius → Package settings → Change visibility →
> Public* once. Later pushes keep that setting.

Build provenance attestation is attempted only when the repository itself is
public — GitHub requires a public repo (or a paid plan) to persist attestations,
and the step runs *after* the image is pushed, so failing it would mark a
successful publish as broken. It starts working on its own once the repository
is made public.

The image is amd64 only. To also build arm64, add
`platforms: linux/amd64,linux/arm64` to the push step — expect a much longer
build, since the arm64 layers are emulated.

## Deploying on SciLifeLab Serve

The image is built for the platform's constraints: it runs as **UID 1000**,
listens on **port 8000** (Serve allows 3000–9999), starts through an executable
**`./start-script.sh` at the WORKDIR** (`/app`), which is how Serve launches a
container, holds no licensed or private data, and needs no persistent volume
unless you enable the homology search.

1. Publish the image and make the GHCR package public, as above.
2. In your project, define a storage mount and place `phobius.model` there —
   Serve's file storage is not public, so the licensed model stays private. Any
   mount path works; `/home/data` and `/mnt/data` are found automatically, and
   anything else can be pointed to with a `phobius.env` file on the same storage.
3. Create an app of type *Other*, point it at
   `ghcr.io/<owner>/phobius:latest`, set the port to 8000, and select the
   mount path.
4. `GET /healthz` returns 503 with a specific message if the model is missing.
   Without it the container refuses to start rather than failing per-request.

Default resources (2 vCPU, 4 GB) are ample. If you raise them, raise
`PHOBIUS_MAX_CONCURRENCY` to match.

## The old address

`phobius.sbc.su.se` has been the published URL since 2004 and appears in cited
methods sections. See [docs/legacy-redirect.md](docs/legacy-redirect.md) for how
to point it at this service — including why a 308 rather than a 301 matters for
scripts that still POST to `/cgi-bin/predict.pl`.

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
skips the homology tests unless `ncbi-blast+` and `kalign` are on `PATH`; the
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
