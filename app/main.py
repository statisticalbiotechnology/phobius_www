"""Phobius web service.

Replaces the ``predict.pl`` CGI script. The prediction engines are unchanged --
this is a wrapper, not a reimplementation of the science -- but everything
around them is different:

* templates escape by default, so a hostile FASTA header cannot inject markup;
* subprocesses are invoked as argument lists, never through a shell;
* the posterior plot is built in memory instead of being written, along with the
  user's sequences, to guessable filenames under the web root;
* one engine invocation serves a whole request rather than one per sequence;
* requests are bounded in size and concurrency instead of serialised behind a
  single global lock.
"""

from __future__ import annotations

import logging
import threading
from contextlib import asynccontextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from . import __version__, engines, features, plot
from .config import (
    available_downloads,
    homology_search_status,
    native_engine_status,
    settings,
)
from .fasta import FastaError, Record, parse, parse_alignment
from .homology import HomologyError, search_and_align
from .models import (
    ApiRequest,
    ConstraintError,
    InputFormat,
    OutputFormat,
    SubmissionError,
    build_constraints,
    check_limits,
    constraint_tokens,
)

log = logging.getLogger("phobius")

# Uvicorn installs handlers for its own loggers only, so our records would print
# without a level and a warning could be scrolled past. Give this logger its own
# handler, formatted like uvicorn's, and stop it propagating to avoid duplicates.
if not log.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter("%(levelname)-8s phobius: %(message)s"))
    log.addHandler(_handler)
    log.propagate = False
    log.setLevel(logging.INFO)

_TEMPLATES = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))

#: Bounds how many predictions run at once. The legacy server took a global
#: exclusive lock for every request, so one large upload stalled the whole site.
_slots = threading.BoundedSemaphore(settings.max_concurrency)
_SLOT_TIMEOUT = 30


class Overloaded(RuntimeError):
    """All prediction slots are busy."""


def _with_slot(fn: Callable, *args, **kwargs):
    if not _slots.acquire(timeout=_SLOT_TIMEOUT):
        raise Overloaded(
            "The server is busy with other predictions. Please try again shortly."
        )
    try:
        return fn(*args, **kwargs)
    finally:
        _slots.release()


@asynccontextmanager
async def lifespan(_: FastAPI):
    """Fail fast on a misconfigured deployment rather than on the first request."""
    problems = settings.check()
    for problem in problems:
        log.error("configuration: %s", problem)
    if problems:
        raise RuntimeError("Phobius cannot start:\n  - " + "\n  - ".join(problems))
    homology_ok, homology_reason = homology_search_status(settings)
    if not homology_ok:
        log.info("homology search disabled: %s", homology_reason)

    native_ok, native_reason = native_engine_status(settings)
    if settings.decodeanhmm and not native_ok:
        log.warning(
            "native fast path requested but unusable, falling back to the Java "
            "engine: %s", native_reason,
        )
    log.info(
        "phobius %s ready (model: %s, homology search: %s, native fast path: %s)",
        __version__,
        settings.model,
        "yes" if homology_ok else "no",
        "yes" if native_ok else "no",
    )
    yield


app = FastAPI(
    title="Phobius",
    version=__version__,
    description="Combined transmembrane topology and signal peptide prediction.",
    docs_url="/docs",
    redoc_url=None,
    lifespan=lifespan,
)
app.mount("/static", StaticFiles(directory=str(Path(__file__).parent / "static")), name="static")


# --------------------------------------------------------------------------
# result assembly
# --------------------------------------------------------------------------

@dataclass
class ResultRow:
    name: str
    long: str = ""
    short: str = ""
    svg: str = ""


def _context(request: Request, page: str, **extra) -> dict:
    return {
        "request": request,
        "page": page,
        "settings": settings,
        "short_header": features.SHORT_HEADER,
        "homology_search": settings.homology_search_available(),
        "homology_db": "Swiss-Prot",
        "host": request.url.netloc,
        **extra,
    }


def _render(request: Request, template: str, page: str, status: int = 200, **extra):
    """Render a template. Starlette 1.x wants the request as the first argument."""
    return _TEMPLATES.TemplateResponse(
        request, template, _context(request, page, **extra), status_code=status
    )


def _error(request: Request, title: str, message: str, back: str = "/", status: int = 400):
    return _render(request, "error.html", "error", status=status,
                   title=title, message=message, back=back)


# --------------------------------------------------------------------------
# pages
# --------------------------------------------------------------------------

@app.get("/", response_class=HTMLResponse)
def page_index(request: Request):
    return _render(request, "index.html", "index")


@app.get("/constrained", response_class=HTMLResponse)
def page_constrained(request: Request):
    return _render(request, "constrained.html", "constrained")


@app.get("/poly", response_class=HTMLResponse)
def page_poly(request: Request):
    return _render(request, "poly.html", "poly")


@app.get("/instructions", response_class=HTMLResponse)
def page_instructions(request: Request):
    return _render(request, "instructions.html", "instructions")


def _describe(files: list[Path]) -> list[dict]:
    return [{"name": f.name, "size_mb": round(f.stat().st_size / 1_048_576, 1)}
            for f in files]


@app.get("/download", response_class=HTMLResponse)
def page_download(request: Request):
    files = available_downloads(settings)
    # Partitioned here rather than in the template: Jinja has no regex test, and
    # this keeps size lookups off the render path.
    return _render(
        request, "download.html", "download",
        phobius_bundles=_describe([f for f in files if f.name.startswith("phobius")]),
        other_downloads=_describe([f for f in files if not f.name.startswith("phobius")]),
    )


@app.get("/download/{name}")
def download(name: str):
    """Serve a bundle from the mounted storage.

    The response is chosen from files *discovered* by scanning the directory, so
    the requested name is only ever compared against that set and never joined
    onto a path. That makes traversal impossible rather than merely filtered.
    """
    offered = {f.name: f for f in available_downloads(settings)}
    target = offered.get(name)
    if target is None:
        raise HTTPException(status_code=404, detail=f"No such download: {name}")
    log.info("download: %s", name)
    return FileResponse(target, filename=name, media_type="application/octet-stream")


@app.get("/api", response_class=HTMLResponse)
def page_api(request: Request):
    return _render(request, "api.html", "api")


@app.get("/healthz", response_class=PlainTextResponse)
def healthz():
    """Liveness probe. Reports configuration health without running an engine."""
    problems = settings.check()
    if problems:
        raise HTTPException(status_code=503, detail="; ".join(problems))
    return "ok"


# --------------------------------------------------------------------------
# prediction
# --------------------------------------------------------------------------

async def _read_input(protseq: str, protfile: UploadFile | None) -> str:
    """Return the submitted text. An uploaded file wins over the textarea.

    That precedence matches the legacy form (predict.pl:670).
    """
    if protfile is not None and protfile.filename:
        # Check the declared size before buffering the body into memory.
        if protfile.size is not None and protfile.size > settings.max_upload_bytes:
            raise SubmissionError(
                f"The uploaded file is larger than "
                f"{settings.max_upload_bytes // 1_000_000} MB."
            )
        raw = await protfile.read(settings.max_upload_bytes + 1)
        if len(raw) > settings.max_upload_bytes:
            raise SubmissionError(
                f"The uploaded file is larger than "
                f"{settings.max_upload_bytes // 1_000_000} MB."
            )
        return raw.decode("utf-8", errors="replace")
    return protseq or ""


def _run_prediction(
    records: list[Record],
    output: OutputFormat,
    constraints: list[str] | None = None,
) -> list[ResultRow]:
    """Plain or constrained prediction over a batch of sequences."""
    predictions = engines.predict(records, constraints)
    rows: list[ResultRow] = []

    want_plot = output is OutputFormat.LONG_WITH_PLOT
    posteriors = engines.posteriors(records, constraints) if want_plot else []

    for index, ((name, labels), record) in enumerate(zip(predictions, records)):
        regions = features.label_runs(labels)
        display_name = record.name
        row = ResultRow(name=display_name)
        if output is OutputFormat.SHORT:
            row.short = features.short_format(display_name, regions)
        else:
            row.long = features.long_format(display_name, regions)
        if want_plot:
            row.svg = plot.render(
                posteriors[index],
                regions,
                f"Posterior label probabilities for {display_name}",
            )
        rows.append(row)
    return rows


def _run_poly(records: list[Record], output: OutputFormat) -> list[ResultRow]:
    """PolyPhobius prediction from an alignment.

    The engine reports one label per alignment *column*. Everything the user sees
    is converted back to their own sequence's coordinates, except the
    "graphics only" mode, which exists precisely to show the alignment.
    """
    name, aligned_labels = engines.predict_alignment(records)
    display_name = records[0].name

    if output is OutputFormat.PLOT_ONLY:
        posterior = engines.posteriors_alignment(records)
        return [ResultRow(
            name=display_name,
            svg=plot.render(
                posterior,
                features.label_runs(aligned_labels),
                "Superimposed posterior label probabilities",
                subtitle="alignment coordinates",
            ),
        )]

    regions = features.label_runs(engines.ungap_labels(aligned_labels))
    row = ResultRow(name=display_name)
    if output is OutputFormat.SHORT:
        row.short = features.short_format(display_name, regions)
    else:
        row.long = features.long_format(display_name, regions)

    if output is OutputFormat.LONG_WITH_PLOT:
        posterior = engines.ungap_posterior(
            engines.posteriors_alignment(records), aligned_labels
        )
        row.svg = plot.render(
            posterior,
            regions,
            f"Posterior label probabilities for {display_name}",
            subtitle="homology-supported",
        )
    return [row]


@app.post("/predict", response_class=HTMLResponse)
@app.post("/cgi-bin/predict.pl", response_class=HTMLResponse, include_in_schema=False)
async def predict(
    request: Request,
    protseq: str = Form(default=""),
    protfile: UploadFile | None = File(default=None),
    format: str = Form(default="plp"),
    poly: str = Form(default=""),
    constrained: str = Form(default=""),
    informat: str = Form(default="align"),
    M: str = Form(default=""),
    i: str = Form(default=""),
    o: str = Form(default=""),
    S: str = Form(default=""),
):
    """Handle a form submission.

    Also mounted at the legacy ``/cgi-bin/predict.pl`` path with the original
    field names, so existing bookmarks and scripts keep working.
    """
    is_poly = bool(poly)
    # The legacy form signalled a constrained run by the presence of the 'M'
    # field (predict.pl:639); the new form sends an explicit flag.
    is_constrained = bool(constrained) or (bool(M) and not is_poly)
    back = "/poly" if is_poly else "/constrained" if is_constrained else "/"

    try:
        output = OutputFormat(format)
    except ValueError:
        return _error(request, "Unknown output format",
                      f"'{format}' is not one of short, nog, plp or aplp.", back)

    try:
        text = await _read_input(protseq, protfile)
    except SubmissionError as exc:
        return _error(request, "Submission too large", str(exc), back)

    constraints_used = []
    try:
        if is_poly and informat == InputFormat.ALIGNMENT.value:
            records = parse_alignment(text)
            check_limits(records)
            rows = await run_in_threadpool(_with_slot, _run_poly, records, output)
            heading = "Homology-supported prediction"
            notice = "Prediction is reported for the first sequence in the alignment."

        elif is_poly:
            if output is OutputFormat.PLOT_ONLY:
                return _error(request, "Wrong input format",
                              "'Graphics only' needs aligned FASTA input.", back)
            if not settings.homology_search_available():
                return _error(request, "Homology search unavailable",
                              "This server has no homology database configured. "
                              "Please supply an alignment instead.", back)
            records = parse(text)
            check_limits(records)
            if len(records) != 1:
                return _error(request, "One sequence at a time",
                              "The homology search accepts a single sequence. "
                              f"You submitted {len(records)}.", back)
            aligned = await run_in_threadpool(_with_slot, search_and_align, records[0])
            rows = await run_in_threadpool(_with_slot, _run_poly, aligned, output)
            heading = "Homology-supported prediction"
            notice = (f"Aligned against {len(aligned) - 1} homologues found in "
                      f"Swiss-Prot with BLAST.")

        elif is_constrained:
            if output is OutputFormat.PLOT_ONLY:
                return _error(request, "Wrong output format",
                              "'Graphics only' needs aligned FASTA input.", back)
            records = parse(text)
            check_limits(records)
            records = records[:1]
            constraints_used = build_constraints(
                len(records[0]),
                membrane=M, cytoplasmic=i, non_cytoplasmic=o, signal_peptide=bool(S),
            )
            rows = await run_in_threadpool(
                _with_slot, _run_prediction, records, output,
                constraint_tokens(constraints_used),
            )
            heading = "Constrained prediction"
            notice = "Only the first sequence is used for a constrained prediction."

        else:
            if output is OutputFormat.PLOT_ONLY:
                return _error(request, "Wrong output format",
                              "'Graphics only' needs aligned FASTA input.", back)
            records = parse(text)
            check_limits(records)
            rows = await run_in_threadpool(_with_slot, _run_prediction, records, output)
            heading = "Prediction"
            notice = ""

    except (FastaError, SubmissionError, ConstraintError) as exc:
        return _error(request, "Could not read the submission", str(exc), back)
    except HomologyError as exc:
        return _error(request, "Homology search failed", str(exc), back)
    except engines.EngineError as exc:
        log.warning("engine error: %s", exc)
        return _error(request, "Prediction failed", str(exc), back)
    except Overloaded as exc:
        return _error(request, "Server busy", str(exc), back, status=503)

    log.info("predicted %d sequence(s), format=%s, poly=%s, constrained=%s",
             len(rows), output.value, is_poly, is_constrained)

    return _render(request, "result.html", "result",
                   results=rows, format=output, heading=heading, notice=notice,
                   constraints=constraints_used, back=back)


# --------------------------------------------------------------------------
# JSON / text API
# --------------------------------------------------------------------------

def _predict_structured(text: str, want_plot: bool) -> list[dict]:
    records = parse(text)
    check_limits(records)
    predictions = engines.predict(records)
    posteriors = engines.posteriors(records) if want_plot else []

    out: list[dict] = []
    for index, ((_, labels), record) in enumerate(zip(predictions, records)):
        regions = features.label_runs(labels)
        item = features.to_dict(record.name, record.sequence, regions)
        if want_plot:
            item["svg"] = plot.render(
                posteriors[index], regions,
                f"Posterior label probabilities for {record.name}",
            )
        out.append(item)
    return out


@app.post("/api/predict")
async def api_predict(body: ApiRequest):
    try:
        predictions = await run_in_threadpool(
            _with_slot, _predict_structured, body.sequence, body.plot
        )
    except (FastaError, SubmissionError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except engines.EngineError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    except Overloaded as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc

    return JSONResponse({"predictions": predictions})


@app.post("/api/predict.txt", response_class=PlainTextResponse)
async def api_predict_text(body: ApiRequest, format: str = "nog"):
    def work() -> str:
        records = parse(body.sequence)
        check_limits(records)
        predictions = engines.predict(records)
        lines: list[str] = []
        if format == "short":
            lines.append(features.SHORT_HEADER)
        for (_, labels), record in zip(predictions, records):
            regions = features.label_runs(labels)
            if format == "short":
                lines.append(features.short_format(record.name, regions))
            else:
                lines.append(features.long_format(record.name, regions).rstrip("\n"))
        return "\n".join(lines) + "\n"

    try:
        return await run_in_threadpool(_with_slot, work)
    except (FastaError, SubmissionError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except engines.EngineError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    except Overloaded as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
