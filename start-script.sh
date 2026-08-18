#!/bin/sh
# Container entry point.
#
# SciLifeLab Serve requires an executable start-up script at the image's
# WORKDIR and invokes it as ./start-script.sh, so this file must stay at
# /app/start-script.sh. See the Dockerfile.
set -eu

MODEL="${PHOBIUS_MODEL:-/mnt/data/phobius.model}"
if [ ! -f "$MODEL" ]; then
    echo "phobius: model file not found at $MODEL" >&2
    echo "phobius: phobius.model is licensed and is deliberately not shipped in" >&2
    echo "phobius: this image. Mount it from the project's file storage and, if" >&2
    echo "phobius: it lives elsewhere, set PHOBIUS_MODEL to its path." >&2
    exit 1
fi

# Serve allows ports 3000-9999; PORT is honoured so the same image runs
# unchanged on platforms that assign one.
exec uvicorn app.main:app \
    --host 0.0.0.0 \
    --port "${PORT:-8000}" \
    --workers 1 \
    --proxy-headers \
    --forwarded-allow-ips="*"
