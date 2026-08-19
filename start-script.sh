#!/bin/sh
# Container entry point.
#
# SciLifeLab Serve requires an executable start-up script at the image's WORKDIR
# and invokes it as ./start-script.sh, so this file must stay at
# /app/start-script.sh. See the Dockerfile.
set -eu

# Directories searched for mounted data. Keep in step with DATA_DIRS and
# DATA_SUBDIRS in app/config.py -- the platform chooses where project storage is
# mounted, and each mount is searched both directly and in its db/
# sub-directory so a multi-file BLAST database can be kept tidy.
DATA_DIRS=""
for _mount in ${PHOBIUS_DATA_DIR:-} /mnt/data /home/data /data /app/data; do
    DATA_DIRS="$DATA_DIRS $_mount $_mount/db"
done

# Serve offers no way to set environment variables for a custom app, so any
# setting can instead be placed in a phobius.env file on the mounted storage.
# Only KEY=VALUE lines are honoured; anything else is ignored rather than
# executed, so this is not a way to run arbitrary code from the volume.
load_env_file() {
    file="$1"
    [ -f "$file" ] || return 0
    echo "phobius: loading settings from $file"
    while IFS= read -r line || [ -n "$line" ]; do
        case "$line" in
            ''|'#'*) continue ;;
            [A-Za-z_]*=*) ;;
            *) continue ;;
        esac
        key="${line%%=*}"
        value="${line#*=}"
        case "$value" in
            \"*\") value="${value#\"}"; value="${value%\"}" ;;
            \'*\') value="${value#\'}"; value="${value%\'}" ;;
        esac
        export "$key=$value"
        echo "phobius:   $key=$value"
    done < "$file"
}

for dir in $DATA_DIRS; do
    load_env_file "$dir/phobius.env"
done

# Locate the licensed model. An explicit PHOBIUS_MODEL always wins; otherwise
# the same directories are searched, so no configuration is needed when the
# model simply sits on the mounted volume.
if [ -z "${PHOBIUS_MODEL:-}" ] || [ ! -f "${PHOBIUS_MODEL}" ]; then
    for dir in $DATA_DIRS; do
        if [ -f "$dir/phobius.model" ]; then
            PHOBIUS_MODEL="$dir/phobius.model"
            export PHOBIUS_MODEL
            break
        fi
    done
fi

if [ -z "${PHOBIUS_MODEL:-}" ] || [ ! -f "${PHOBIUS_MODEL}" ]; then
    echo "phobius: phobius.model not found." >&2
    echo "phobius: It is licensed and is deliberately not shipped in this image." >&2
    echo "phobius: Searched:" >&2
    for dir in $DATA_DIRS; do
        if [ -d "$dir" ]; then
            echo "phobius:   $dir/phobius.model  (directory exists, contents: $(ls -A "$dir" 2>/dev/null | tr '\n' ' ' | sed 's/ $//'))" >&2
        else
            echo "phobius:   $dir/phobius.model  (no such directory)" >&2
        fi
    done
    echo "phobius: Mount the model at one of those paths, or set PHOBIUS_MODEL" >&2
    echo "phobius: (or PHOBIUS_DATA_DIR) in a phobius.env file on the same storage." >&2
    exit 1
fi

echo "phobius: using model $PHOBIUS_MODEL"

# Serve allows ports 3000-9999; PORT is honoured so the same image runs
# unchanged on platforms that assign one.
exec uvicorn app.main:app \
    --host 0.0.0.0 \
    --port "${PORT:-8000}" \
    --workers 1 \
    --proxy-headers \
    --forwarded-allow-ips="*"
