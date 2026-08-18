# Phobius web service.
#
# The image contains no licensed material: the prediction model (phobius.model)
# and the optional native decodeanhmm binary are mounted at runtime from a
# non-public volume. That keeps the image freely publishable, which SciLifeLab
# Serve requires, while the licensed parts stay under their own terms.
FROM python:3.12-slim-bookworm

# openjdk    - runs homologhmm.jar (GPL) and biojava.jar (LGPL)
# diamond    - homology search for PolyPhobius, replacing legacy blastall
# kalign     - multiple alignment, replacing the bundled 32-bit Kalign
# curl       - container healthcheck
RUN apt-get update && apt-get install -y --no-install-recommends \
        openjdk-17-jre-headless \
        diamond-aligner \
        kalign \
        curl \
    && rm -rf /var/lib/apt/lists/*

# SciLifeLab Serve requires the container to run as a non-root user with UID
# 1000. Any other UID needs approval from the platform team.
RUN useradd --uid 1000 --create-home --shell /usr/sbin/nologin phobius

WORKDIR /app

COPY pyproject.toml ./
COPY app ./app
RUN pip install --no-cache-dir --no-compile . && rm -rf /root/.cache

COPY engine/homologhmm.jar engine/biojava.jar ./engine/
COPY scripts ./scripts

# SciLifeLab Serve launches the container by running ./start-script.sh from the
# working directory, so it must sit at WORKDIR and be executable. chmod runs
# here, while we are still root and before USER below.
COPY start-script.sh ./start-script.sh
RUN chmod +x ./start-script.sh

# Mount point for the licensed model and the DIAMOND database.
RUN mkdir -p /mnt/data && chown phobius:phobius /mnt/data

ENV PHOBIUS_MODEL=/mnt/data/phobius.model \
    PHOBIUS_ENGINE_DIR=/app/engine \
    PHOBIUS_DIAMOND_DB=/mnt/data/swissprot.dmnd \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

USER 1000

# Serve only allows ports 3000-9999.
EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD curl -fsS http://127.0.0.1:8000/healthz || exit 1

ENTRYPOINT ["./start-script.sh"]
