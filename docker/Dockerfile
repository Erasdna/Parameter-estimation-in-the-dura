FROM ubuntu:22.04

ENV DEB_PYTHON_INSTALL_LAYOUT=deb_system

# Install pip and git with apt
RUN apt-get update && \
    apt-get install -y python3-pip git wget unzip rsync

# We upgrade pip and setuptools
RUN python3 -m pip install --no-cache-dir pip setuptools --upgrade
# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

ARG REPO_BRANCH="main"

WORKDIR /repo

COPY . .

ENV UV_COMPILE_BYTECODE=0

# Make environment
RUN uv sync --verbose

ENV CORES="10"
ENV N_QMC="1000"

ENTRYPOINT [ "/bin/bash", "code/run.sh"]
CMD [ "10", "1000"  ]
