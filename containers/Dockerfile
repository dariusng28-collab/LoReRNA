FROM mambaorg/micromamba:latest

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    python -m pip install --no-cache-dir MisER==1.0.0 && \
    micromamba clean --all --yes

WORKDIR /workspace
