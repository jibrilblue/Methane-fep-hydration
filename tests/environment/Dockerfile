FROM condaforge/miniforge3:24.3.0-0

ENV DEBIAN_FRONTEND=noninteractive \
    OPENMM_CPU_THREADS=4

RUN mamba install -y -n base -c conda-forge \
        python=3.10 \
        openmm=8.1.1 \
        openmmtools=0.24.0 \
        pymbar=4.0.3 \
        numpy=1.26.4 \
    && mamba clean -afy

WORKDIR /app
COPY ./inputs/ /app/inputs/
RUN mkdir -p /app/output
