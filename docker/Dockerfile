FROM condaforge/miniforge3:26.1.1-3

LABEL org.opencontainers.image.source=https://github.com/andrewgull/MGERT2

RUN mamba install --yes -c conda-forge -c bioconda \
        python=3.11 \
        biopython=1.86 \
        "pandas=3.0.1" \
        "matplotlib=3.10.8" \
        "seaborn=0.13.2" \
        "numpy=2.4.1" \
        "blast>=2.16,<3" \
        "emboss=6.6.0" \
    && mamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
