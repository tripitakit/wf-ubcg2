FROM patrickdemarta/ubcg2:latest

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    procps \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba for dependency management
RUN wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local bin/micromamba

# Install MAFFT, Python, and Biopython using micromamba
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    mafft=7.520 \
    python=3.11 \
    biopython=1.81 \
    && micromamba clean --all --yes

# Update PATH to include micromamba environment
ENV PATH="/root/.local/share/mamba/bin:${PATH}"

# Verify installations
RUN java -version && \
    mafft --version && \
    python3 --version && \
    python3 -c "import Bio; print('Biopython version:', Bio.__version__)"

# Stay as root user for nextflow execution
USER root
