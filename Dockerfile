FROM patrickdemarta/ubcg2:latest

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    procps \
    wget \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba
RUN wget -qO /tmp/micromamba.tar.bz2 https://micro.mamba.pm/api/micromamba/linux-64/latest && \
    tar -xvj -C /usr/local bin/micromamba -f /tmp/micromamba.tar.bz2 && \
    rm /tmp/micromamba.tar.bz2

# Set up micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/bin:$PATH

# Initialize micromamba and create base environment with all dependencies
RUN micromamba create -y -n base -c conda-forge -c bioconda \
    mafft=7.520 \
    python=3.11 \
    biopython=1.81 \
    && micromamba clean --all --yes

# Activate environment by default
RUN echo "source /usr/local/bin/micromamba shell hook --shell bash && micromamba activate base" >> /etc/bash.bashrc
ENV PATH=/opt/conda/bin:$PATH
ENV PYTHONPATH=/opt/conda/lib/python3.11/site-packages:$PYTHONPATH

# Create activation script for non-interactive shells
RUN micromamba shell init -s bash && \
    echo "micromamba activate base" >> ~/.bashrc

# Verify installations
RUN /opt/conda/bin/python3 --version && \
    /opt/conda/bin/python3 -c "import Bio; print('Biopython version:', Bio.__version__)" && \
    /opt/conda/bin/mafft --version && \
    java -version

# Set python3 as default
RUN ln -sf /opt/conda/bin/python3 /usr/bin/python3 && \
    ln -sf /opt/conda/bin/python3 /usr/bin/python

# Stay as root user for nextflow execution
USER root
