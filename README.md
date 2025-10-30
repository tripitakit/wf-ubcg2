# wf-ubcg2: Bacterial Phylogenomics Workflow

A comprehensive Nextflow workflow for bacterial phylogenomic analysis using Universal Bacterial Core Genes (UBCG2).

## Overview

This workflow combines two independent stages for bacterial phylogenomic analysis:

- **Stage 1: UCG Extraction** - Extract universal bacterial core genes from genome assemblies using UBCG2
- **Stage 2: Alignment & Concatenation** - Align and concatenate core genes from multiple samples to create a phylogenetic supermatrix

Both stages can be run independently or together, providing flexibility for different analysis scenarios.

## Features

- ✅ Extract 81 universal bacterial core genes using UBCG2
- ✅ Multiple sequence alignment with MAFFT
- ✅ Automatic trimming of alignment overhangs
- ✅ Concatenation into phylogenetic supermatrix
- ✅ Partition file generation for phylogenetic software
- ✅ Paralog detection and reporting
- ✅ Docker/Singularity container support
- ✅ EPI2ME Labs compatible

## Quick Start

### Prerequisites

- Nextflow >= 23.04.2
- Docker or Singularity
- Bacterial genome assembly (FASTA format) for Stage 1
- UCG files from multiple samples for Stage 2

### Installation

```bash
# Clone the repository
git clone https://github.com/microbion/wf-ubcg2.git
cd wf-ubcg2

# Build Docker container (optional if using pre-built image)
docker build -t patrickdemarta/wf-ubcg2:v1.0.0 .
```

### Usage Examples

#### Stage 1 Only: Extract UCG from a genome

```bash
nextflow run main.nf \
    --fasta genome.fasta \
    --out_dir output_stage1
```

**Output:**
- `output_stage1/ucg_output/[sample]/*.ucg` - UCG JSON files

#### Stage 2 Only: Create phylogenetic supermatrix

```bash
nextflow run main.nf \
    --ucg_dir path/to/ucg_files/ \
    --mafft_threads 8 \
    --out_dir output_stage2
```

**Output:**
- `output_stage2/concatenated_supermatrix.fasta` - Phylogenetic supermatrix
- `output_stage2/concatenated_supermatrix.partitions.txt` - Partition file for RAxML/IQ-TREE
- `output_stage2/alignments/*.aln` - Individual gene alignments
- `output_stage2/trimmed_alignments/*.trimmed.fasta` - Trimmed alignments

#### Both Stages: Complete phylogenomic pipeline

```bash
# First, extract UCG from multiple genomes
for genome in genomes/*.fasta; do
    nextflow run main.nf \
        --fasta $genome \
        --out_dir output/stage1
done

# Then, create supermatrix from all UCG files
nextflow run main.nf \
    --ucg_dir output/stage1/ucg_output \
    --mafft_threads 8 \
    --out_dir output/stage2
```

## Parameters

### Stage 1: UCG Extraction

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--fasta` | file | `null` | Path to bacterial genome assembly (FASTA format) |

### Stage 2: Alignment & Concatenation

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--ucg_dir` | directory | `null` | Directory containing UCG files (*.ucg) |
| `--mafft_threads` | integer | `4` | Number of threads for MAFFT alignment |

### Common Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--out_dir` | directory | `output` | Output directory for results |
| `--disable_ping` | boolean | `false` | Disable telemetry ping |

## Workflow Details

### Stage 1: UCG Extraction

1. **Input:** Single bacterial genome assembly (FASTA)
2. **Process:** UBCG2 extracts 81 universal bacterial core genes
3. **Output:** UCG JSON file containing gene sequences and metadata

**Key Process:**
- `extractUCG` - Runs UBCG2 Java tool
  - Resources: 4 CPUs, 8GB RAM
  - Runtime: 5-10 minutes per genome

### Stage 2: Alignment & Concatenation

1. **Parse UCG files** - Extract gene sequences from all samples
2. **Align genes** - MAFFT multiple sequence alignment per gene
3. **Trim alignments** - Remove overhanging sequences
4. **Concatenate** - Merge all genes into supermatrix

**Key Processes:**
- `parseUCG` - Extract sequences, detect paralogs
- `alignGene` - MAFFT alignment (parallelized per gene)
- `trimAlignment` - Remove overhangs
- `concatenateAlignments` - Create supermatrix + partition file

## Output Files

### Stage 1 Outputs

```
output/
├── ucg_output/
│   └── [sample]/
│       ├── *.ucg                    # UCG JSON file
│       └── ubcg2.log               # UBCG2 log file
└── versions_extraction.txt          # Software versions
```

### Stage 2 Outputs

```
output/
├── parsed_genes/
│   ├── *.fasta                      # Individual gene FASTA files
│   └── gene_selection_report.txt   # Paralog information
├── alignments/
│   └── *.aln                        # MAFFT alignments
├── trimmed_alignments/
│   └── *.trimmed.fasta             # Trimmed alignments
├── concatenated_supermatrix.fasta   # Final supermatrix
├── concatenated_supermatrix.partitions.txt  # Partition file
└── versions_alignment.txt           # Software versions
```

## Using the Supermatrix for Phylogenetic Analysis

The output supermatrix can be used with various phylogenetic inference tools:

### RAxML-NG

```bash
raxml-ng --all \
    --msa concatenated_supermatrix.fasta \
    --model concatenated_supermatrix.partitions.txt \
    --prefix phylogeny \
    --threads 8 \
    --bs-trees 100
```

### IQ-TREE

```bash
iqtree -s concatenated_supermatrix.fasta \
    -p concatenated_supermatrix.partitions.txt \
    -m MFP \
    -bb 1000 \
    -nt 8
```

### MrBayes

Convert partition file format and use with MrBayes for Bayesian inference.

## Container Information

The workflow uses a custom Docker container that extends `patrickdemarta/ubcg2:latest` with additional tools:

- **UBCG2 v2.0** - Core gene extraction
- **MAFFT v7.520** - Multiple sequence alignment
- **Python 3.11** - Data processing
- **Biopython 1.81** - Sequence manipulation

### Building the Container

```bash
docker build -t patrickdemarta/wf-ubcg2:v1.0.0 .
docker push patrickdemarta/wf-ubcg2:v1.0.0
```

## System Requirements

### Recommended

- **CPUs:** 4
- **Memory:** 8GB
- **Storage:** 10GB for typical analysis

### Minimum

- **CPUs:** 2
- **Memory:** 4GB

## Runtime

- **Stage 1:** 5-10 minutes per genome
- **Stage 2:** 10-30 minutes depending on sample count and gene count

## Execution Profiles

The workflow supports multiple execution profiles:

```bash
# Docker (default)
nextflow run main.nf --fasta genome.fasta

# Singularity
nextflow run main.nf --fasta genome.fasta -profile singularity

# AWS Batch
nextflow run main.nf --fasta genome.fasta -profile awsbatch \
    --aws_queue my-queue \
    --aws_image_prefix my-ecr-prefix

# Local execution (development)
nextflow run main.nf --fasta genome.fasta -profile local
```

## Troubleshooting

### No UCG files found

**Error:** `No UCG files (*.ucg) found in <directory>`

**Solution:** Ensure the directory contains files with `.ucg` extension from Stage 1 output.

### Out of memory error

**Solution:** Increase memory allocation in `nextflow.config` or reduce `mafft_threads`.

### Paralog warnings

**Info:** The workflow detects genes with multiple copies (paralogs) and reports them in `gene_selection_report.txt`. By default, the first copy is used for alignment.

## Citation

If you use this workflow, please cite:

- **UBCG2:** Na SI et al. (2018) UBCG: Up-to-date bacterial core gene set and pipeline for phylogenomic tree reconstruction. Journal of Microbiology 56:280-285. [DOI: 10.1007/s12275-018-8014-6](https://doi.org/10.1007/s12275-018-8014-6)

- **MAFFT:** Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7. Molecular Biology and Evolution 30:772-780. [DOI: 10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

## License

This workflow is distributed under the MIT License. See LICENSE file for details.

## Support

For issues, questions, or contributions:
- **Issues:** [GitHub Issues](https://github.com/microbion/wf-ubcg2/issues)
- **Contact:** Patrick De Marta - Microbion Srl

## Acknowledgments

This workflow combines and extends:
- Original UBCG2 tool by Na et al.
- EPI2ME Labs workflow framework
- Nextflow DSL2 best practices

## Version History

- **v1.0.0** (2025-10-30)
  - Initial release
  - Merged wf-ucg and wf-ucg-aln-concat workflows
  - Two independent stages with flexible execution
  - Custom Docker container with all dependencies
