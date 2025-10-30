#!/usr/bin/env python3
"""
Parse UCG JSON files and extract gene sequences.
For each gene found across all samples, create a FASTA file.
"""
import json
import sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_ucg_file(ucg_path):
    """Parse a single UCG JSON file and extract sample name and gene sequences."""
    with open(ucg_path, 'r') as f:
        data = json.load(f)

    sample_name = data['genome_info']['label']
    genes = {}
    paralog_info = {}  # Track paralog information

    # Extract DNA sequences for each gene
    for gene_name, gene_data_list in data['data'].items():
        if gene_data_list and len(gene_data_list) > 0:
            # Record paralog information if multiple copies exist
            n_copies = len(gene_data_list)
            if n_copies > 1:
                paralog_info[gene_name] = {
                    'n_copies': n_copies,
                    'selected_bitscore': gene_data_list[0].get('bitscore', 'N/A'),
                    'selected_evalue': gene_data_list[0].get('eValue', 'N/A'),
                    'discarded_copies': n_copies - 1
                }

            # Take the first entry (best match)
            gene_entry = gene_data_list[0]
            if 'dna' in gene_entry:
                genes[gene_name] = gene_entry['dna']

    return sample_name, genes, paralog_info


def write_detailed_report(output_dir, all_genes_per_sample, all_paralogs,
                          genes_in_all, missing_genes, ucg_files):
    """Write a detailed report file with gene and paralog statistics."""
    report_file = output_dir / "gene_selection_report.txt"

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("UCG ALIGNMENT AND CONCATENATION - GENE SELECTION REPORT\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Number of Samples: {len(ucg_files)}\n")
        f.write(f"Samples: {', '.join(sorted(all_genes_per_sample.keys()))}\n\n")

        # Section 1: Gene Statistics
        f.write("-" * 80 + "\n")
        f.write("1. GENE STATISTICS\n")
        f.write("-" * 80 + "\n")
        all_genes = set()
        for genes in all_genes_per_sample.values():
            all_genes.update(genes)
        f.write(f"Total unique genes detected: {len(all_genes)}\n")
        f.write(f"Genes present in ALL samples: {len(genes_in_all)}\n")
        f.write(f"Genes missing in ≥1 sample: {len(missing_genes)}\n\n")

        # Section 2: Missing Genes Detail
        if missing_genes:
            f.write("-" * 80 + "\n")
            f.write("2. MISSING GENES (EXCLUDED FROM CONCATENATION)\n")
            f.write("-" * 80 + "\n")
            for gene in sorted(missing_genes):
                samples_with = [s for s, g in all_genes_per_sample.items() if gene in g]
                samples_without = [s for s, g in all_genes_per_sample.items() if gene not in g]
                f.write(f"\nGene: {gene}\n")
                f.write(f"  Present in {len(samples_with)}/{len(ucg_files)} samples: {', '.join(samples_with)}\n")
                f.write(f"  Missing in {len(samples_without)}/{len(ucg_files)} samples: {', '.join(samples_without)}\n")

        # Section 3: Paralog Information
        f.write("\n" + "-" * 80 + "\n")
        f.write("3. PARALOG DETECTION AND SELECTION\n")
        f.write("-" * 80 + "\n")

        total_paralogs = sum(len(p) for p in all_paralogs.values())
        if total_paralogs > 0:
            f.write(f"Total paralog events detected: {total_paralogs}\n\n")

            for sample, paralog_info in sorted(all_paralogs.items()):
                if paralog_info:
                    f.write(f"\nSample: {sample}\n")
                    f.write(f"  Genes with paralogs: {len(paralog_info)}\n")
                    for gene, info in sorted(paralog_info.items()):
                        f.write(f"    - {gene}:\n")
                        f.write(f"        Total copies detected: {info['n_copies']}\n")
                        f.write(f"        Selected copy (first/best): bitscore={info['selected_bitscore']}, E-value={info['selected_evalue']}\n")
                        f.write(f"        Discarded copies: {info['discarded_copies']}\n")
        else:
            f.write("No paralogs detected (all genes are single-copy).\n")

        # Section 4: Genes Included in Concatenation
        f.write("\n" + "-" * 80 + "\n")
        f.write("4. GENES INCLUDED IN FINAL SUPERMATRIX\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total genes in concatenation: {len(genes_in_all)}\n")
        f.write(f"Genes: {', '.join(sorted(genes_in_all))}\n\n")

        f.write("=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

    print(f"✓ Detailed report written to: {report_file}", file=sys.stderr)


def main():
    if len(sys.argv) < 3:
        print("Usage: parse_ucg.py <output_dir> <ucg_file1> <ucg_file2> ...", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(sys.argv[1])
    output_dir.mkdir(parents=True, exist_ok=True)

    ucg_files = sys.argv[2:]

    # Dictionary to store gene sequences: {gene_name: {sample_name: sequence}}
    gene_sequences = defaultdict(dict)

    # Track all genes across all samples
    all_genes_per_sample = {}  # {sample_name: set(gene_names)}
    all_paralogs = {}  # {sample_name: {gene_name: paralog_info}}

    # Parse all UCG files
    for ucg_file in ucg_files:
        sample_name, genes, paralog_info = parse_ucg_file(ucg_file)
        all_genes_per_sample[sample_name] = set(genes.keys())
        all_paralogs[sample_name] = paralog_info

        # Report paralogs immediately
        if paralog_info:
            print(f"\n⚠️  Sample {sample_name}: {len(paralog_info)} gene(s) with paralogs:", file=sys.stderr)
            for gene, info in sorted(paralog_info.items()):
                print(f"   - {gene}: {info['n_copies']} copies (selected copy: bitscore={info['selected_bitscore']}, E-value={info['selected_evalue']})", file=sys.stderr)

        print(f"✓ Parsed sample: {sample_name} with {len(genes)} genes", file=sys.stderr)

        for gene_name, sequence in genes.items():
            gene_sequences[gene_name][sample_name] = sequence

    # Identify missing genes
    all_gene_names = set()
    for genes in all_genes_per_sample.values():
        all_gene_names.update(genes)

    genes_in_all_samples = set(all_gene_names)
    for sample_genes in all_genes_per_sample.values():
        genes_in_all_samples &= sample_genes

    missing_genes = all_gene_names - genes_in_all_samples

    # Report missing genes
    if missing_genes:
        print(f"\n⚠️  {len(missing_genes)} gene(s) missing in one or more samples (will be EXCLUDED from concatenation):", file=sys.stderr)
        for gene in sorted(missing_genes):
            samples_with_gene = [s for s, genes in all_genes_per_sample.items() if gene in genes]
            samples_without_gene = [s for s, genes in all_genes_per_sample.items() if gene not in genes]
            print(f"   - {gene}:", file=sys.stderr)
            print(f"     Present in: {', '.join(samples_with_gene)}", file=sys.stderr)
            print(f"     Missing in: {', '.join(samples_without_gene)}", file=sys.stderr)

    # Write FASTA file for each gene
    genes_written = 0
    genes_excluded = 0
    for gene_name, samples in gene_sequences.items():
        # Only write genes that are present in all samples
        if len(samples) == len(ucg_files):
            output_file = output_dir / f"{gene_name}.fasta"
            records = []
            for sample_name, sequence in samples.items():
                record = SeqRecord(
                    Seq(sequence),
                    id=sample_name,
                    description=f"{gene_name}"
                )
                records.append(record)

            SeqIO.write(records, output_file, "fasta")
            genes_written += 1
            print(f"Wrote {gene_name} with {len(records)} sequences", file=sys.stderr)
        else:
            genes_excluded += 1

    # Enhanced summary
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"SUMMARY", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Total samples: {len(ucg_files)}", file=sys.stderr)
    print(f"Total unique genes found: {len(all_gene_names)}", file=sys.stderr)
    print(f"Genes in ALL samples: {genes_written} (included in concatenation)", file=sys.stderr)
    print(f"Genes missing in ≥1 sample: {genes_excluded} (excluded from concatenation)", file=sys.stderr)

    # Report total paralogs
    total_paralog_events = sum(len(p) for p in all_paralogs.values())
    if total_paralog_events > 0:
        print(f"Total paralog events: {total_paralog_events}", file=sys.stderr)

    print(f"Output directory: {output_dir}", file=sys.stderr)
    print(f"{'='*60}\n", file=sys.stderr)

    # Write detailed report
    write_detailed_report(output_dir, all_genes_per_sample, all_paralogs,
                         genes_in_all_samples, missing_genes, ucg_files)


if __name__ == "__main__":
    main()
