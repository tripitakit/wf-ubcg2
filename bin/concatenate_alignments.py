#!/usr/bin/env python3
"""
Concatenate multiple gene alignments into a single supermatrix.
Creates a partition file indicating the positions of each gene.
"""
import sys
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_alignment(fasta_file):
    """Read a FASTA alignment file and return sequences as dict."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def main():
    if len(sys.argv) < 3:
        print("Usage: concatenate_alignments.py <output_prefix> <alignment1.fasta> <alignment2.fasta> ...", file=sys.stderr)
        sys.exit(1)

    output_prefix = sys.argv[1]
    alignment_files = [Path(f) for f in sys.argv[2:]]

    # Sort alignment files by gene name for consistent order
    alignment_files = sorted(alignment_files, key=lambda x: x.stem)

    # Dictionary to store concatenated sequences: {sample_name: sequence}
    concatenated_seqs = defaultdict(str)

    # List to store partition information
    partitions = []
    current_pos = 1

    print(f"Concatenating {len(alignment_files)} gene alignments...", file=sys.stderr)

    # Process each alignment file
    for aln_file in alignment_files:
        gene_name = aln_file.stem

        # Read alignment
        sequences = read_alignment(aln_file)

        if not sequences:
            print(f"Warning: No sequences in {aln_file}", file=sys.stderr)
            continue

        # Get alignment length (all sequences should be same length)
        aln_length = len(next(iter(sequences.values())))

        # Verify all sequences have same length
        if not all(len(seq) == aln_length for seq in sequences.values()):
            print(f"Error: Sequences in {aln_file} have different lengths", file=sys.stderr)
            sys.exit(1)

        # Concatenate sequences for each sample
        for sample_name, sequence in sequences.items():
            concatenated_seqs[sample_name] += sequence

        # Record partition
        end_pos = current_pos + aln_length - 1
        partitions.append(f"{gene_name} = {current_pos}-{end_pos}")
        current_pos = end_pos + 1

        print(f"  {gene_name}: {len(sequences)} samples, {aln_length} bp", file=sys.stderr)

    # Write concatenated supermatrix
    output_fasta = f"{output_prefix}.fasta"
    records = []
    for sample_name, sequence in sorted(concatenated_seqs.items()):
        record = SeqRecord(
            Seq(sequence),
            id=sample_name,
            description=f"concatenated_supermatrix length={len(sequence)}"
        )
        records.append(record)

    SeqIO.write(records, output_fasta, "fasta")
    print(f"\nWrote supermatrix to: {output_fasta}", file=sys.stderr)
    print(f"  Samples: {len(records)}", file=sys.stderr)
    print(f"  Total length: {len(records[0].seq)} bp", file=sys.stderr)

    # Write partition file
    partition_file = f"{output_prefix}.partitions.txt"
    with open(partition_file, 'w') as f:
        for partition in partitions:
            f.write(partition + '\n')

    print(f"Wrote partition file to: {partition_file}", file=sys.stderr)
    print(f"  Number of genes: {len(partitions)}", file=sys.stderr)


if __name__ == "__main__":
    main()
