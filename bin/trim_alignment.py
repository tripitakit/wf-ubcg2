#!/usr/bin/env python3
"""
Trim overhanging sequences from alignment ends to make all sequences the same length.
Removes gaps at 5' and 3' ends that are not present in all sequences.
"""
import sys
from pathlib import Path
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def trim_alignment(alignment):
    """
    Trim alignment by removing columns that have gaps at the ends.
    Find the leftmost and rightmost positions where all sequences have non-gap characters.
    """
    if len(alignment) == 0:
        return alignment

    aln_length = alignment.get_alignment_length()
    num_seqs = len(alignment)

    # Find first position where all sequences have non-gap characters
    start_pos = 0
    for i in range(aln_length):
        column = alignment[:, i]
        if all(base != '-' for base in column):
            start_pos = i
            break

    # Find last position where all sequences have non-gap characters
    end_pos = aln_length - 1
    for i in range(aln_length - 1, -1, -1):
        column = alignment[:, i]
        if all(base != '-' for base in column):
            end_pos = i
            break

    # Create trimmed alignment
    trimmed_records = []
    for record in alignment:
        trimmed_seq = str(record.seq)[start_pos:end_pos + 1]
        trimmed_record = SeqRecord(
            Seq(trimmed_seq),
            id=record.id,
            description=record.description
        )
        trimmed_records.append(trimmed_record)

    return trimmed_records


def main():
    if len(sys.argv) != 3:
        print("Usage: trim_alignment.py <input_alignment.fasta> <output_trimmed.fasta>", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Read alignment
    alignment = AlignIO.read(input_file, "fasta")
    original_length = alignment.get_alignment_length()

    print(f"Original alignment length: {original_length}", file=sys.stderr)

    # Trim alignment
    trimmed_records = trim_alignment(alignment)

    if trimmed_records:
        trimmed_length = len(trimmed_records[0].seq)
        print(f"Trimmed alignment length: {trimmed_length}", file=sys.stderr)
        print(f"Removed: {original_length - trimmed_length} positions", file=sys.stderr)

        # Write trimmed alignment
        SeqIO.write(trimmed_records, output_file, "fasta")
        print(f"Wrote trimmed alignment to: {output_file}", file=sys.stderr)
    else:
        print("Warning: No sequences in trimmed alignment", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
