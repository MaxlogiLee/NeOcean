#!/usr/bin/env python3
"""
reverse_fasta.py
Generate decoy sequences by reversing the input FASTA sequences.
Usage:
    python reverse_fasta.py -i input.fa -o output.fa
"""

from Bio import SeqIO
import argparse


def reverse_fasta_sequences(input_file, output_file):
    with open(input_file, "r") as in_handle:
        with open(output_file, "w") as out_handle:
            for record in SeqIO.parse(in_handle, "fasta"):
                reversed_seq = record.seq[::-1]
                reversed_record = record
                reversed_record.seq = reversed_seq
                SeqIO.write(reversed_record, out_handle, "fasta")


def main():
    parser = argparse.ArgumentParser(
        description="Reverse FASTA sequences to generate decoy database."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input FASTA file"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output decoy FASTA file"
    )
    args = parser.parse_args()

    reverse_fasta_sequences(args.input, args.output)
    print(f"Decoy FASTA written to: {args.output}")


if __name__ == "__main__":
    main()
