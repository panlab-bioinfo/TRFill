#!/bin/env python
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import sys

def read_fasta(file_path):
    """Read FASTA file and return the sequence string"""
    with open(file_path) as f:
        seq = ''.join([str(record.seq) for record in SeqIO.parse(f, "fasta")])
    return seq

def clip_sequence(seq, clip_start, clip_end):
    """Clip the specified region of the sequence"""
    if clip_start < 0 or clip_end > len(seq):
        raise ValueError("Clip coordinates exceed sequence length")
    return seq[clip_start:clip_end]

def process_insert(gap_seq, direction):
    """Process the insert sequence based on direction (+/-)"""
    if direction == "+":
        return gap_seq
    elif direction == "-":
        return gap_seq.reverse_complement()
    else:
        raise ValueError("Direction must be '+' or '-'")

def main():
    parser = argparse.ArgumentParser(description="Manually insert a gap sequence into a chromosome")
    parser.add_argument("-c", "--chromosome", required=True, help="Chromosome FASTA file path")
    parser.add_argument("-g", "--gap", required=True, help="Gap FASTA file path to be inserted")
    parser.add_argument("-s", "--insert_start", type=int, required=True, help="Insert start position (chromosome coordinate)")
    parser.add_argument("-e", "--insert_end", type=int, required=True, help="Insert end position (chromosome coordinate)")
    parser.add_argument("-C", "--clip_start", type=int, default=0, help="Clip start position of gap sequence (0-based)")
    parser.add_argument("-E", "--clip_end", type=int, default=0, help="Clip end position of gap sequence (0-based)")
    parser.add_argument("-d", "--direction", choices=["+", "-"], default="+", help="Insertion direction (+ or -)")
    parser.add_argument("-o", "--output", help="Output file path (default: standard output)")

    args = parser.parse_args()

    # Read chromosome sequence
    chrom_seq = read_fasta(args.chromosome)
    chrom_len = len(chrom_seq)
    print(f"[INFO] Chromosome length: {chrom_len} bp", file=sys.stderr)

    # Read gap sequence
    gap_seq = Seq(read_fasta(args.gap))
    gap_len = len(gap_seq)
    print(f"[INFO] Gap sequence length: {gap_len} bp", file=sys.stderr)

    # Clip gap sequence
    if args.clip_end == 0:
        args.clip_end = gap_len
    clipped_gap = clip_sequence(gap_seq, args.clip_start, args.clip_end)
    print(f"[INFO] Clipped gap length: {len(clipped_gap)} bp", file=sys.stderr)

    # Process insertion direction
    processed_gap = process_insert(clipped_gap, args.direction)
    print(f"[INFO] Insertion direction: {args.direction}", file=sys.stderr)

    # Insert into chromosome
    if args.insert_start >= chrom_len or args.insert_end > chrom_len:
        raise ValueError("Insert position exceeds chromosome length")
    if args.insert_start > args.insert_end:
        raise ValueError("Insert start must be <= insert end")

    left_part = chrom_seq[:args.insert_start]
    right_part = chrom_seq[args.insert_end:]
    merged_seq = left_part + str(processed_gap) + right_part

    # Output result
    if args.output:
        with open(args.output, "w") as f:
            f.write(f">merged_chromosome\n{merged_seq}\n")
        print(f"[INFO] Results saved to file: {args.output}", file=sys.stderr)
    else:
        print(f">merged_chromosome\n{merged_seq}")

if __name__ == "__main__":
    main()