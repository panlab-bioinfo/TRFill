#!/usr/bin/env python
import sys
from Bio import SeqIO

def read_config(config_file):
    """
    Reads configuration from the given file.
    :param config_file: Path to the configuration file.
    :return: Dictionary containing configuration parameters.
    """
    config = {}
    with open(config_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            key, value = line.split('=', 1)
            # Handle lists like chrs, starts, ends
            if key in ['chrs', 'starts', 'ends']:
                value = list(map(lambda x: int(x) if x.lstrip('-').isdigit() else x, value.strip('()').split()))
            config[key.strip()] = value.strip() if isinstance(value, str) else value
    return config

def extract_sequences(reference_fa, chrs, starts, ends):
    """
    Extracts sequences from the reference FASTA file based on chromosome names and positions.
    :param reference_fa: Path to the reference FASTA file.
    :param chrs: List of chromosome names.
    :param starts: List of start positions (1-based).
    :param ends: List of end positions (inclusive).
    :return: Dictionary containing extracted sequences.
    """
    sequences = {}
    records = SeqIO.index(reference_fa, "fasta")
    
    for chr_name, start, end in zip(chrs, starts, ends):
        chr_name = chr_name.replace("\"", "")
        if chr_name in records:
            seq = str(records[chr_name].seq[start-1:end])
            sequences[chr_name] = seq
        else:
            print(f"Warning: Chromosome {chr_name} not found in reference.")
    return sequences

def write_sequences_to_fasta(sequences, output_file):
    """
    Writes sequences to a FASTA file.
    :param sequences: Dictionary containing chromosome names and their sequences.
    :param output_file: File to save the output.
    """
    with open(output_file, 'w') as f:
        for chr_name, seq in sequences.items():
            f.write(f">{chr_name}\n{seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: exact_ref.py <reference.fa> <config_file> <output.fa>")
        sys.exit(1)

    reference_fa = sys.argv[1]
    config_file = sys.argv[2]
    output_file = sys.argv[3]

    config = read_config(config_file)
    chrs = config['chrs']
    starts = config['starts']
    ends = config['ends']

    # Extract sequences
    sequences = extract_sequences(reference_fa, chrs, starts, ends)

    # Write sequences to output FASTA file
    write_sequences_to_fasta(sequences, output_file)

    print(f"Sequences slices have been written to {output_file}")