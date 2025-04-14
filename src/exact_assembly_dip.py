#!/usr/bin/env python
import sys
from Bio import SeqIO
import os

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
            if key in ['chrs', 'mat_starts', 'mat_ends', 'pat_starts', 'pat_ends']:
                value = list(map(lambda x: int(x) if x.lstrip('-').isdigit() else x, value.strip('()').split()))
            config[key.strip()] = value.strip() if isinstance(value, str) else value
            if key in ['assembly_pat', 'assembly_pat']:
                value = value.strip()
                config[key.strip()] = value
    return config

def write_sequences_to_fasta(sequences, output_file): 
    """
    Writes sequences to a FASTA file.
    :param sequences: Dictionary containing chromosome names and their sequences.
    :param output_file: File to save the output.
    """
    with open(output_file, 'w') as f:
        for chr_name, seq in sequences.items():
            f.write(f">{chr_name}\n{seq}\n")

def extract_sequences(reference_fa, chrs, starts, ends, output):
    """
    Extracts sequences from the reference FASTA file based on chromosome names and positions.
    :param reference_fa: Path to the reference FASTA file.
    :param chrs: List of chromosome names.
    :param starts: List of start positions (1-based).
    :param ends: List of end positions (inclusive).
    :return: Dictionary containing extracted sequences.
    """
    # chr_visited = {}
    records = SeqIO.index(reference_fa, "fasta")
    print(chrs, starts, ends)
    print(records)
    for chr_name, start, end in zip(chrs, starts, ends):
        # if chr_name not in chr_visited:
        #     chr_visited[chr_name] = 1
        # else:
        #     chr_visited[chr_name] +=1
        chr_name = chr_name.replace("\"", "")
        print(chr_name)
        if chr_name in records:
            sequences = {}
            chr_name_l = chr_name+"_l"
            chr_name_r = chr_name+"_r"
            total_len = len(records[chr_name].seq)
            seq_l = str(records[chr_name].seq[max(start-10000, 0):start])
            seq_r = str(records[chr_name].seq[end: min(end+10000, total_len)])
            sequences[chr_name_l] = seq_l
            sequences[chr_name_r] = seq_r
            write_sequences_to_fasta(sequences, output+"/"+chr_name+".shores.fa")
        else:
            print(f"Warning: Chromosome {chr_name} not found in reference.")
            sys.exit()
        # gap_fa = SeqIO.read("chr_name/scaffolding/hifi_paf_link.available.fa")
        # sequences[chr_name+"_gap"] = str(gap_fa.seq)
    return 1



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: exact_assembly_dip.py <config_file> <output_path>")
        sys.exit(1)

    config_file = sys.argv[1]
    output = sys.argv[2]

    config = read_config(config_file)
    chrs = config['chrs']
    mat_starts = config['mat_starts']
    mat_ends = config['mat_ends']
    pat_starts= config["pat_starts"]
    pat_ends= config["pat_ends"]

    # Extract sequences
    os.mkdir(f"{output}/mat")
    os.mkdir(f"{output}/pat")
    extract_sequences(config["assembly_mat"], chrs, mat_starts, mat_ends, output+'/mat')
    extract_sequences(config["assembly_pat"], chrs, pat_starts, pat_ends, output+"/pat")
    # Write sequences to output FASTA file

    print(f"shores and gap sequence have been written to {output}")