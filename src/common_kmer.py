#!/usr/bin/env python
import sys

def read_kmers(file_path):
    kmers = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            kmer = parts[0]
            count = int(parts[1])
            kmers[kmer] = count
    return kmers

def find_common_kmers(kmer_dict_1, kmer_file_2, output_file):
    with open(output_file, 'w') as out:
        with open(kmer_file_2, 'r') as file:
            for line in file:
                parts = line.strip().split()
                kmer = parts[0]
                if kmer in kmer_dict_1:
                    out.write(f'{kmer} {kmer_dict_1[kmer]}\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: script.py <kmer1_dump> <kmer2_dump> <output_file>")
        sys.exit(1)
    
    kmer1_file = sys.argv[1]
    kmer2_file = sys.argv[2]
    output_file = sys.argv[3]

    kmer1_dict = read_kmers(kmer1_file)
    find_common_kmers(kmer1_dict, kmer2_file, output_file)

    print(f"Common kmers have been written to {output_file}")