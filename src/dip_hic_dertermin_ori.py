#!/usr/bin/env python
"""

Usage:
    python3 script.py <ref_file> <read1_file> <read2_file>
"""

import sys
from Bio import SeqIO
def load_ref_sequences(ref_file):
    ref_dict = {}
    current_seq = None
    kmer_dict = {}
    with open(ref_file, 'r') as f:
        file = f.read().split("\n")
        for line in file:
            line = line.strip()
            if line.startswith('@'):
                current_seq = line[1:]  
                ref_dict[current_seq] = []
            else:
                parts = line.split('\t')
                if len(parts) >= 2:
                    position = int(parts[0], 16)
                    kmer = int(parts[1], 16)
                    kmer_dict[kmer] = (position, current_seq)
    return kmer_dict

def parse_reads(reads_file, kmer_dict, reflen):
    """
    return:
        reads_info: {ref_id, l/r}
    """
    reads_info = {} 
    with open(reads_file, 'r') as f:
        reads = f.read()[1:].split("@")
        for read in reads:
            read = read.strip().split("\n")
            name = read[0].replace("/2", "/1")
            refs = {}
            l, r = 0, 0
            # print(read)
            for kmers in read[1:]:
                # print (kmers)
                pos, kmer, ori = kmers.split("\t")
                
                ref2read = kmer_dict[int(kmer, 16)]
                if ref2read[0] < reflen[ref2read[1]] / 2:
                    l += 1
                else:
                    r += 1
                if ref2read not in refs:
                    refs[ref2read] = 0
                refs[ref2read] += 1
            ref = max(refs.items(), key=lambda item: item[1])
            flag = 'l' if l >= r else 'r'
            reads_info[name] = (ref[0], flag)
    return reads_info


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"{sys.argv[0]} <fasta> <ref_file> <read1_file> <read2_file>")
        sys.exit(1)
    ref_fa = sys.argv[1]    
    ref_file = sys.argv[2]
    read1_file = sys.argv[3]
    read2_file = sys.argv[4]

    reflen = {}
    kmer_dict = load_ref_sequences(ref_file)
    # print(kmer_dict)
    handle = SeqIO.parse(ref_fa, "fasta")

    for rec in handle:
        reflen[rec.id]=len(rec.seq)
    reads1 = parse_reads(read1_file, kmer_dict, reflen)

    reads2 = parse_reads(read2_file, kmer_dict, reflen)

    print(reads1)
    reflink = {
        'cen000001l-mat000001l': [0, 0, 0], 
        'cen000001l-mat000002l': [0, 0, 0],
        'cen000001l-pat000001l': [0, 0, 0],
        'cen000001l-pat000002l': [0, 0, 0],
        'cen000002l-mat000001l': [0, 0, 0], 
        'cen000002l-mat000002l': [0, 0, 0],
        'cen000002l-pat000001l': [0, 0, 0],
        'cen000002l-pat000002l': [0, 0, 0]
    }

    for r in reads1:
        if r not in reads2:
            continue
        ref1 = reads1[r][0][1]
        strand1 = reads1[r][1]
        ref2 = reads2[r][0][1]
        strand1 = reads2[r][1]
        key1 = ref1+'-'+ref2
        key2 = ref2+'-'+ref1
        if key1 in reflink:
            # print(key1)
            if strand1 == 'l':
                reflink[key1][0] += 1
                reflink[key1][1] += 1
            else:
                reflink[key1][0] += 1
                reflink[key1][2] += 1
        
        elif key2 in reflink:
            if reads2[r][1] == 'l':
                reflink[key2][0] += 1
                reflink[key2][1] += 1
            else:
                reflink[key2][0] += 1
                reflink[key2][2] += 1
    del(reads1)
    del(reads2)
    del(kmer_dict)
    cen1_mat_cen2_pat = reflink['cen000001l-mat000001l'][0] + reflink['cen000001l-mat000002l'][0] + \
                        reflink['cen000002l-pat000001l'][0] + reflink['cen000002l-pat000002l'][0]
    cen1_pat_cen2_mat = reflink['cen000001l-pat000001l'][0] + reflink['cen000001l-pat000002l'][0] + \
                        reflink['cen000002l-mat000001l'][0] + reflink['cen000002l-mat000002l'][0]
    
    with open('result.log', "w") as file:
        if cen1_mat_cen2_pat >= cen1_pat_cen2_mat:
            cen1_l = reflink['cen000001l-mat000001l'][1] + reflink['cen000001l-mat000002l'][2]
            cen1_r = reflink['cen000001l-mat000001l'][2] + reflink['cen000001l-mat000002l'][1]
            cen2_l = reflink['cen000002l-pat000001l'][1] + reflink['cen000002l-pat000002l'][2]
            cen2_r = reflink['cen000002l-pat000001l'][2] + reflink['cen000002l-pat000002l'][1]
            if cen1_l >= cen1_r:
                file.write("cen000001l\tmat\t+\n")
            else:
                file.write("cen000001l\tmat\t-\n")
            if cen2_l >= cen2_r:
                file.write("cen000002l\tpat\t+\n")
            else:
                file.write("cen000002l\tpat\t-\n")
        else:
            cen1_l = reflink['cen000001l-pat000001l'][1] + reflink['cen000001l-pat000002l'][2]
            cen1_r = reflink['cen000001l-pat000001l'][2] + reflink['cen000001l-pat000002l'][1]
            cen2_l = reflink['cen000002l-mat000001l'][1] + reflink['cen000002l-mat000002l'][2]
            cen2_r = reflink['cen000002l-mat000001l'][2] + reflink['cen000002l-mat000002l'][1]
            if cen1_l >= cen1_r:
                file.write("cen000001l\tpat\t+\n")
            else:
                file.write("cen000001l\tpat\t-\n")
            if cen2_l >= cen2_r:
                file.write("cen000002l\tmat\t+\n")
            else:
                file.write("cen000002l\tmat\t-\n")
    print("link_id\ttotal_link\tl_lint\tr_link")
    for item in reflink.items():
        print(f"{item[0], item[1]}")
    