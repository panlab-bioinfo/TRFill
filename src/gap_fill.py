#!/bin/env python
from Bio.SeqIO import parse as fasta_parse, read as fasta_read
from Bio.Seq import Seq
import sys
import argparse
from collections import defaultdict


def process_paf(paf_file, gap_len):
    records = defaultdict(list)
    with open(paf_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            try:
                mapq = int(parts[11])
            except ValueError:
                continue
            if mapq < 30:
                continue
            is_primary = False
            for opt in parts[12:]:
                if opt.startswith('tp:A:'):
                    if opt == 'tp:A:P':
                        is_primary = True
                    break
            if not is_primary:
                continue
            qstart = int(parts[2])
            qend = int(parts[3])
            strand = parts[4]
            tname = parts[5]
            tstart = int(parts[7])
            tend = int(parts[8])
            # 规范化负链：翻转qstart/qend到正链等效位置
            if strand == '-':
                qstart, qend = gap_len - qend, gap_len - qstart
                strand = '+'  # 规范化后视为正链
            records[tname].append({
                'tname': tname,
                'qstart': qstart,
                'qend': qend,
                'strand': strand,
                'tstart': tstart,
                'tend': tend,
                'mapq': mapq
            })
    return records

def determine_direction(shore, gap_len):
    l_half = gap_len / 2
    count = sum(1 for i in shore if i["qstart"] < l_half)
    if count >= len(shore) / 2:
        return "l"
    else:
        return "r"

def get_optimal_gap_alignment(gap_records, chr_name, gap_len, default_ori):
    start_trim = 0
    end_trim = gap_len
    ori = default_ori
    shore_l = gap_records.get(chr_name + "_l", [])
    shore_r = gap_records.get(chr_name + "_r", [])
    
    ori_l = determine_direction(shore_l, gap_len) if shore_l else None
    ori_r = determine_direction(shore_r, gap_len) if shore_r else None
    
    if ori_l and ori_r:
        if ori_l == 'l' and ori_r == 'r':
            if shore_l:
                start_trim = max(x['qend'] for x in shore_l)
            if shore_r:
                end_trim = min(x['qstart'] for x in shore_r)
            ori = "+"
        elif ori_l == "r" and ori_r == 'l':
            if shore_r:
                start_trim = max(x['qend'] for x in shore_r)
            if shore_l:
                end_trim = min(x['qstart'] for x in shore_l)
            ori = "-"
        else:
            # 不匹配时，回退到default_ori，并使用所有shore的聚合边界
            all_shore = shore_l + shore_r
            if all_shore:
                start_trim = max(0, max(x['qend'] for x in shore_l) if shore_l else 0)
                end_trim = min(gap_len, min(x['qstart'] for x in shore_r) if shore_r else gap_len)
            # 如果default_ori == '-', 后续会reverse
    elif ori_l:
        if ori_l == 'l':
            start_trim = max(x['qend'] for x in shore_l) if shore_l else 0
            ori = '+'
        elif ori_l == "r":
            end_trim = min(x['qstart'] for x in shore_l) if shore_l else gap_len
            ori = '-'
    elif ori_r:
        if ori_r == 'l':
            start_trim = max(x['qend'] for x in shore_r) if shore_r else 0
            ori = '-'
        elif ori_r == "r":
            end_trim = min(x['qstart'] for x in shore_r) if shore_r else gap_len
            ori = '+'
    
    # 确保start_trim < end_trim
    if start_trim >= end_trim:
        print(f"Warning: Invalid trim for {chr_name}, using full sequence with default ori.")
        start_trim = 0
        end_trim = gap_len
        ori = default_ori
    
    return start_trim, end_trim, ori

def main():
    parser = argparse.ArgumentParser(description="Fill chromosome gaps with Hi-C support")
    parser.add_argument("paf", help="gap mapping to shores")
    parser.add_argument("assembly_fasta", help="Path to reference FASTA file")
    parser.add_argument("gap_fasta", help="Path to gap assembly FASTA file")
    parser.add_argument("chromosome", help="Chromosome name")
    parser.add_argument("start", type=int, help="Gap start position (1-based)")
    parser.add_argument("end", type=int, help="Gap end position (1-based)")
    parser.add_argument("out", help="Output file name")
    parser.add_argument("hic_ori", help="Orientation of HiC")
    args = parser.parse_args()
    
    gap_seq = fasta_read(args.gap_fasta, "fasta")
    gap_len = len(gap_seq.seq)
    records = process_paf(args.paf, gap_len)
    start_trim, end_trim, ori = get_optimal_gap_alignment(records, args.chromosome, gap_len, args.hic_ori)
    
    trimed_seq = gap_seq.seq[start_trim:end_trim]
    if start_trim == 0:
        trimed_seq = Seq('N' * 100 + str(trimed_seq))
    if end_trim == gap_len:
        trimed_seq = Seq(str(trimed_seq) + 'N' * 100)
    
    if ori == '-':
        trimed_seq = trimed_seq.reverse_complement()
        print(f"The orientation of gap filled back for {args.chromosome} is: -")
    else:
        print(f"The orientation of gap filled back for {args.chromosome} is: +")
    
    print(f"Trimmed length: {len(trimed_seq)}")
    
    found = False
    for chr_record in fasta_parse(args.assembly_fasta, "fasta"):
        if chr_record.id == args.chromosome:
            chr_seq = str(chr_record.seq)
            # 修正切片：替换1-based start到end (inclusive)
            chromosome = chr_seq[:args.start - 1] + str(trimed_seq) + chr_seq[args.end:]
            with open(args.out, "w") as file:
                file.write(f">{chr_record.id}\n")
                file.write(f"{chromosome}\n")
            found = True
            break
    if not found:
        print(f"Error: Chromosome {args.chromosome} not found in assembly_fasta.", file=sys.stderr)

if __name__ == "__main__":
    main()