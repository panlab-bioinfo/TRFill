#!/bin/env python
from Bio import SeqIO, Seq
import sys
import argparse
from collections import defaultdict
import statistics


def process_paf(paf_file):
    records = defaultdict(list)
    with open(paf_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:  # 确保包含基本字段
                continue

            # 解析基础字段
            try:
                mapq = int(parts[11])  # 第12列是MAPQ
            except ValueError:
                continue
            is_primary = False
            for opt in parts[12:]:
                if opt.startswith('tp:A:'):
                    if opt == 'tp:A:P':
                        is_primary = True
                    break
            if mapq < 30:
                continue
            if not is_primary:
                continue
            qstart = int(parts[2])
            qend = int(parts[3])
            strand = parts[4]
            tname = parts[5]
            tstart = int(parts[7])
            tend = int(parts[8])
            mapq = int(parts[11].split(':')[-1])
            records[tname].append({
                'tname' : tname,
                'qstart': qstart,
                'qend': qend,
                'strand': strand,
                'tstart': tstart,
                'tend': tend,
                'mapq': mapq
            })
    return records

def determine_direction(shore, gap_len):
    l_half = gap_len/2
    count = 0
    for i in shore:
        if i["qstart"] < l_half:
            count += 1
    if count >= len(shore)/2:
        ori = "l"
    else:
        ori = "r"
    return ori


def get_start_end(shore_l, shore_r):
    all_qstarts = [x['qstart'] for x in shore_l + shore_r]
    all_qends = [x['qend'] for x in shore_l + shore_r]
    start = min(all_qstarts)
    end = max(all_qends)
    return start, end

def get_optimal_gap_alignment(gap_records, chr_name, gap_len, default_ori):
    """
    gap_records: dict of align in paf
    chr_name: current chromsome
    """
    start_trim = 0
    end_trim = gap_len
    ori = default_ori
    if chr_name+"_l" in gap_records:
        shore_l = gap_records[chr_name+"_l"]
        shore_l.sort(key=lambda x: x['qstart'])
        ori_l = determine_direction(shore_l, gap_len)

    if chr_name+'_r' in gap_records:
        shore_r = gap_records[chr_name+'_r']
        shore_r.sort(key = lambda x: x['qstart'])
        ori_r = determine_direction(shore_r, gap_len)

    
    if (ori_l and ori_r):
        if (ori_l == 'l') and (ori_r == 'r'):
            start_trim = shore_l[0]["qend"]
            end_trim = shore_r[-1]["qstart"]
            ori = "+"
        elif (ori_l == "r") and (ori_r== 'l'):
            start_trim = shore_r[0]['qstart']
            end_trim = shore_l[-1]["qend"]
            ori = "-"
        else:
            all_shore = shore_l + shore_r
            all_shore.sort(key=lambda x: x['qstart'])
            n = len(all_shore)-1
            for i in range(0, len(all_shore)):
                tname1 = all_shore[i]["tname"]
                tname2 = all_shore[n-i]["tname"] 
                if tname1 != tname2:
                    start_trim = all_shore[i]["qend"]
                    end_trim = all_shore[n-i]["qstart"]
                    if tname1 == chr_name+"_l":
                        ori = '+'
                    else:
                        ori = "-"
    elif ori_l:
        if (ori_l == 'l'):
            start_trim = shore_l[0]["qend"]
            ori = '+'
        elif (ori_l == "r"):
            end_trim = shore_l[-1]["qstart"]
            ori = '-'
    elif ori_r:
        if (ori_r == 'l'):
            start_trim = shore_l[0]["qend"]
            ori = '-'
        elif (ori_r == "r"):
            end_trim = shore_l[-1]["qstart"]
            ori = '+'
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
    # output = args.chromosome+".fasta"
    records = process_paf(args.paf)
    gap_seq = SeqIO.read(args.gap_fasta, "fasta")
    start_trim, end_trim, ori = get_optimal_gap_alignment(records, args.chromosome, len(gap_seq.seq), args.hic_ori)
    trimed_seq = gap_seq.seq[start_trim : end_trim]
    if start_trim == 0:
        trimed_seq =Seq.Seq('N'*100 + str(trimed_seq))
    if end_trim == len(gap_seq):
        trimed_seq = Seq.Seq(str(trimed_seq)+ 'N'*100)

    if ori == '-':
        trimed_seq = trimed_seq.reverse_complement()
        print(f"The orientation of gap filled back for {args.chromosome} is: -")
    else:
        print(f"The orientation of gap filled back for {args.chromosome} is: +")
    # print(f"the trimmed coordation of gap for {args.chromosome} is: {start_trim}\t{end_trim}")
    assembly = SeqIO.parse(args.assembly_fasta, "fasta")
    for chr in assembly:
        if chr.id == args.chromosome:
            chr_seq = str(chr.seq)
            chromosome = chr_seq[:args.start] + str(trimed_seq) + chr_seq[args.end:]
            with open(args.out, "w") as file:
                file.write(f">{str(chr.id)}\n")
                file.write(f"{chromosome}\n")
    print(len(trimed_seq))

if __name__ == "__main__":
    main()