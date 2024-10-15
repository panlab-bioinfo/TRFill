#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def insert_sequence(main_seq_file, insert_seq, main_seq_id, start, end, output_file):
    end_dis = 0
    with open(output_file, "w") as output_handle:
        records = SeqIO.parse(main_seq_file, "fasta")
        for record in records:
            # print(record.id)
            # print(main_seq_id)
            if (record.id == main_seq_id):
                main_seq = str(record.seq)
                new_seq = main_seq[:start] + insert_seq + main_seq[end:]
                record.seq = Seq(new_seq)
                end_dis = start + len(insert_seq) + 1 - end
                SeqIO.write(record, output_handle, "fasta")
                break
    return end_dis

if __name__ == "__main__":
    mat_seq_file = sys.argv[1]
    pat_seq_file = sys.argv[2]
    main_seq_id = sys.argv[3]
    mat_start = int(sys.argv[4])
    mat_end = int(sys.argv[5])
    pat_start = int(sys.argv[6])
    pat_end = int(sys.argv[7])
    mat_output_file = main_seq_id + "_mat_filled.fasta"
    pat_output_file = main_seq_id + "_pat_filled.fasta"

    cen1, cen2 = list(SeqIO.parse("to_be_phased_centromere.fa", "fasta"))
    cen1_seq, cen2_seq = cen1.seq, cen2.seq

    with open('result.log', 'r') as f:
        lines = f.readlines()
        last_two_lines = lines[-2:]
        mat2cen1 = int(last_two_lines[0].split()[-1])
        mat2cen2 = int(last_two_lines[1].split()[-1])
    # print(mat2cen1, mat2cen2)
    if mat2cen1 >= mat2cen2:
        mat_end_dis = insert_sequence(mat_seq_file, cen1_seq, main_seq_id, mat_start, mat_end, mat_output_file)
        pat_end_dis = insert_sequence(pat_seq_file, cen2_seq, main_seq_id, pat_start, pat_end, pat_output_file)
        print(f"{mat_end_dis} {pat_end_dis}")
    else:
        mat_end_dis = insert_sequence(mat_seq_file, cen2_seq, main_seq_id, mat_start, mat_end, mat_output_file)
        pat_end_dis = insert_sequence(pat_seq_file, cen1_seq, main_seq_id, pat_start, pat_end, pat_output_file)
        print(f"{mat_end_dis} {pat_end_dis}")