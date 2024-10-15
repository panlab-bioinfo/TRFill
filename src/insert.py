#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def insert_sequence(main_seq_file, insert_seq, main_seq_id, start, end, output_file):
    end_dis = 0
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(main_seq_file, "fasta"):
            if (record.id == main_seq_id):
                main_seq = str(record.seq)
                new_seq = main_seq[:start] + insert_seq + main_seq[end:]
                record.seq = Seq(new_seq)
                end_dis = start + len(insert_seq) + 1 - end
                SeqIO.write(record, output_handle, "fasta")
                break
    return end_dis



if __name__ == "__main__":
    main_seq_file = sys.argv[1]
    insert_seq_file = sys.argv[2]
    main_seq_id = sys.argv[3]
    start = int(sys.argv[4])
    end = int(sys.argv[5])
    output_file = main_seq_id + "_filled.fasta"
    insert_seq_record = SeqIO.read(insert_seq_file, "fasta")
    insert_seq = str(insert_seq_record.seq)
    end_dis = insert_sequence(main_seq_file, insert_seq, main_seq_id, start, end, output_file)
    print(end_dis)
