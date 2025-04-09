from Bio import SeqIO
import sys
import os 

def exact_shores(fa, chr, start, end, outfile):
    handle=SeqIO.index(fa)
    
