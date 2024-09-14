#!/usr/bin/env python
import sys

chr = sys.argv[1]
mat_start = int(sys.argv[2])
mat_end = int(sys.argv[3])
pat_start = int(sys.argv[4])
pat_end = int(sys.argv[5])
# fa_file = sys.argv[6]

"""flag = 0
with open(fa_file) as file:
    for line in file:
        line = line.strip()
        if(line[1:] == chr):
            flag = 1
        elif(flag):
            print(">shore0001l")
            print(line[:start])
            print(">shore0002l")
            print(line[end + 1:])
            break"""


flag = 0
writer = open(f"mat_shores.fa","w")
with open("/data/wenhuaming/data/HG002/assembly/HG002.mat.chr.fa") as file:
    for line in file:
        line = line.strip()
        if(line[1:] == chr):
            flag = 1
        elif(flag == 1):
            print(">mat000001l",file = writer)
            print(line[:mat_start],file = writer)
            print(">mat000002l",file = writer)
            print(line[mat_end + 1:],file = writer)
            break
writer.close()

flag = 0
writer = open(f"pat_shores.fa","w")
with open("/data/wenhuaming/data/HG002/assembly/HG002.pat.chr.fa") as file:
    for line in file:
        line = line.strip()
        if(line[1:] == chr):
            flag = 1
        elif(flag == 1):
            print(">pat000001l",file = writer)
            print(line[:pat_start],file = writer)
            print(">pat000002l",file = writer)
            print(line[pat_end + 1:],file = writer)
            break
writer.close()