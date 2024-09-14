#!/usr/bin/env python
from collections import defaultdict
import copy
import itertools

pairs_matchs = defaultdict(list)
transtable = str.maketrans('ATGC','TACG')
trans_strand = str.maketrans('+-',"-+")

def pair_process(query_name,contig1,info1,contig2,info2):
    query_len1 = info1[0]
    query_start1 = info1[1]
    query_end1 = info1[2]
    strand1 = info1[3]
    ref_len1 = info1[4]
    ref_start1 = info1[5]
    ref_end1 = info1[6]
    query_len2 = info2[0]
    query_start2 = info2[1]
    query_end2 = info2[2]
    strand2 = info2[3]
    ref_len2 = info2[4]
    ref_start2 = info2[5]
    ref_end2 = info2[6]
    if(strand1 != strand2):
        if(strand1 == "-"):
            temp = ref_start1
            ref_start1 = ref_len1 - ref_end1
            ref_end1 = ref_len1 - temp
        else:
            temp = ref_start2
            ref_start2 = ref_len2 - ref_end2
            ref_end2 = ref_len2 - temp
    if(strand1 == strand2 and strand1 == "-"):
        temp = query_start1
        query_start1 = query_len1 - query_end1
        query_end1 = query_len1 - temp
        temp = query_start2
        query_start2 = query_len2 - query_end2
        query_end2 = query_len2 - temp
        strand1 = "+"
        strand2 = "+"
    if(query_start1 - ref_start1 < query_start2 - ref_start2 and query_end1 + ref_len1 - ref_end1 < query_end2 + ref_len2 - ref_end2):
        match = (query_end1 + ref_len1 - ref_end1) - (query_start2 - ref_start2)
        if(contig2 + strand2.translate(trans_strand) + contig1 + strand1.translate(trans_strand) in pairs_matchs):
            pairs_matchs[contig2 + strand2.translate(trans_strand) + contig1 + strand1.translate(trans_strand)].append(match)
        else:
            pairs_matchs[contig1 + strand1 + contig2 + strand2].append(match)
        #print(query_name,contig1,strand1,contig2,strand2,match,1)
    elif(query_start1 - ref_start1 > query_start2 - ref_start2 and query_end1 + ref_len1 - ref_end1 > query_end2 + ref_len2 - ref_end2):
        match = (query_end2 + ref_len2 - ref_end2) - (query_start1 - ref_start1)
        if(contig1 + strand1.translate(trans_strand) + contig2 + strand2.translate(trans_strand) in pairs_matchs):
            pairs_matchs[contig1 + strand1.translate(trans_strand) + contig2 + strand2.translate(trans_strand)].append(match)
        else:
            pairs_matchs[contig2 + strand2 + contig1 + strand1].append(match)
        #print(query_name,contig2,strand2,contig1,strand1,match,2)
    #else:
        #print(query_name,contig1,strand1,contig2,strand2,"contained")

PAF = defaultdict(list)
name = ""
with open("hifitohifi_paf_link.available.paf") as file:
    for line in file:
        line = line.strip().split()
        if(line[0] != name and name != ""):
            # process
            for key1,key2 in itertools.combinations(PAF.keys(),2):
                for value1 in PAF[key1]:
                    for value2 in PAF[key2]:
                        pair_process(name,key1,value1,key2,value2)
            PAF.clear()
        name = line[0]
        PAF[line[5]].append([int(line[1]),int(line[2]),int(line[3]),line[4],int(line[6]),int(line[7]),int(line[8])])
    for key1,key2 in itertools.combinations(PAF.keys(),2):
        for value1 in PAF[key1]:
            for value2 in PAF[key2]:
                pair_process(name,key1,value1,key2,value2)
    PAF.clear()

pairs_match = defaultdict(list)
for key,values in pairs_matchs.items():
    """key1 = key[0:10]
    key2 = key[11:21]
    buck = defaultdict(list)
    for match in values:
        buck[match // 10].append(match)
    tempbuck = copy.copy(buck)
    max_key,_ = max(tempbuck.items(),key = lambda items : len(buck[items[0]] + buck[items[0] + 1] + buck[items[0] - 1]))
    temp = buck[max_key - 1] + buck[max_key] + buck[max_key + 1]
    temp.sort()
    pairs_match[key] = temp"""
    print(key,len(values))