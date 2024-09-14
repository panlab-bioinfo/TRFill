#!/usr/bin/env python
from collections import defaultdict
import itertools
import sys
import copy
import random
import numpy as np
import networkx as nx

chr = sys.argv[1]
gap_start = int(sys.argv[2])
gap_end = int(sys.argv[3])
combination_file = sys.argv[4]
gfa_file = sys.argv[5]
available_fa_file = sys.argv[6]
mat_pat_available_fa_file = sys.argv[7]
hic_link_file = sys.argv[8]
hifi_link_file = sys.argv[9]
transtable = str.maketrans('ATGC','TACG')
trans_strand = str.maketrans('+-',"-+")
scaffold_scale = 4700
scaffold_penalty_scale = 12
len_scale = 50
haploid_penalty = 100

def calculate_score(phase):
    mat = []
    pat = []
    for i,num in enumerate(phase):
        if(num == 1):
            mat.append(index_reverse[i])
        elif(num == 2):
            pat.append(index_reverse[i])
        elif(num == 3):
            mat.append(index_reverse[i])
            pat.append(index_reverse[i])

    mat = sorted(mat,key = lambda scaffold : (scaffolds_info[scaffold][1] + scaffolds_info[scaffold][2]) / 2)
    pat = sorted(pat,key = lambda scaffold : (scaffolds_info[scaffold][1] + scaffolds_info[scaffold][2]) / 2)
    # if connected components are contained, reorder them in different manners
    for target in set(isolate_node) & set(mat):
        min_index1 = scaffolds_info[target][1]
        max_index1 = scaffolds_info[target][2]
        temp = [set(component) & set(mat) for component in connected_components]
        for component in temp:
            if(len(component) > 1):
                min_index2 = min([scaffolds_info[scaffold][1] for scaffold in component])
                max_index2 = max([scaffolds_info[scaffold][2] for scaffold in component])
                if(min_index2 < min_index1 < max_index1 < max_index2): # mid
                    mat.remove(target)
                    min_insert_index = min([mat.index(scaffold) for scaffold in component])
                    max_insert_index = max([mat.index(scaffold) for scaffold in component]) + 1
                    if(min_index1 + max_index1 < min_index2 + max_index2): 
                        mat.insert(min_insert_index,target)
                    else:
                        mat.insert(max_insert_index,target)
                elif(min_index1 < min_index2 < max_index2 < max_index1): # head
                    mat.remove(target)
                    min_insert_index = min([mat.index(scaffold) for scaffold in component])
                    mat.insert(min_insert_index,target)
    for target in set(isolate_node) & set(pat):
        min_index1 = scaffolds_info[target][1]
        max_index1 = scaffolds_info[target][2]
        temp = [set(component) & set(pat) for component in connected_components]
        for component in temp:
            if(len(component) > 1):
                min_index2 = min([scaffolds_info[scaffold][1] for scaffold in component])
                max_index2 = max([scaffolds_info[scaffold][2] for scaffold in component])
                if(min_index2 < min_index1 < max_index1 < max_index2): # mid
                    pat.remove(target)
                    min_insert_index = min([pat.index(scaffold) for scaffold in component])
                    max_insert_index = max([pat.index(scaffold) for scaffold in component]) + 1
                    if(min_index1 + max_index1 < min_index2 + max_index2): 
                        pat.insert(min_insert_index,target)
                    else:
                        pat.insert(max_insert_index,target)
                elif(min_index1 < min_index2 < max_index2 < max_index1): # head
                    pat.remove(target)
                    min_insert_index = min([pat.index(scaffold) for scaffold in component])
                    pat.insert(min_insert_index,target)

    shore_score = 0
    for shore in shores_order:
        hap = mat if(shore[:3] == "mat") else pat
        if(len(hap) > 0):
            scaffold = hap[0] if(shore[-2] == "1") else hap[-1]
            shore_score += shore_link[shore][scaffold][2]

    scaffold_score = 0
    """for i,scaffold1 in enumerate(mat):
        if(i != len(mat) - 1):
            for scaffold2 in mat[i + 1 : i + 2]: #  min(i + 3,len(mat))
                scaffold_score += scaffold_link[scaffold1][scaffold2][2]
    for i,scaffold1 in enumerate(pat):
        if(i != len(pat) - 1):
            for scaffold2 in pat[i + 1 : i + 2]: #  min(i + 3,len(pat))
                scaffold_score += scaffold_link[scaffold1][scaffold2][2]"""
    for hap in [mat,pat]: # one step window between components, full length window within components
        temp = [set(component) & set(hap) for component in connected_components_plus if(set(component) & set(hap))]
        for component in temp:
            if(len(component) > 1):
                for scaffold1,scaffold2 in itertools.combinations(component,2):
                    scaffold_score += scaffold_link[scaffold1][scaffold2][2]
        for i,scaffold1 in enumerate(hap):
            if(i != len(hap) - 1):
                scaffold2 = hap[i + 1]
                if(scaffold1 not in connected_components_index or scaffold2 not in connected_components_index or connected_components_index[scaffold1] != connected_components_index[scaffold2]):
                    scaffold_score += scaffold_link[scaffold1][scaffold2][2]
    for scaffold1,scaffold2 in itertools.product(*(mat,pat)):
        if(any(scaffold1 in component and scaffold2 in component and (connected_components_len[i] < centromere_len or len(component) <= 4) for i,component in enumerate(connected_components))): # in haploid component, try to keep scaffolds together
            scaffold_score -= haploid_penalty

    len_score = abs(sum([length[scaffold] for scaffold in mat]) - sum([length[scaffold] for scaffold in pat])) / centromere_len
    len_score = len_scale * len_score if(len_score > 0.4) else 0 # penalty if length differ too much
    for hap in [mat,pat]:
        temp = [set(component) & set(hap) for component in connected_components_plus if(set(component) & set(hap))]
        for component1,component2 in itertools.combinations(temp,2):
            min_index1 = min([scaffolds_info[scaffold][1] for scaffold in component1])
            max_index1 = max([scaffolds_info[scaffold][2] for scaffold in component1])
            min_index2 = min([scaffolds_info[scaffold][1] for scaffold in component2])
            max_index2 = max([scaffolds_info[scaffold][2] for scaffold in component2])
            min_index = min(min_index1,min_index2)
            max_index = max(max_index1,max_index2)
            overlap = max_index - min_index - (max_index1 - min_index1 + 1) - (max_index2 - min_index2 + 1)
            if(overlap < 0 and abs(overlap) / max(max_index1 - min_index1 + 1,max_index2 - min_index2 + 1) > 0.5): # different components in one hap should not overlap too much
                len_score += 10 * len_scale * abs(overlap) / max(max_index1 - min_index1 + 1,max_index2 - min_index2 + 1)

    return mat,pat,shore_score,scaffold_score,len_score,shore_score + scaffold_score - len_score # 

centromere_len = gap_end - gap_start + 1

GFA = nx.DiGraph()
with open(gfa_file) as file:
    for line in file:
        line = line.strip().split()
        if(line[0] == "S"):
            GFA.add_node(line[1],length = len(line[2]),seq = line[2])
        elif(line[0] == "L" and line[1] != line[3]):
            GFA.add_edge(line[1],line[3],strand1 = line[2],strand2 = line[4],match = int(line[5][:-1]))

goal_combination = []
with open(combination_file) as file:
    for line in file:
        line = line.strip().split()
        goal_combination.append(line)

scaffolds_info = {}
for interval in goal_combination:
    scaffolds_info[interval[0]] = [interval[1],int(interval[2]),int(interval[3])]
scaffolds_info["mat000001l"] = ["+",gap_start,gap_start]
scaffolds_info["mat000002l"] = ["+",gap_end,gap_end]
scaffolds_info["pat000001l"] = ["+",gap_start,gap_start]
scaffolds_info["pat000002l"] = ["+",gap_end,gap_end]

available_scaffolds = set()
kmers = defaultdict(set)
homologous = {}
flag = 0
name = ""
with open(available_fa_file) as file:
    for line in file:
        line = line.strip()
        if(line[1:] in scaffolds_info):
            name = line[1:]
            available_scaffolds.add(name)
            if(scaffolds_info[name][0] == "+"):
                flag = 1
            else:
                flag = 2
        elif(flag == 1):
            temp = set()
            for i in range(0, len(line) - 21 + 1):
                temp.add(line[i : i + 21])
            kmers[name] = temp
            flag = 0
        elif(flag == 2):
            line = line.translate(transtable)[::-1]
            temp = set()
            for i in range(0, len(line) - 21 + 1):
                temp.add(line[i : i + 21])
            kmers[name] = temp
            flag = 0
for kmer1,kmer2 in itertools.combinations(kmers.items(),2):
    # match = 0
    if((kmer1[0],kmer2[0]) in GFA.edges()):
        # match = GFA[kmer1[0]][kmer2[0]]["match"]
        homologous[" ".join((sorted([kmer1[0],kmer2[0]])))] = 0
    else:
        homologous[" ".join((sorted([kmer1[0],kmer2[0]])))] = (len(kmer1[1] & kmer2[1]) / min(len(kmer1[1]),len(kmer2[1]))) # * (1 - match / min(GFA.nodes[kmer1[0]]["length"],GFA.nodes[kmer2[0]]["length"]))
    print(kmer1[0],kmer2[0])
    print(homologous[" ".join((sorted([kmer1[0],kmer2[0]])))])
homologous_coefficient = sum(homologous.values()) / len([value for value in homologous.values() if(value != 0)])
# homologous_coefficient = sum([(GFA.nodes[key.split()[0]]["length"] * GFA.nodes[key.split()[1]]["length"]) * value for key,value in homologous.items() if(value != 0)]) / sum([(GFA.nodes[key.split()[0]]["length"] * GFA.nodes[key.split()[1]]["length"]) for key,value in homologous.items() if(value != 0)])
print("homologous_coefficient:\t",homologous_coefficient)

length = {}
name = ""
with open(mat_pat_available_fa_file) as file:
    for line in file:
        line = line.strip()
        if(len(line) > 0 and line[0] == ">"):
            name = line[1:]
        else:
            length[name] = len(line)

sources = set()
for edge in GFA.edges(): # get those whose indegree is 0
    src = edge[0]
    dst = edge[1]
    flag = 0
    for node in GFA[dst].keys():
        if(GFA[src][dst]["strand2"] == GFA[dst][node]["strand1"]):
            flag = 1
            break
    if(flag == 0):
        sources.add(dst)
connected_components = []
for node in sources:
    if(any(node in component for component in connected_components)):
        continue
    temp = nx.bfs_tree(GFA,node).nodes()
    connected_components.append(list(temp))
for component in connected_components:
    for scaffold in list(component):
        if(scaffold not in available_scaffolds):
            component.remove(scaffold)
connected_components_len = [sum([length[scaffold] for scaffold in component]) for component in connected_components]
connected_components_index = {}
for scaffold in available_scaffolds:
    for i,component in enumerate(connected_components):
        if(scaffold in component):
            connected_components_index[scaffold] = i

connected_components_plus = copy.copy(connected_components)
isolate_node = []
for node in GFA.nodes():
    if(GFA.in_degree[node] == 0 and GFA.out_degree[node] == 0):
        isolate_node.append(node)
        connected_components_plus.append([node])

avaliable_length_sum = sum([length[scaffold] for scaffold in available_scaffolds])
# if diploid component is large enough, we should consider gain more weight for shore score
if(any(connected_components_len[i] >= centromere_len and len(component) > 4 and connected_components_len[i] / avaliable_length_sum > 0.8 for i,component in enumerate(connected_components))):
    scaffold_scale /= scaffold_penalty_scale

shore_link = defaultdict(dict)
with open(hic_link_file) as file:
    for line in file:
        line = line.strip().split()
        scaffold1,strand1,scaffold2,strand2 = "","","",""
        if(line[0] not in scaffolds_info or line[1] not in scaffolds_info):
            continue
        # fit the orientation
        if((scaffolds_info[line[0]][1] + scaffolds_info[line[0]][2]) / 2 <= (scaffolds_info[line[1]][1] + scaffolds_info[line[1]][2]) / 2):
            scaffold1 = line[0]
            strand1 = scaffolds_info[scaffold1][0]
            scaffold2 = line[1]
            strand2 = scaffolds_info[scaffold2][0]
        else:
            scaffold1 = line[1]
            strand1 = scaffolds_info[scaffold1][0]
            scaffold2 = line[0]
            strand2 = scaffolds_info[scaffold2][0]

        if((("mat" not in scaffold1 and "pat" not in scaffold1 and ("mat" in scaffold2 or "pat" in scaffold2)) # scaffold,shore
            or
            ("mat" not in scaffold2 and "pat" not in scaffold2 and ("mat" in scaffold1 or "pat" in scaffold1))) # shore,scaffold
        and ((strand1 == strand2 and line[2] == "same") or (strand1 != strand2 and line[2] == "diff"))): # fit the order
            if(scaffold1 in ["mat000001l","pat000001l"] or scaffold2 in ["mat000001l","pat000001l"]):
                if(1 and scaffolds_info[scaffold1][1] + 100000 < gap_start or scaffolds_info[scaffold2][1] + 100000 < gap_start): # overflow
                    shore_link[scaffold1][scaffold2] = [strand1,strand2,int(line[3]) / (length[scaffold1] * length[scaffold2] * abs((scaffolds_info[scaffold1][1] + scaffolds_info[scaffold1][2]) / 2 - (scaffolds_info[scaffold2][1] + scaffolds_info[scaffold2][2]) / 2))]
                    shore_link[scaffold2][scaffold1] = [strand2.translate(trans_strand),strand1.translate(trans_strand),int(line[3]) / (length[scaffold1] * length[scaffold2] * abs((scaffolds_info[scaffold1][1] + scaffolds_info[scaffold1][2]) / 2 - (scaffolds_info[scaffold2][1] + scaffolds_info[scaffold2][2]) / 2))]
                else:
                    shore_link[scaffold1][scaffold2] = [strand1,strand2,int(line[3]) / (length[scaffold1] * length[scaffold2] * abs(scaffolds_info[scaffold1][1] - scaffolds_info[scaffold2][1]))]
                    shore_link[scaffold2][scaffold1] = [strand2.translate(trans_strand),strand1.translate(trans_strand),int(line[3]) / (length[scaffold1] * length[scaffold2] * abs(scaffolds_info[scaffold1][1] - scaffolds_info[scaffold2][1]))]
            elif(scaffold1 in ["mat000002l","pat000002l"] or scaffold2 in ["mat000002l","pat000002l"]):
                if(1 and scaffolds_info[scaffold1][2] > gap_end + 100000 or scaffolds_info[scaffold2][2] > gap_end + 100000): # overflow
                    shore_link[scaffold1][scaffold2] = [strand1,strand2,int(line[3]) / (length[scaffold1] * length[scaffold2] * abs((scaffolds_info[scaffold1][1] + scaffolds_info[scaffold1][2]) / 2 - (scaffolds_info[scaffold2][1] + scaffolds_info[scaffold2][2]) / 2))]
                    shore_link[scaffold2][scaffold1] = [strand2.translate(trans_strand),strand1.translate(trans_strand),int(line[3]) / (length[scaffold1] * length[scaffold2] * abs((scaffolds_info[scaffold1][1] + scaffolds_info[scaffold1][2]) / 2 - (scaffolds_info[scaffold2][1] + scaffolds_info[scaffold2][2]) / 2))]
                else:
                    shore_link[scaffold1][scaffold2] = [strand1,strand2,int(line[3]) / (length[scaffold1] * length[scaffold2] * abs(scaffolds_info[scaffold1][2] - scaffolds_info[scaffold2][2]))]
                    shore_link[scaffold2][scaffold1] = [strand2.translate(trans_strand),strand1.translate(trans_strand),int(line[3]) / (length[scaffold1] * length[scaffold2] * abs(scaffolds_info[scaffold1][2] - scaffolds_info[scaffold2][2]))]
print(shore_link)

shores_order = ["mat000001l","mat000002l","pat000001l","pat000002l"]

for shore in shores_order:
    for scaffold in length.keys():
        if("mat" not in scaffold and "pat" not in scaffold and scaffold in scaffolds_info and scaffold not in shore_link[shore]):
            if((scaffolds_info[scaffold][1] + scaffolds_info[scaffold][2]) / 2 >= (scaffolds_info[shore][1] + scaffolds_info[shore][2]) / 2):
                shore_link[shore][scaffold] = ["+",scaffolds_info[scaffold][0],0]
                shore_link[scaffold][shore] = [scaffolds_info[scaffold][0].translate(trans_strand),"-",0]
            else:
                shore_link[scaffold][shore] = [scaffolds_info[scaffold][0],"+",0]
                shore_link[shore][scaffold] = ["-",scaffolds_info[scaffold][0].translate(trans_strand),0]

temp = {}
for shore in shores_order:
    print(shore,sorted(shore_link[shore].items(),key = lambda item : item[1][2],reverse = True))
    for scaffold,value in shore_link[shore].items():
        temp[shore + scaffold] = value[2]
print(temp)
mean = sum(temp.values()) / len(temp.values())
std = np.std(list(temp.values()))
for shore in shores_order:
    for scaffold2,value in shore_link[shore].items():
        value[2] = (value[2] - mean) / std
    print(shore,sorted(shore_link[shore].items(),key = lambda item : item[1][2],reverse = True))

scaffold_link = defaultdict(dict)
with open(hifi_link_file) as file:
    for line in file:
        line = line.strip().split()
        scaffold1 = line[0][ : 10]
        strand1 = line[0][10]
        scaffold2 = line[0][11 : 21]
        strand2 = line[0][21]
        # keep links that
        # scaffolds fit the order and orientation in global combination
        if((scaffold1 in scaffolds_info and scaffold2 in scaffolds_info) 
        and ((strand1 == scaffolds_info[scaffold1][0] and strand2 == scaffolds_info[scaffold2][0] and (scaffolds_info[scaffold1][1] + scaffolds_info[scaffold1][2]) / 2 <= (scaffolds_info[scaffold2][1] + scaffolds_info[scaffold2][2]) / 2) 
            or (strand1.translate(trans_strand) == scaffolds_info[scaffold1][0] and strand2.translate(trans_strand) == scaffolds_info[scaffold2][0] and (scaffolds_info[scaffold1][1] + scaffolds_info[scaffold1][2]) / 2 >= (scaffolds_info[scaffold2][1] + scaffolds_info[scaffold2][2]) / 2))):
            if(scaffold1 not in scaffold_link):
                scaffold_link[scaffold1] = {}
            if(scaffold2 not in scaffold_link):
                scaffold_link[scaffold2] = {}
            scaffold_link[scaffold1][scaffold2] = [strand1,strand2,int(line[1]) / (length[scaffold1] * length[scaffold2])]
            scaffold_link[scaffold2][scaffold1] = [strand2.translate(trans_strand),strand1.translate(trans_strand),int(line[1]) / (length[scaffold1] * length[scaffold2])]

for scaffold1,scaffold2 in itertools.combinations([scaffold for scaffold in available_scaffolds if(scaffold in scaffolds_info)],2):
    if(scaffold2 not in scaffold_link[scaffold1]):
        if((scaffolds_info[scaffold1][1] + scaffolds_info[scaffold1][2]) / 2 >= (scaffolds_info[scaffold2][1] + scaffolds_info[scaffold2][2]) / 2):
            scaffold_link[scaffold2][scaffold1] = [scaffolds_info[scaffold2][0],scaffolds_info[scaffold1][0],0]
            scaffold_link[scaffold1][scaffold2] = [scaffolds_info[scaffold1][0].translate(trans_strand),scaffolds_info[scaffold2][0].translate(trans_strand),0]
        else:
            scaffold_link[scaffold1][scaffold2] = [scaffolds_info[scaffold1][0],scaffolds_info[scaffold2][0],0]
            scaffold_link[scaffold2][scaffold1] = [scaffolds_info[scaffold2][0].translate(trans_strand),scaffolds_info[scaffold1][0].translate(trans_strand),0]

temp = {}
for scaffold1,values in scaffold_link.items():
    print(scaffold1,sorted(values.items(),key = lambda item : item[1][2],reverse = True))
    for scaffold2,value in values.items():
        temp_list = [scaffold1,scaffold2]
        temp["".join((sorted(temp_list)))] = value[2]
if(len(temp) > 1):
    print(temp)
    mean = sum(temp.values()) / len(temp.values())
    std = np.std(list(temp.values()))
    min_value = sys.maxsize
    max_value = -sys.maxsize - 1
    for scaffold1,values in scaffold_link.items():
        for scaffold2,value in values.items():
            value[2] = (value[2] - mean) / std
            if(value[2] > max_value):
                max_value = value[2]
            if(value[2] < min_value):
                min_value = value[2]
        print(scaffold1,sorted(values.items(),key = lambda item : item[1][2],reverse = True))
    for scaffold1,values in scaffold_link.items():
        for scaffold2,value in values.items():
            value[2] = (value[2] - min_value) / (max_value - min_value)
        print(scaffold1,sorted(values.items(),key = lambda item : item[1][2],reverse = True))
    for scaffold1,values in scaffold_link.items():
        for scaffold2,value in values.items():
            value[2] = scaffold_scale * (homologous_coefficient - homologous[" ".join((sorted([scaffold1,scaffold2])))]) * value[2] # 1
        print(scaffold1,sorted(values.items(),key = lambda item : item[1][2],reverse = True))

index = {}
index_reverse = {}
value_ranges = []
for i,scaffold in enumerate(available_scaffolds):
    index[scaffold] = i
    index_reverse[i] = scaffold
    # 1,2 stands for each haplotype, 3 for both, 0 for none, nobody wants you :(
    """if(GFA.nodes[scaffold]["length"] < 100000):
        if(scaffold in connected_components_index and scaffold not in sources and connected_components_len[connected_components_index[scaffold]] >= centromere_len and len(connected_components[connected_components_index[scaffold]]) > 4):
            value_ranges.append([-1,0,1,2])
        else:
            value_ranges.append([-1,1,2])
    else:
        value_ranges.append([-1,1])"""
    if(scaffold in connected_components_index and connected_components_len[connected_components_index[scaffold]] >= centromere_len and len(connected_components[connected_components_index[scaffold]]) > 4):
        if(scaffold in sources or len(GFA[scaffold].keys()) >= 4):
            value_ranges.append([0,1,2,3])
        else:
            value_ranges.append([1,2])
    else:
        value_ranges.append([0,1,2])

whole_set = defaultdict(list)
global_phase = []   
for value_range in value_ranges:
    global_phase.append(random.choice(value_range))

"""if(chr == "chr1"):
    global_phase[index["sfd000001l"]] = 1
    global_phase[index["sfd000005l"]] = 1
    global_phase[index["sfd000006l"]] = 1
    global_phase[index["utg000168l"]] = 0
    global_phase[index["sfd000003l"]] = -1
    global_phase[index["sfd000004l"]] = -1
    global_phase[index["sfd000002l"]] = -1
    global_phase[index["utg000078l"]] = -1
    global_phase[index["utg000355l"]] = -1
elif(chr == "chr4"):
    global_phase[index["utg000009l"]] = 1
    global_phase[index["utg000064l"]] = 1
    global_phase[index["sfd000002l"]] = 1
    global_phase[index["sfd000001l"]] = 1
    global_phase[index["sfd000012l"]] = 1
    global_phase[index["utg000023l"]] = -1
    global_phase[index["utg000005l"]] = -1
    global_phase[index["sfd000008l"]] = -1
    global_phase[index["sfd000003l"]] = -1
    global_phase[index["utg000099l"]] = -1
    global_phase[index["utg000023l"]] = -1
    global_phase[index["sfd000005l"]] = -1
    global_phase[index["sfd000011l"]] = 2
    global_phase[index["sfd000007l"]] = -1
elif(chr == "chr10"):
    global_phase[index["utg000003l"]] = -1
    global_phase[index["sfd000001l"]] = -1
    global_phase[index["sfd000002l"]] = 1
    global_phase[index["sfd000004l"]] = 1
    global_phase[index["utg000032l"]] = -1
    global_phase[index["utg000037l"]] = 2
    global_phase[index["sfd000005l"]] = 1
    global_phase[index["sfd000003l"]] = 2
elif(chr == "chr11"):
    global_phase[index["sfd000001l"]] = 1
    global_phase[index["utg000004l"]] = 1
    global_phase[index["utg000008l"]] = 1
    global_phase[index["sfd000002l"]] = -1
    global_phase[index["sfd000003l"]] = -1
    global_phase[index["sfd000005l"]] = -1
    global_phase[index["sfd000007l"]] = 2
    global_phase[index["sfd000008l"]] = 2
    global_phase[index["sfd000006l"]] = -1
    global_phase[index["sfd000009l"]] = -1
elif(chr == "chr12"):
    global_phase[index["utg000009l"]] = -1
    global_phase[index["sfd000007l"]] = -1
    global_phase[index["sfd000011l"]] = -1
    global_phase[index["utg000031l"]] = 0
    global_phase[index["utg000001l"]] = -1
    global_phase[index["sfd000003l"]] = 1
    global_phase[index["utg000040l"]] = 1
    global_phase[index["utg000047l"]] = 2
    global_phase[index["sfd000010l"]] = 1
    global_phase[index["utg000027l"]] = 1
    global_phase[index["utg000046l"]] = 2
    global_phase[index["sfd000006l"]] = 2
    global_phase[index["utg000012l"]] = 2
    global_phase[index["utg000035l"]] = 2
    global_phase[index["utg000050l"]] = 2
    global_phase[index["utg000045l"]] = 2
elif(chr == "chr20"):
    global_phase[index["sfd000001l"]] = 1
    global_phase[index["sfd000012l"]] = 1
    global_phase[index["sfd000013l"]] = -1
    global_phase[index["sfd000014l"]] = -1
    global_phase[index["sfd000015l"]] = -1
    global_phase[index["utg000036l"]] = -1
    global_phase[index["sfd000007l"]] = -1
    global_phase[index["sfd000016l"]] = -1
    global_phase[index["utg000012l"]] = -1
    global_phase[index["utg000186l"]] = -1
    global_phase[index["utg000131l"]] = -1
    global_phase[index["utg000126l"]] = -1
    global_phase[index["utg000115l"]] = 2
    global_phase[index["utg000145l"]] = -1
    global_phase[index["utg000153l"]] = 2
    global_phase[index["utg000100l"]] = 2
    global_phase[index["utg000146l"]] = -1
    global_phase[index["sfd000009l"]] = -1
    global_phase[index["sfd000010l"]] = 2"""

global_objfunction = calculate_score(global_phase)
local_objfunction = global_objfunction
local_phase = copy.copy(global_phase)
temp_objfunction = global_objfunction
temp_phase = copy.copy(local_phase)
for i in range(1000):
    print(i,global_objfunction)
    j = 0
    while(j < 100):
        flip_position = random.sample(available_scaffolds, 1)[0]
        flip_position_index = index[flip_position]
        temp_phase[flip_position_index] = random.choice([value for value in value_ranges[flip_position_index] if value != temp_phase[flip_position_index]])
        temp_objfunction = calculate_score(temp_phase)
        if(temp_objfunction[5] > local_objfunction[5]):
            local_objfunction = temp_objfunction
            local_phase = copy.copy(temp_phase)
            j = 0
        else:
            temp_objfunction = local_objfunction
            temp_phase = copy.copy(local_phase)
            j += 1

    whole_set[" ".join(local_objfunction[0]) + "|" + " ".join(local_objfunction[1])] = local_objfunction[2:]  # type: ignore

    if(local_objfunction[5] > global_objfunction[5]):
        global_objfunction = local_objfunction
        global_phase = copy.copy(local_phase)
    else:
        local_objfunction = global_objfunction
        local_phase = copy.copy(global_phase)
    
    disturb = random.sample(available_scaffolds, len(available_scaffolds) // 2)
    for unitig in disturb:
        flip_position_index = index[unitig]
        local_phase[flip_position_index] = random.choice([value for value in value_ranges[flip_position_index] if value != local_phase[flip_position_index]])
    local_objfunction = calculate_score(local_phase)

mat_order,pat_order,_,_,_,_ = calculate_score(global_phase)

for key,value in sorted(whole_set.items(), key = lambda item : item[1][3], reverse = True):
    print(key,value)

mat_available = np.ones(len(mat_order))
pat_available = np.ones(len(pat_order))
for i,scaffold in enumerate(mat_order):
    if(mat_available[i]):
        # former two are homologous
        if(i > 1 and (mat_available[i - 1] and (mat_order[i - 1],scaffold) in GFA.edges()) and (mat_available[i - 2] and (mat_order[i - 2],scaffold) in GFA.edges()) and homologous[" ".join((sorted([mat_order[i - 1],mat_order[i - 2]])))] > homologous_coefficient):
            # keep latter one
            print(mat_order[i],mat_order[i - 1])
            mat_available[i - 2] = 0
        # latter two are homologous
        if(i < len(mat_order) - 2 and (mat_available[i + 1] and (scaffold,mat_order[i + 1]) in GFA.edges()) and (mat_available[i + 2] and (scaffold,mat_order[i + 2]) in GFA.edges()) and homologous[" ".join((sorted([mat_order[i + 1],mat_order[i + 2]])))] > homologous_coefficient):
            # keep latter one
            print(mat_order[i],mat_order[i + 2])
            mat_available[i + 1] = 0
        # got edge with node witnin two or three radius, and do not have edge with next node, then discard nodes in between
        if(i < len(mat_order) - 2 and (mat_available[i + 2] and (scaffold,mat_order[i + 2]) in GFA.edges()) and (mat_available[i + 1] and (scaffold,mat_order[i + 1]) not in GFA.edges() and (mat_order[i + 1],mat_order[i + 2]) not in GFA.edges())):
            print(mat_order[i],mat_order[i + 1])
            mat_available[i + 1] = 0
        if(i < len(mat_order) - 3 and (mat_available[i + 3] and (scaffold,mat_order[i + 3]) in GFA.edges()) and (mat_available[i + 1] and mat_available[i + 2] and (scaffold,mat_order[i + 1]) not in GFA.edges() and (mat_order[i + 2],mat_order[i + 3]) not in GFA.edges())):
            print(mat_order[i],mat_order[i + 1])
            print(mat_order[i],mat_order[i + 2])
            mat_available[i + 1] = 0
            mat_available[i + 2] = 0
        # deal with small isolate scaffolds
        if(GFA.nodes[scaffold]["length"] < 100000 and mat_available[i]):
            left = np.where(mat_available[ : i] == 1)[0]
            right = np.where(mat_available[i + 1 :] == 1)[0] + i + 1
            left_most = max(left) if(len(left) > 0) else -1
            right_most = min(right) if(len(right) > 0) else -1
            print(scaffold,left,left_most,right,right_most)
            if(
                (left_most != -1 and right_most != -1 and (mat_order[left_most],scaffold) not in GFA.edges() and (scaffold,mat_order[right_most]) not in GFA.edges())
                or
                (left_most != -1 and right_most == -1 and (mat_order[left_most],scaffold) not in GFA.edges())
                or
                (left_most == -1 and right_most != -1 and (scaffold,mat_order[right_most]) not in GFA.edges())
            ):
                print(mat_order[i])
                mat_available[i] = 0
for i,scaffold in enumerate(pat_order):
    if(pat_available[i]):
        # former two are homologous
        if(i > 1 and (pat_available[i - 1] and (pat_order[i - 1],scaffold) in GFA.edges()) and (pat_available[i - 2] and (pat_order[i - 2],scaffold) in GFA.edges()) and homologous[" ".join((sorted([pat_order[i - 1],pat_order[i - 2]])))] > homologous_coefficient):
            # keep latter one
            print(pat_order[i],pat_order[i - 1])
            pat_available[i - 2] = 0
        # latter two are homologous
        if(i < len(pat_order) - 2 and (pat_available[i + 1] and (scaffold,pat_order[i + 1]) in GFA.edges()) and (pat_available[i + 2] and (scaffold,pat_order[i + 2]) in GFA.edges()) and homologous[" ".join((sorted([pat_order[i + 1],pat_order[i + 2]])))] > homologous_coefficient):
            # keep latter one
            print(pat_order[i],pat_order[i + 2])
            pat_available[i + 1] = 0
        # got edge with node witnin two or three radius, and do not have edge with next node, then discard nodes in between
        if(i < len(pat_order) - 2 and (pat_available[i + 2] and (scaffold,pat_order[i + 2]) in GFA.edges()) and (pat_available[i + 1] and (scaffold,pat_order[i + 1]) not in GFA.edges() and (pat_order[i + 1],pat_order[i + 2]) not in GFA.edges())):
            print(pat_order[i],pat_order[i + 1])
            pat_available[i + 1] = 0
        if(i < len(pat_order) - 3 and (pat_available[i + 3] and (scaffold,pat_order[i + 3]) in GFA.edges()) and (pat_available[i + 1] and pat_available[i + 2] and (scaffold,pat_order[i + 1]) not in GFA.edges() and (pat_order[i + 2],pat_order[i + 3]) not in GFA.edges())):
            print(pat_order[i],pat_order[i + 1])
            print(pat_order[i],pat_order[i + 2])
            pat_available[i + 1] = 0
            pat_available[i + 2] = 0
        # deal with small isolate scaffolds
        if(GFA.nodes[scaffold]["length"] < 100000 and pat_available[i]):
            left = np.where(pat_available[ : i] == 1)[0]
            right = np.where(pat_available[i + 1 :] == 1)[0] + i + 1
            left_most = max(left) if(len(left) > 0) else -1
            right_most = min(right) if(len(right) > 0) else -1
            print(scaffold,left,left_most,right,right_most)
            if(
                (left_most != -1 and right_most != -1 and (pat_order[left_most],scaffold) not in GFA.edges() and (scaffold,pat_order[right_most]) not in GFA.edges())
                or
                (left_most != -1 and right_most == -1 and (pat_order[left_most],scaffold) not in GFA.edges())
                or
                (left_most == -1 and right_most != -1 and (scaffold,pat_order[right_most]) not in GFA.edges())
            ):
                print(pat_order[i])
                pat_available[i] = 0

mat_order = [mat_order[i] for i,available in enumerate(mat_available) if(available)]
pat_order = [pat_order[i] for i,available in enumerate(pat_available) if(available)]
print(mat_order,pat_order)
writer = open("to_be_phased_centromere.fa","w")
mat_sequence = ""
pat_sequence = ""
match = 0
for i,scaffold in enumerate(mat_order):
    if(scaffolds_info[scaffold][0] == "+"):
        mat_sequence += GFA.nodes[scaffold]["seq"][match : ]
    else:
        mat_sequence += GFA.nodes[scaffold]["seq"].translate(transtable)[::-1][match : ]
    if(i != len(mat_order) - 1):
        if((scaffold,mat_order[i + 1]) in GFA.edges()):
            match = GFA[scaffold][mat_order[i + 1]]["match"]
        else:
            match = 0
            mat_sequence += "N" * 500
match = 0
for i,scaffold in enumerate(pat_order):
    if(scaffolds_info[scaffold][0] == "+"):
        pat_sequence += GFA.nodes[scaffold]["seq"][match : ]
    else:
        pat_sequence += GFA.nodes[scaffold]["seq"].translate(transtable)[::-1][match : ]
    if(i != len(pat_order) - 1):
        if((scaffold,pat_order[i + 1]) in GFA.edges()):
            match = GFA[scaffold][pat_order[i + 1]]["match"]
        else:
            match = 0
            pat_sequence += "N" * 500
writer.write(">cen000001l\n")
writer.write(mat_sequence + "\n")
writer.write(">cen000002l\n")
writer.write(pat_sequence + "\n")
writer.close()