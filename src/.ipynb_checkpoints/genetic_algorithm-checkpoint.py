#!/usr/bin/env python
from collections import defaultdict
from scipy import signal
import sys
import copy
import math
import random
import itertools
import numpy as np
import networkx as nx

chr = sys.argv[1]
gap_start = int(sys.argv[2])
gap_end = int(sys.argv[3])
paf_file = sys.argv[4]
gfa_file = sys.argv[5]
fa_file = sys.argv[6]
ploid = sys.argv[7]
output_fa_file = sys.argv[8]
output_combination_file = sys.argv[9]
transtable = str.maketrans('ATGC','TACG')
trans_strand = str.maketrans('+-',"-+")
adaptivity_shift = -1
index_shift = 10000
interval_middle_shift = 10000
global_length_threshold = 0.8

def is_concord(individual):
    print(individual)
    length = sum([GFA.nodes[key]["length"] for key in individual.keys()])
    edge_concord = {key : 0 for key in individual.keys()}
    for key1,key2 in itertools.combinations(individual.keys(),2): # look at neighbors
        if(
            (key1,key2) in GFA.edges()
                and
            (
                ((individual[key1][1] + individual[key1][2]) / 2 <= (individual[key2][1] + individual[key2][2]) / 2 and GFA[key1][key2]["strand1"] == individual[key1][0] and GFA[key1][key2]["strand2"] == individual[key2][0])
                    or
                ((individual[key2][1] + individual[key2][2]) / 2 <= (individual[key1][1] + individual[key1][2]) / 2 and GFA[key1][key2]["strand1"].translate(trans_strand) == individual[key1][0] and GFA[key1][key2]["strand2"].translate(trans_strand) == individual[key2][0])
            )
        ):
            edge_concord[key1] += 1
            edge_concord[key2] += 1
    print(edge_concord)
    component_concord = {key1 : {key2 : 0 for key2 in individual.keys() if(key1 != key2)} for key1 in individual.keys()}
    component_disconcord = {key1 : {key2 : 0 for key2 in individual.keys() if(key1 != key2)} for key1 in individual.keys()}
    for key1,key2 in itertools.combinations(individual.keys(),2): # look at components
        if(any(key1 in component and key2 in component for component in connected_components) and homologous["".join(sorted([key1,key2])) + (" same" if(individual[key1] == individual[key2]) else " diff")] < 0.8):
            if(key1 + key2 in strand_src_dst):
                strand1 = strand_src_dst[key1 + key2][0]
                strand2 = strand_src_dst[key1 + key2][1]
            elif(key2 + key1 in strand_src_dst):
                strand1 = strand_src_dst[key2 + key1][1].translate(trans_strand)
                strand2 = strand_src_dst[key2 + key1][0].translate(trans_strand)
            else:
                continue

            if(
                ((individual[key1][1] + individual[key1][2]) / 2 <= (individual[key2][1] + individual[key2][2]) / 2 and strand1 == individual[key1][0] and strand2 == individual[key2][0])
                    or
                ((individual[key2][1] + individual[key2][2]) / 2 <= (individual[key1][1] + individual[key1][2]) / 2 and strand1.translate(trans_strand) == individual[key1][0] and strand2.translate(trans_strand) == individual[key2][0])
            ):
                component_concord[key1][key2] += 1
                component_concord[key2][key1] += 1
            else:
                component_disconcord[key1][key2] += 1
                component_disconcord[key2][key1] += 1    
    # print(component_concord)
    # print(component_disconcord)
    new_individual = {}
    discard = {}
    for key in individual.keys():
        if(
            (
                edge_concord[key] != 0 
                    and 
                (GFA.nodes[key]["length"] > 100000 or (GFA.nodes[key]["length"] <= 100000 and edge_concord[key] == len(GFA[key].keys())))
            )
                or 
            (
                edge_concord[key] == 0 
                    and
                ((GFA.in_degree[key] == 0 and GFA.out_degree[key] == 0) or all(edge_concord[neighbor] == 0 for neighbor in GFA[key].keys() if(neighbor in edge_concord)) or GFA.nodes[key]["length"] > 100000)
            )
        ): # concord(small nodes require completly concord, big nodes do not), isolate, all neighbors are disconcord
            new_individual[key] = individual[key]
        else:
            discard[key] = individual[key]
    for key in sorted(discard.keys(),key = lambda key : (-edge_concord[key],-GFA.nodes[key]["length"])):
        if(edge_concord[key] != 0 and edge_concord[key] == len([neighbor for neighbor in GFA[key].keys() if(neighbor in new_individual)])):
            new_individual[key] = individual[key]
            del discard[key]
    
    # print(new_individual)

    while(any(component_disconcord[key1][key2] > 0 for key1,key2 in itertools.combinations(new_individual.keys(),2))): # any disconcord exist
        target = sorted(new_individual.keys(),key = lambda key1 : (-sum([value for key2,value in component_disconcord[key1].items() if(key2 in new_individual)]),GFA.nodes[key]["length"]))[0] # delete the most disconcord (1) while shortest(2) node to reach concord
        discard[target] = new_individual[target]
        del new_individual[target]
    print(new_individual)
    for target in sorted(discard.keys(), key = lambda key : GFA.nodes[key]["length"],reverse = True): # try to retrieve discarded node, and place it to appropriate position according to surrounding
        for component in connected_components:
            if(target in component):
                left_node = ""
                left_min_distance = sys.maxsize
                right_node = ""
                right_min_distance = sys.maxsize
                strand1 = ""
                strand2 = ""
                temp = -1
                for key in component & new_individual.keys():
                    if(target + key in strand_src_dst):
                        strand1 = strand_src_dst[target + key][0]
                        strand2 = strand_src_dst[target + key][1]
                    elif(key + target in strand_src_dst):
                        strand1 = strand_src_dst[key + target][1].translate(trans_strand)
                        strand2 = strand_src_dst[key + target][0].translate(trans_strand)
                    else:
                        continue
                    if(strand1 == individual[target][0] and strand2 == individual[key][0]):
                        # key on right side of target
                        temp = nx.shortest_path_length(GFA,target,key)
                        if(temp < right_min_distance or (temp == right_min_distance and new_individual[key][1] + new_individual[key][2] < new_individual[right_node][1] + new_individual[right_node][2])):
                            right_min_distance = temp
                            right_node = key
                    elif(strand1.translate(trans_strand) == individual[target][0] and strand2.translate(trans_strand) == individual[key][0]):
                        # key on left side of target
                        temp = nx.shortest_path_length(GFA,target,key)
                        if(temp < left_min_distance or (temp == left_min_distance and new_individual[key][1] + new_individual[key][2] > new_individual[left_node][1] + new_individual[left_node][2])):
                            left_min_distance = temp
                            left_node = key
                    else:
                        continue
                if(len(gene_banks[target]) == 1 and GFA.nodes[target]["length"] > 100000):
                    if(left_node != "" and right_node != ""):
                        temp = (new_individual[left_node][1] + new_individual[left_node][2] + new_individual[right_node][1] + new_individual[right_node][2]) / 4            
                    elif(left_node != ""):
                        temp = (new_individual[left_node][1] + new_individual[left_node][2]) / 2 + interval_middle_shift
                    elif(right_node != ""):
                        temp = (new_individual[right_node][1] + new_individual[right_node][2]) / 2 - interval_middle_shift
                    else:
                        break
                    # print(target,left_node,right_node)
                    new_individual[target] = discard[target]
                    new_individual[target][1] = int(temp - GFA.nodes[target]["length"] / 2)
                    new_individual[target][2] = int(temp + GFA.nodes[target]["length"] / 2)
                    del discard[target]      
                    break
                else:
                    for gene in gene_banks[target]:
                        if(left_node != "" and right_node != ""):
                            strand = strand_src_dst[left_node + target][1] if(left_node + target in strand_src_dst) else strand_src_dst[target + left_node][0].translate(trans_strand)
                            if(strand == gene[0] and new_individual[left_node][1] + new_individual[left_node][2] < gene[1] + gene[2] < new_individual[right_node][1] + new_individual[right_node][2]):
                                # print(target,left_node,right_node)
                                new_individual[target] = gene
                                del discard[target]
                                break
                        elif(left_node != ""):
                            strand = strand_src_dst[left_node + target][1] if(left_node + target in strand_src_dst) else strand_src_dst[target + left_node][0].translate(trans_strand)
                            if(strand == gene[0] and new_individual[left_node][1] + new_individual[left_node][2] < gene[1] + gene[2]):
                                # print(target,left_node,right_node)
                                new_individual[target] = gene
                                del discard[target]
                                break
                        elif(right_node != ""):
                            strand = strand_src_dst[target + right_node][0] if(target + right_node in strand_src_dst) else strand_src_dst[right_node + target][1].translate(trans_strand)
                            if(strand == gene[0] and gene[1] + gene[2] < new_individual[right_node][1] + new_individual[right_node][2]):
                                # print(target,left_node,right_node)
                                new_individual[target] = gene
                                del discard[target]
                                break
                break

    component_concord = {key1 : {key2 : 0 for key2 in new_individual.keys() if(key1 != key2)} for key1 in new_individual.keys()}
    for key1,key2 in itertools.combinations(new_individual.keys(),2): # look at components
        if(any(key1 in component and key2 in component for component in connected_components) and homologous["".join(sorted([key1,key2])) + (" same" if(new_individual[key1] == new_individual[key2]) else " diff")] < 0.8):
            if(key1 + key2 in strand_src_dst):
                strand1 = strand_src_dst[key1 + key2][0]
                strand2 = strand_src_dst[key1 + key2][1]
            elif(key2 + key1 in strand_src_dst):
                strand1 = strand_src_dst[key2 + key1][1].translate(trans_strand)
                strand2 = strand_src_dst[key2 + key1][0].translate(trans_strand)
            else:
                continue

            if(
                ((new_individual[key1][1] + new_individual[key1][2]) / 2 <= (new_individual[key2][1] + new_individual[key2][2]) / 2 and strand1 == new_individual[key1][0] and strand2 == new_individual[key2][0])
                    or
                ((new_individual[key2][1] + new_individual[key2][2]) / 2 <= (new_individual[key1][1] + new_individual[key1][2]) / 2 and strand1.translate(trans_strand) == new_individual[key1][0] and strand2.translate(trans_strand) == new_individual[key2][0])
            ):
                component_concord[key1][key2] += 1
                component_concord[key2][key1] += 1

    if(sum([GFA.nodes[key]["length"] for key in new_individual.keys()]) < length * global_length_threshold): # length control
        return 0,0,{}
    else:
        return 1,sum([sum(component_concord[key1][key2] for key2 in component_concord[key1] if(key1 != key2 and key2 in new_individual)) for key1 in component_concord if(key1 in new_individual)]),new_individual
    
def assess_individual(id,individual):
    depth = sum([value[3] for value in individual.values()])
    # temp_depth = 0
    # for key,value in individual.items():
    #     temp_depth += sum(depth[key][value[0]][1][value[1] - depth[key][value[0]][0] : value[2] - depth[key][value[0]][0] + 1])
    # temp_depth /= len(individual)

    # distance = np.var([(value[1] + value[2]) / 2 for value in individual.values()])
    distance = 0
    for key1,key2 in itertools.combinations(individual.keys(),2):
        distance += abs((individual[key1][1] + individual[key1][2]) / 2 - (individual[key2][1] + individual[key2][2]) / 2)
    # distance /= len(individual) * (len(individual) - 1)

    coverage = 0
    min_index = min(min(value[1] for value in individual.values()),gap_start) # gap_start
    max_index = max(max(value[2] for value in individual.values()),gap_end) # gap_end
    array = np.zeros(max_index - min_index + 1)
    for component in connected_components_plus:
        if(any(node in individual for node in component)):
            left_most = max(min_index,min([individual[node][1] for node in component if(node in individual)])) - min_index
            right_most = min(max([individual[node][2] for node in component if(node in individual)]),max_index) - min_index
            temp_array = np.zeros(max_index - min_index + 1)
            for node in component:
                if(node in individual):
                    temp_array[max(0,individual[node][1] - min_index) : min(individual[node][2] - min_index + 1,max_index + 1)] += 1
                    array[max(0,individual[node][1] - min_index) : min(individual[node][2] - min_index + 1,max_index + 1)] += 1
            mask = (temp_array[left_most : right_most + 1] == 0)
            array[left_most : right_most + 1][mask] -= 1
    # for value in individual.values():
    #     array[max(0,value[1] - min_index) : min(value[2] - min_index + 1,max_index + 1)] += 1
    print(id, len(np.where(array >= 2)[0]) , len(np.where(array == 1)[0]) , len(np.where(array == 0)[0]) , len(np.where(array > 2)[0]) , len(np.where(array < 0)[0]) , len(np.where(array >= 2)[0]) * 2 - len(np.where(array == 1)[0]) - len(np.where(array == 0)[0]) * 2 - len(np.where(array < 0)[0]) * 2) #   - len(np.where(array > 2)[0]) * 0.5
    if(ploid == "haploid"):
        coverage = len(np.where(array == 1)[0]) * 2 + len(np.where(array == 2)[0]) - len(np.where(array == 0)[0]) * 2 - len(np.where(array < 0)[0]) * 2.5 - len(np.where(array > 2)[0]) * 2  #   - len(np.where(array > 2)[0]) * 0.5
    elif(ploid == "diploid"):
        coverage = len(np.where(array >= 2)[0]) * 2 - len(np.where(array == 1)[0]) - len(np.where(array == 0)[0]) * 2 - len(np.where(array < 0)[0]) * 2.5  #   - len(np.where(array > 2)[0]) * 0.5

    return depth,distance,coverage

GFA = nx.DiGraph()
with open(gfa_file) as file:
    for line in file:
        line = line.strip().split()
        if(line[0] == "S"):
            GFA.add_node(line[1],length = len(line[2]),seq = line[2])
        elif(line[0] == "L" and line[1] != line[3]):
            GFA.add_edge(line[1],line[3],strand1 = line[2],strand2 = line[4],match = int(line[5][:-1]))

PAF_info = {}
with open(paf_file) as file:
    for line in file:
        line = line.strip().split()
        if(line[0] not in PAF_info):
            PAF_info[line[0]] = defaultdict(list)
        PAF_info[line[0]][line[5]].append(line) 

PAF = {}
for scaffold,values in PAF_info.items():
    if(chr in values.keys()):
        array = np.zeros(GFA.nodes[scaffold]["length"])
        for value in values[chr]:
            array[int(value[2]) : int(value[3]) + 1] += 1
        temp = [value for value in values[chr] if(max(gap_start,int(value[7])) <= min(gap_end,int(value[8])))]
        if((len(temp) > 2 * (sum(map(len,values.values())) - len(temp)) or sum([int(value[3]) - int(value[2]) for value in temp]) > 2 * sum([sum([int(value[3]) - int(value[2]) for value in values[key] if(key != chr or max(gap_start,int(value[7])) > min(gap_end,int(value[8])))]) for key in values.keys()])) and len(np.where(array > 0)[0]) > 0.1 * GFA.nodes[scaffold]["length"]):
            # print(len(temp),(sum(map(len,values.values())) - len(temp)),sum([int(value[3]) - int(value[2]) for value in temp]),sum([sum([int(value[3]) - int(value[2]) for value in values[key] if(key != chr or max(gap_start,int(value[7])) > min(gap_end,int(value[8])))]) for key in values.keys()]))
            PAF[scaffold] = defaultdict(list)
            for info in values[chr]:
                if(max(gap_start,int(info[7])) <= min(gap_end,int(info[8]))):
                    PAF[scaffold][info[4]].append([int(info[1]),int(info[2]),int(info[3]),int(info[7]),int(info[8])]) # query_len,query_start,query_end,ref_start,ref_end
if(sum(GFA.nodes[scaffold]["length"] for scaffold in PAF.keys()) < 0.9 * (gap_end - gap_start)):
    for scaffold,values in PAF_info.items():
        if(chr in values.keys()):
            array = np.zeros(GFA.nodes[scaffold]["length"])
            for value in values[chr]:
                if(max(gap_start,int(value[7])) <= min(gap_end,int(value[8]))):
                    array[int(value[2]) : int(value[3]) + 1] += 1
            if(len(np.where(array > 0)[0]) > 0.4 * GFA.nodes[scaffold]["length"]):
                PAF[scaffold] = defaultdict(list)
                for info in values[chr]:
                    if(max(gap_start,int(info[7])) <= min(gap_end,int(info[8]))):
                        PAF[scaffold][info[4]].append([int(info[1]),int(info[2]),int(info[3]),int(info[7]),int(info[8])]) # query_len,query_start,query_end,ref_start,ref_end
if(len(PAF) == 0 or sum(GFA.nodes[scaffold]["length"] for scaffold in PAF.keys()) < 0.2 * (gap_end - gap_start)):
    raise Exception("Heterozygosity between reference and query is too high, can't find much similarity!\n")

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
    connected_components.append(temp)

connected_components_plus = copy.copy(connected_components)
for node in GFA.nodes():
    if(GFA.in_degree[node] == 0 and GFA.out_degree[node] == 0):
        connected_components_plus.append([node])

strand_src_dst = {}
for component in connected_components:
    for key1,key2 in itertools.combinations(component,2):
        found_flag = 0
        for strand_src in ["+","-"]:
            if(found_flag):
                break
            queue = [strand_src + key1]
            visited = set()
            while(len(queue) > 0):
                strand = queue[0][0]
                node = queue[0][1:]
                if(node == key2):
                    found_flag = 1
                    strand_src_dst[key1 + key2] = [strand_src,strand]
                    break
                queue.pop(0)
                while("+" + node in queue):
                    queue.remove("+" + node)
                while("-" + node in queue):
                    queue.remove("-" + node)  
                visited.add(node)
                for neighbor in GFA[node].keys():
                    if(strand == GFA[node][neighbor]["strand1"] and neighbor not in visited):
                        queue.append(GFA[node][neighbor]["strand2"] + neighbor)

gene_banks = defaultdict(list)
# get intervals basing on dynamic programming
for scaffold,values in PAF.items():
    # print(scaffold)
    for strand,infos in values.items():
        qry_id = {}
        ref_id = {}
        intervals = []
        min_index = -1
        max_index = -1
        # print(strand)
        length = infos[0][0]
        if(strand == "+"):
            info_sorted = sorted(infos,key = lambda info : int(info[1]))
            min_index = min([info[3] - info[1] for info in info_sorted]) - index_shift
            max_index = max([info[4] + info[0] - info[2] - 1 for info in info_sorted]) + index_shift
            id = 0
            count = 0
            for info in info_sorted:
                lower = info[3] - info[1]
                upper = max(info[4] + info[0] - info[2] - 1,lower + length - 1)
                if(str(info[1]) + " " + str(info[2]) in qry_id):
                    id = qry_id[str(info[1]) + " " + str(info[2])]
                else:
                    id = count
                    count += 1
                qry_id[str(info[1]) + " " + str(info[2])] = id
                ref_id[str(info[3]) + " " + str(info[4])] = [id,lower,upper]
                # print(info[0],info[1],info[2],info[3],info[4])
        else:
            info_sorted = sorted(infos,key = lambda info : int(info[2]),reverse = True)
            min_index = min([info[3] - (info[0] - info[2] - 1) for info in info_sorted]) - index_shift
            max_index = max([info[4] + info[1] for info in info_sorted]) + index_shift
            id = 0
            count = 0
            for info in info_sorted:
                lower = info[3] - (info[0] - info[2] - 1)
                upper = max(info[4] + info[1],lower + length - 1)
                if(str(info[1]) + " " + str(info[2]) in qry_id):
                    id = qry_id[str(info[1]) + " " + str(info[2])]
                else:
                    id = count
                    count += 1
                qry_id[str(info[1]) + " " + str(info[2])] = id
                ref_id[str(info[3]) + " " + str(info[4])] = [id,lower,upper]
                # print(info[0],info[1],info[2],info[3],info[4])
        zipped = list(zip(*sorted(ref_id.items(),key = lambda item : int(item[0].split()[0]))))
        intervals = zipped[0]
        ids = list(zip(*zipped[1]))[0]
        """for key,value in sorted(ref_id.items(),key = lambda item : int(item[0].split()[0])):
            print(key,value)"""
        results = {}
        for index_start in range(len(ids)):
            ref_start = ref_id[intervals[index_start]][1]
            ref_end = ref_id[intervals[index_start]][2]
            index_end = index_start
            for index,interval in enumerate(intervals[index_start + 1 :]):
                interval = interval.split()
                interval_start = int(interval[0])
                interval_end = int(interval[1])
                if(interval_start <= ref_end):
                    index_end = index + index_start + 1
                elif(ref_end < interval_start):
                    break
            local_ids = ids[index_start : index_end + 1]
            dp = [1] * len(local_ids)
            for i in range(1,len(local_ids)):
                for j in range(i):
                    if(local_ids[j] <= local_ids[i]):
                        dp[i] = max(dp[i],dp[j] + 1)
            results[str(ref_start) + " " + str(ref_end)] = max(dp)
        # print(sorted(results.items(),key = lambda item : -item[1]))
        # according to dp decending order, choose non-overlaping interval in top 7 dp value as candidate intervals
        top = 0
        intervals = []
        last = -1
        count = 0
        array = np.zeros(max_index - min_index + 1,dtype = int)
        for key,value in sorted(results.items(),key = lambda item : item[1],reverse = True):
            # print(key,value)
            if(value != last):
                if(count != 0):
                    max_peaks = signal.find_peaks(array,distance = length)[0]
                    for index in max_peaks:
                        if(all(index + min_index < interval[0] or interval[1] < index + min_index for interval in intervals)):
                            start = max(0,index - length)
                            coverage = sum(array[start : start + length])
                            max_start = start
                            max_coverage = coverage
                            for i in range(start + 1,min(index + 1,len(array) - length + 1)):
                                coverage = coverage + array[i + length - 1] - array[i - 1]
                                if(coverage > max_coverage):
                                    max_coverage = coverage
                                    max_start = i
                            if(all(max_start + min_index + length - 1 < interval[0] or max_start + min_index > interval[1] for interval in intervals)): # non overlap
                                intervals.append([max_start + min_index,max_start + min_index + length - 1,max_coverage])
                    array = np.zeros(max_index - min_index + 1,dtype = int)
                    count = 0
                    top += 1
            if(top == 7):
                break
            key = key.split()
            lower = int(key[0]) - min_index
            upper = int(key[1]) - min_index
            array[lower : upper + 1] += 1
            last = value
            count += 1
        if(top < 7):
            max_peaks = signal.find_peaks(array,distance = length)[0]
            for index in max_peaks:
                if(all(index + min_index < interval[0] or interval[1] < index + min_index for interval in intervals)):
                    start = max(0,index - length)
                    coverage = sum(array[start : start + length])
                    max_start = start
                    max_coverage = coverage
                    for i in range(start + 1,min(index + 1,len(array) - length + 1)):
                        coverage = coverage + array[i + length - 1] - array[i - 1]
                        if(coverage > max_coverage):
                            max_coverage = coverage
                            max_start = i
                    if(all(max_start + min_index + length - 1 < interval[0] or max_start + min_index > interval[1] for interval in intervals)): # non overlap
                        intervals.append([max_start + min_index,max_start + min_index + length - 1,max_coverage])
        for interval in intervals:
            gene_banks[scaffold].append([strand,interval[0],interval[1],interval[2]])
            
gene_banks_combinations_count = 1
for key,value in gene_banks.items():
    print(key,value)
    gene_banks_combinations_count *= len(value)
if(len(gene_banks.keys()) > 200 or math.log2(gene_banks_combinations_count) > len(gene_banks.keys()) / 2):
    global_length_threshold = 0.5

kmers_pos = defaultdict(set)
kmers_neg = defaultdict(set)
homologous = {}
flag = 0
name = ""
with open(fa_file) as file:
    for line in file:
        line = line.strip()
        if(line[0] == ">"):
            name = line[1:]
        else:
            temp = set()
            for i in range(0, len(line) - 21 + 1):
                temp.add(line[i : i + 21])
            kmers_pos[name] = temp

            line = line.translate(transtable)[::-1]
            temp = set()
            for i in range(0, len(line) - 21 + 1):
                temp.add(line[i : i + 21])
            kmers_neg[name] = temp

for kmer1,kmer2 in itertools.combinations(kmers_pos.items(),2):
    if((kmer1[0],kmer2[0]) in GFA.edges()):
        homologous["".join(sorted([kmer1[0],kmer2[0]])) + " same"] = 0
    else:
        homologous["".join(sorted([kmer1[0],kmer2[0]])) + " same"] = (len(kmer1[1] & kmer2[1]) / min(len(kmer1[1]),len(kmer2[1])))

for key1,key2 in itertools.combinations(kmers_pos.keys(),2):
    if((key1,key2) in GFA.edges()):
        homologous["".join(sorted([key1,key2])) + " diff"] = 0
    else:
        homologous["".join(sorted([key1,key2])) + " diff"] = (len(kmers_pos[key1] & kmers_neg[key2]) / min(len(kmers_pos[key1]),len(kmers_neg[key2])))


class Individual:
    def __init__(self, id, genes, depth, distance, coverage, concord):
        self.id = id
        self.genes = genes
        self.depth = depth
        self.distance = distance
        self.coverage = coverage
        self.concord = concord
        self.adaptivity = 0

class Population:
    def __init__(self, individuals = [], parents = []):
        self.individuals = individuals
        self.parents = parents


Group = Population()
"""if(chr == "chr20"):
    individual = {}
    individual["sfd000001l"] = ['+', 26383741, 28495664, 2111924]
    individual["sfd000012l"] = ['+', 28364123, 29327870, 963748]
    individual["sfd000013l"] = ['+', 26377183, 27204361, 827179]
    individual["sfd000014l"] = ['+', 27594111, 27728487, 134377]
    individual["sfd000015l"] = ['+', 28361675, 28447677, 86003]
    individual["sfd000007l"] = ['-', 28381032, 28585941, 454994]
    individual["sfd000005l"] = ['-', 28347922, 28749648, 754512]
    individual["sfd000016l"] = ['-', 27602543, 28817274, 1214732]
    individual['utg000012l'] = ['-', 28443376, 28773414, 789743]
    concord_flag,concord,new_individual = is_concord(individual)

    if(concord_flag == 1):
        temp_depth,distance,coverage = assess_individual(len(Group.individuals),new_individual) # assess
        individual = Individual(len(Group.individuals),new_individual,temp_depth,distance,coverage,concord)
        print(new_individual,temp_depth,distance,coverage,individual.concord)
        # print()
        Group.individuals.append(individual)"""

loop_count = 0
while(len(Group.individuals) < 500): # randomly generate enough concord individual
    if(loop_count >= gene_banks_combinations_count * 4): # if all the combinations are not valid under given threshold, try to lower the threshold to jump out of endless loop
        global_length_threshold -= 0.2
        if(global_length_threshold <= 0):
            raise Exception("Heterozygosity between reference and query is too high, causing mapping imprecise!\n")
        loop_count = 0
    individual = {}
    for gene in gene_banks.keys():
        individual[gene] = copy.copy(random.choice(gene_banks[gene]))

    concord_flag,concord,new_individual = is_concord(individual)

    if(concord_flag == 0):
        loop_count += 1
        continue
    else:
        loop_count = 0
        temp_depth,distance,coverage = assess_individual(len(Group.individuals),new_individual) # assess
        individual = Individual(len(Group.individuals),new_individual,temp_depth,distance,coverage,concord)
        print(new_individual,temp_depth,distance,coverage,individual.concord)
        print()
        Group.individuals.append(individual)

iter = 0
while(iter < 200):
    print(iter)

    # standardize
    break_flag = 0
    mean = sum([individual.depth for individual in Group.individuals]) / len(Group.individuals)
    std = np.std([individual.depth for individual in Group.individuals])
    if(std == 0):
        for individual in Group.individuals:
            individual.depth = 0
        break_flag = 1
    else:
        for individual in Group.individuals:
            individual.depth = (individual.depth - mean) / std
    mean = sum([individual.distance for individual in Group.individuals]) / len(Group.individuals)
    std = np.std([individual.distance for individual in Group.individuals])
    if(std == 0):
        for individual in Group.individuals:
            individual.distance = 0
        break_flag = 1
    else:
        for individual in Group.individuals:
            individual.distance = (individual.distance - mean) / std
    mean = sum([individual.concord for individual in Group.individuals]) / len(Group.individuals)
    std = np.std([individual.concord for individual in Group.individuals])
    if(std == 0):
        for individual in Group.individuals:
            individual.concord = 0
        break_flag = 1
    else:
        for individual in Group.individuals:
            individual.concord = (individual.concord - mean) / std
    mean = sum([individual.coverage for individual in Group.individuals]) / len(Group.individuals)
    std = np.std([individual.coverage for individual in Group.individuals])
    if(std == 0):
        for individual in Group.individuals:
            individual.coverage = 0
        break_flag = 1
    else:
        for individual in Group.individuals:
            individual.coverage = 2 * (individual.coverage - mean) / std
    for individual in Group.individuals:
        individual.adaptivity = individual.coverage + individual.distance # + individual.depth + individual.concord 

    if(len(set([individual.concord for individual in sorted(Group.individuals,key = lambda individual : individual.adaptivity,reverse = True)[:10]])) == 1): # topN have same concord degree
        for individual in sorted(Group.individuals,key = lambda individual : individual.adaptivity,reverse = True)[:10]:
            print([[key,value[0],value[1],value[2],(value[1] + value[2]) / 2] for key,value in individual.genes.items()],individual.id,individual.depth,individual.distance,individual.coverage,individual.concord,individual.adaptivity)
        print()
        break

    for individual in sorted(Group.individuals,key = lambda individual : individual.adaptivity,reverse = True):
        print([[key,value[0],value[1],value[2],(value[1] + value[2]) / 2] for key,value in individual.genes.items()],individual.id,individual.depth,individual.distance,individual.coverage,individual.concord,individual.adaptivity)
    
    if(break_flag == 1):
        break

    Group.parents = []
    min_score = min([individual.adaptivity for individual in Group.individuals]) - adaptivity_shift
    sum_score = sum([individual.adaptivity - min_score for individual in Group.individuals])
    probability = [(individual.adaptivity - min_score) / sum_score for individual in Group.individuals]
    while(len(Group.parents) < 60): # choose parents according to adaptivity percentage
        temp = random.choices(Group.individuals,probability)[0]
        if(temp not in Group.parents):
            Group.parents.append(temp)

    print()
    for individual in sorted(Group.parents,key = lambda individual : individual.adaptivity,reverse = True):
        print([[key,value[0],value[1],value[2],(value[1] + value[2]) / 2] for key,value in individual.genes.items()],individual.id,individual.depth,individual.distance,individual.coverage,individual.concord,individual.adaptivity)
    

    if(all(len([key for key in set(individual1.genes.keys()) & set(individual2.genes.keys()) if(individual1.genes[key] != individual2.genes[key])]) == 0 for individual1,individual2 in itertools.combinations(Group.parents,2))):
        Group.individuals = Group.parents
        break

    Group.individuals = []
    min_score = min([individual.adaptivity for individual in Group.parents]) - adaptivity_shift
    sum_score = sum([individual.adaptivity - min_score for individual in Group.parents])
    probability = [(individual.adaptivity - min_score) / sum_score for individual in Group.parents]
    while(len(Group.individuals) < 300): # choose parents according to adaptivity percentage to generate offspring
        mat = random.choices(Group.parents,probability)[0]
        pat = random.choices(Group.parents,probability)[0]
        if(mat != pat):
            common_gene = [key for key in set(mat.genes.keys()) & set(pat.genes.keys()) if(mat.genes[key] != pat.genes[key])]
            if(len(common_gene) == 0): # TODO: endless loop sometimes
                continue
            mat_gene = []
            pat_gene = []
            switch_gene = ""
            while(mat_gene == pat_gene): # different gene
                switch_gene = random.choice(list(common_gene)) # switch only one gene
                mat_gene = copy.copy(mat.genes[switch_gene])
                pat_gene = copy.copy(pat.genes[switch_gene])
            child1_genes = copy.deepcopy(mat.genes)
            child2_genes = copy.deepcopy(pat.genes)
            child1_genes[switch_gene] = pat_gene # switch
            child1_genes[switch_gene] = mat_gene # switch
            for child_genes in [child1_genes,child2_genes]:
                concord_flag,concord,new_individual = is_concord(child_genes)
                if(concord_flag):
                    for gene in new_individual.keys(): # mutation
                        if(len(gene_banks[gene]) > 1 and random.random() < 0.05): # there are rooms for mutation, mutation probability = 0.05
                            temp = copy.copy(random.choice(gene_banks[gene]))
                            while(temp == new_individual[gene]): # if mutated, it must mutate to a different location
                                temp = copy.copy(random.choice(gene_banks[gene]))
                            new_individual[gene] = temp
                    for gene in set(gene_banks.keys()) - set(new_individual.keys()): # adding gene
                        if(random.random() < 0.05):
                            new_individual[gene] = copy.copy(random.choice(gene_banks[gene]))
                    concord_flag,concord,new_individual = is_concord(new_individual)
                    if(concord_flag):
                        temp_depth,distance,coverage = assess_individual(len(Group.individuals),new_individual)
                        individual = Individual(len(Group.individuals),new_individual,temp_depth,distance,coverage,concord)
                        print(new_individual,temp_depth,distance,coverage,individual.concord)
                        print()
                        Group.individuals.append(individual)
        
    iter += 1

# print()
origin_score = [0] * 10
for i,individual in enumerate(sorted(Group.individuals,key = lambda individual : individual.adaptivity,reverse = True)[:10]):
    for key,value in individual.genes.items():
        if(any(gene[1] == value[1] and gene[2] == value[2] for gene in gene_banks[key])):
            origin_score[i] += 1
max_index = origin_score.index(max(origin_score))

writer = open(output_combination_file,"w")
goal_combination = sorted(Group.individuals,key = lambda individual : individual.adaptivity,reverse = True)[max_index].genes
goal_combination = {key : value for key,value in goal_combination.items() if(max(value[1],gap_start) <= min(value[2],gap_end))} # delete those lie totally outside of the gap

if(ploid == "haploid"):
    goal_combination = {key : value1 for key,value1 in goal_combination.items() if(all(value2[1] > value1[1] or value1[2] > value2[2] for value2 in goal_combination.values() if(value1 != value2)))} # delete those surrounded by others

    G = nx.DiGraph()
    for node in sorted(goal_combination.keys(), key = lambda node : goal_combination[node][1]):
        G.add_node(node)
        for temp_node in G.nodes():
            if(temp_node != node and max(goal_combination[node][1],goal_combination[temp_node][1]) <= min(goal_combination[node][2],goal_combination[temp_node][2])):
                G.add_edge(temp_node,node)
    sources = [node for node in G.nodes() if(G.in_degree[node] == 0)]
    connected_components = []
    for node in sources:
        if(any(node in component for component in connected_components)):
            continue
        temp = nx.bfs_tree(G,node).nodes()
        connected_components.append(list(temp))
    for node in G.nodes():
        if(G.in_degree[node] == 0 and G.out_degree[node] == 0):
            G.nodes[node]["isolated"] = 1
        else:
            G.nodes[node]["isolated"] = 0
    for (src,dst) in list(G.edges()): # remove edges that overlap too much with others
        if(min(goal_combination[src][2],goal_combination[dst][2]) - max(goal_combination[src][1],goal_combination[dst][1]) + 1 > 0.9 * min(GFA.nodes[src]["length"],GFA.nodes[dst]["length"])):
            print(src,dst)
            G.remove_edge(src,dst)
    deleted_node = set()
    for node in list(G.nodes()): # remove nodes whose edges are depleted in last step
        if(G.in_degree[node] == 0 and G.out_degree[node] == 0 and G.nodes[node]["isolated"] == 0):
            G.remove_node(node)
            deleted_node.add(node)
    for component in connected_components: # in case that some components are depleted
        if(all(node not in G.nodes() for node in component)):
            G.add_node(sorted(deleted_node & set(component), key = lambda node : GFA.nodes[node]["length"], reverse = True)[0])

    sources = [node for node in G.nodes() if(G.in_degree[node] == 0)]
    connected_components = []
    for node in sources:
        if(any(node in component for component in connected_components)):
            continue
        temp = nx.bfs_tree(G,node).nodes()
        connected_components.append(list(temp))
    final_sequence = ""
    temp = sorted(connected_components, key = lambda component : min([goal_combination[node][1] for node in component]))
    for i,component in enumerate(temp):
        if(len(component) == 1):
            final_sequence += GFA.nodes[component[0]]["seq"]
            print((" ").join([component[0]] + list(map(str,goal_combination[component[0]][ : 3]))), file = writer)
        else:
            src = sorted(component, key = lambda node : goal_combination[node][1])[0]
            dst = sorted(component, key = lambda node : goal_combination[node][1])[-1]
            path = sorted([path for path in nx.all_shortest_paths(G,src,dst)], key = lambda path : sum([GFA.nodes[node]["length"] for node in path]), reverse = True)[0] # if several, pick the longest
            print(path)
            for j,node in enumerate(path):
                print((" ").join([node] + list(map(str,goal_combination[node][ : 3]))), file = writer)
                final_sequence += GFA.nodes[node]["seq"]
                if(j != len(path) - 1):
                    final_sequence += 'N' * 25
        if(i != len(temp) - 1):
            print(file = writer)
            final_sequence += 'N' * 50
    writer.close()

    writer = open(output_fa_file,"w")
    writer.write(">cen000001l\n" + final_sequence + "\n")
    writer.close()
elif(ploid == "diploid"):
    for key,value in goal_combination.items():
        writer.write((" ").join([key] + list(map(str,value[ : 3]))) + "\n")
    writer.close()

    writer = open(output_fa_file,"w")
    flag = 0
    with open(fa_file) as file:
        for line in file:
            line = line.strip()
            if(line[0] == ">" and line[1:] in goal_combination):
                writer.write(line + "\n")
                flag = 1
            elif(flag):
                writer.write(line + "\n")
                flag = 0
    writer.close()
else:
    raise Exception("Ploid option must be haploid or diploid!\n")
