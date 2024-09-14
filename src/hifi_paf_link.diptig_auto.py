#!/usr/bin/env python
from collections import defaultdict,Counter
from Bio import pairwise2
import matplotlib.pyplot as plt
import networkx as nx
import itertools
import copy
import sys

fa_file = sys.argv[1]
gfa_file = sys.argv[2]
paf_file = sys.argv[3]
output_file = sys.argv[4]

pairs_matchs = defaultdict(list)
raw_pairs_matchs = defaultdict(list)
transtable = str.maketrans('ATGC','TACG')
trans_strand = str.maketrans('+-',"-+")

def pair_process(query_name,contig1,info1,contig2,info2,GFA):
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
    assert query_len1 == query_len2
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
        if(match < query_len1 * 3):
            if(contig2 + strand2.translate(trans_strand) + contig1 + strand1.translate(trans_strand) in pairs_matchs):
                pairs_matchs[contig2 + strand2.translate(trans_strand) + contig1 + strand1.translate(trans_strand)].append(match)
            elif((contig1,contig2) in GFA.edges() and GFA[contig1][contig2]['strand1'] == strand1 and GFA[contig1][contig2]['strand2'] == strand2):
                pairs_matchs[contig1 + strand1 + contig2 + strand2].append(match)
            elif(contig2 + strand2.translate(trans_strand) + contig1 + strand1.translate(trans_strand) in raw_pairs_matchs):
                raw_pairs_matchs[contig2 + strand2.translate(trans_strand) + contig1 + strand1.translate(trans_strand)].append(match)
            else:
                raw_pairs_matchs[contig1 + strand1 + contig2 + strand2].append(match)
            #print(query_len1,contig1,strand1,contig2,strand2,match,1)
    elif(query_start1 - ref_start1 > query_start2 - ref_start2 and query_end1 + ref_len1 - ref_end1 > query_end2 + ref_len2 - ref_end2):
        match = (query_end2 + ref_len2 - ref_end2) - (query_start1 - ref_start1)
        if(match < query_len1 * 3):
            if(contig1 + strand1.translate(trans_strand) + contig2 + strand2.translate(trans_strand) in pairs_matchs):
                pairs_matchs[contig1 + strand1.translate(trans_strand) + contig2 + strand2.translate(trans_strand)].append(match)
            elif((contig2,contig1) in GFA.edges() and GFA[contig2][contig1]['strand1'] == strand2 and GFA[contig2][contig1]['strand2'] == strand1):
                pairs_matchs[contig2 + strand2 + contig1 + strand1].append(match)
            elif(contig1 + strand1.translate(trans_strand) + contig2 + strand2.translate(trans_strand) in raw_pairs_matchs):
                raw_pairs_matchs[contig1 + strand1.translate(trans_strand) + contig2 + strand2.translate(trans_strand)].append(match)
            else:
                raw_pairs_matchs[contig2 + strand2 + contig1 + strand1].append(match)
            #print(query_len1,contig2,strand2,contig1,strand1,match,2)
    #else:
        #print(query_name,contig1,strand1,contig2,strand2,"contained")

def homologous_seq(seqs):
    kmers = []
    flag = 0
    with open(fa_file) as file:
        for line in file:
            line = line.strip()
            if(line[1:] in seqs):
                if(seqs[line[1:]] == "+"):
                    flag = 1
                else:
                    flag = 2
            elif(flag == 1):
                temp = set()
                for i in range(0, len(line) - 21 + 1):
                    temp.add(line[i : i + 21])
                kmers.append(temp)
                flag = 0
            elif(flag == 2):
                line = line.translate(transtable)[::-1]
                temp = set()
                for i in range(0, len(line) - 21 + 1):
                    temp.add(line[i : i + 21])
                kmers.append(temp)
                flag = 0
    for (kmer1,kmer2) in itertools.combinations(kmers,2):
        if(len(kmer1 & kmer2) / min(len(kmer1),len(kmer2)) >= 0.8):
            return 1
    return 0

def structure_detect(node,target,next_hops,prior): # small bubble detection
    for next_hop in next_hops.keys():# see if target can reach next_hop through 0/1/2/3 middle node
        if(next_hop != target):
            if((target,next_hop) in GFA.edges() and GFA[node][next_hop]["strand2"] == GFA[target][next_hop]["strand2"]): # 0 middle node
                return 1
            else: # 1 middle node
                for middle_node1 in GFA[target].keys():
                    if(middle_node1 != node and GFA[node][target]["strand2"] == GFA[target][middle_node1]["strand1"] and (middle_node1,next_hop) in GFA.edges() and GFA[middle_node1][next_hop]["strand2"] == GFA[node][next_hop]["strand2"]):
                        prior[middle_node1] += 1
                        return 1
                    else: # 2 middle node
                        for middle_node2 in GFA[middle_node1].keys():
                            if(middle_node2 != node and middle_node2 != target and GFA[node][target]["strand2"] == GFA[target][middle_node1]["strand1"] and GFA[target][middle_node1]["strand2"] == GFA[middle_node1][middle_node2]["strand1"] and (middle_node2,next_hop) in GFA.edges() and GFA[middle_node2][next_hop]["strand2"] == GFA[node][next_hop]["strand2"]):
                                prior[middle_node1] += 2
                                prior[middle_node2] += 1
                                return 1
                            else: # 3 middle node
                                for middle_node3 in GFA[middle_node2].keys():
                                    if(GFA[node][target]["strand2"] == GFA[target][middle_node1]["strand1"] and GFA[target][middle_node1]["strand2"] == GFA[middle_node1][middle_node2]["strand1"] and GFA[middle_node1][middle_node2]["strand2"] == GFA[middle_node2][middle_node3]["strand1"] and (middle_node3,next_hop) in GFA.edges() and GFA[middle_node3][next_hop]["strand2"] == GFA[node][next_hop]["strand2"]):
                                        prior[middle_node1] += 3
                                        prior[middle_node2] += 2
                                        prior[middle_node3] += 1
                                        return 1
    return 0

def loop_next_hop(src,next_hop): # would next_hop induce loop?
    temp = []
    count = 0
    for node in GFA[next_hop].keys():
        if(node != src and GFA[src][next_hop]["strand2"] == GFA[next_hop][node]["strand1"]):
            count += 1
            if(len(GFA[node].keys()) // 2 > 1 and GFA.nodes[node]["visit_count"] <= 1):
                temp.append(node)
            elif(len(GFA[node].keys()) // 2 <= 1 and GFA.nodes[node]["visit_count"] <= 0):
                temp.append(node)
    if(len(temp) > 0 and len(temp) == count):
        return 1
    else:
        return 0

def get_source_node():
    tips = set()
    queue = set()
    for edge in GFA.edges(): # get those whose indegree is 0, ignore isolate node(recollect them later)
        src = edge[0]
        dst = edge[1]
        flag = 0
        for node in GFA[dst].keys():
            if(GFA[src][dst]["strand2"] == GFA[dst][node]["strand1"]):
                flag = 1
                break
        if(flag == 0):
            queue.add(dst)

    connected_components = []
    for node in queue:
        if(any(node in component for component in connected_components)):
            continue
        temp = nx.bfs_tree(GFA,node).nodes()
        connected_components.append(temp)

    homologous = {}
    for component in connected_components:
        """if(sum([GFA.nodes[node]["length"] for node in component]) < 100000): # control connected component size # 300000
            continue"""
        temp = queue & component
        for (node1,node2) in itertools.combinations(temp,2):
            if(homologous_seq({node1:GFA[node1][list(GFA[node1].keys())[0]]["strand1"],node2:GFA[node2][list(GFA[node2].keys())[0]]["strand1"]})):
                homologous[node1] = node2
                homologous[node2] = node1

    temp = copy.copy(queue) # remove tip
    for node1 in temp:
        for node2 in GFA[node1].keys():
            if(GFA.nodes[node1]["length"] < 40000 and any(node3 != node1 and node3 not in GFA[node1].keys() and GFA[node3][node2]["strand2"] == GFA[node1][node2]["strand2"] for node3 in GFA[node2].keys())): # 25000
                queue.remove(node1)
                tips.add(node1)
                break
    print(tips)

    for tip in list(tips): # resume homologous
        if(tip in tips and tip in homologous and homologous[tip] in tips):
            tips.remove(tip)
            tips.remove(homologous[tip])
            queue.add(tip)
            queue.add(homologous[tip])

    temp = copy.copy(tips)
    for (node1,node2) in itertools.combinations(temp,2):
        if(node1 in tips and node2 in tips and list(GFA[node1].keys())[0] == list(GFA[node2].keys())[0]):
            tips.remove(max([node1,node2],key = lambda node : GFA.nodes[node]["coverage"]))

    print(tips)
    for tip in list(tips):
        for dst in [node for node in GFA[tip].keys()]:
            for src in [node for node in GFA[dst].keys() if(node != tip and GFA[tip][dst]["strand2"] == GFA[node][dst]["strand2"])]:
                if(len(GFA[tip].keys()) == 1 and (src + GFA[src][dst]["strand1"] + tip + GFA[tip][dst]["strand1"] in raw_pairs_matchs or tip + GFA[tip][dst]["strand1"].translate(trans_strand) + src + GFA[src][dst]["strand1"].translate(trans_strand) in raw_pairs_matchs)):
                    match = 0
                    max_rate = 0
                    for temp_match in raw_pairs_matchs[src + GFA[src][dst]["strand1"] + tip + GFA[tip][dst]["strand1"]]:
                        rate = test_match(src,GFA[src][dst]["strand1"],tip,GFA[tip][dst]["strand1"],temp_match)
                        if(rate > max_rate):
                            match = temp_match
                            max_rate = rate
                    for temp_match in raw_pairs_matchs[tip + GFA[tip][dst]["strand1"].translate(trans_strand) + src + GFA[src][dst]["strand1"].translate(trans_strand)]:
                        rate = test_match(tip,GFA[tip][dst]["strand1"].translate(trans_strand),src,GFA[src][dst]["strand1"].translate(trans_strand),temp_match)
                        if(rate > max_rate):
                            match = temp_match
                            max_rate = rate
                    if(match <= 30000):
                        if(src + GFA[src][dst]["strand1"] + dst + GFA[src][dst]["strand2"] in pairs_matchs):
                            pairs_matchs.pop(src + GFA[src][dst]["strand1"] + dst + GFA[src][dst]["strand2"])
                            pairs_match.pop(src + GFA[src][dst]["strand1"] + dst + GFA[src][dst]["strand2"])
                        elif(dst + GFA[src][dst]["strand2"].translate(trans_strand) + src + GFA[src][dst]["strand1"].translate(trans_strand) in pairs_matchs):
                            pairs_matchs.pop(dst + GFA[src][dst]["strand2"].translate(trans_strand) + src + GFA[src][dst]["strand1"].translate(trans_strand))
                            pairs_match.pop(dst + GFA[src][dst]["strand2"].translate(trans_strand) + src + GFA[src][dst]["strand1"].translate(trans_strand))
                        pairs_matchs[src + GFA[src][dst]["strand1"] + tip + GFA[tip][dst]["strand1"]] = raw_pairs_matchs[src + GFA[src][dst]["strand1"] + tip + GFA[tip][dst]["strand1"]] + raw_pairs_matchs[tip + GFA[tip][dst]["strand1"].translate(trans_strand) + src + GFA[src][dst]["strand1"].translate(trans_strand)]
                        pairs_match[src + GFA[src][dst]["strand1"] + tip + GFA[tip][dst]["strand1"]] = raw_pairs_match[src + GFA[src][dst]["strand1"] + tip + GFA[tip][dst]["strand1"]] + raw_pairs_match[tip + GFA[tip][dst]["strand1"].translate(trans_strand) + src + GFA[src][dst]["strand1"].translate(trans_strand)]
                        print(src + GFA[src][dst]["strand1"] + tip + GFA[tip][dst]["strand1"],match,sorted(pairs_matchs[src + GFA[src][dst]["strand1"] + tip + GFA[tip][dst]["strand1"]]))
                        strand1 = GFA[src][dst]["strand1"]
                        GFA.remove_edge(src,dst)
                        GFA.remove_edge(dst,src)
                        GFA.add_edge(src,tip,strand1 = strand1,strand2 = GFA[tip][dst]["strand1"],match = match)
                        GFA.add_edge(tip,src,strand1 = GFA[tip][dst]["strand1"].translate(trans_strand),strand2 = strand1.translate(trans_strand),match = match)
                        tips.remove(tip)

    result = set()
    for component in connected_components: # leave as many source node as possible in a connected component, in order to get all the route
        """if(sum([GFA.nodes[node]["length"] for node in component]) < 100000): # control connected component size # 200000
            continue"""
        temp = queue & set(component)
        if(len(temp & homologous.keys()) != 0):
            result |= temp & homologous.keys()
        if(len(temp) != 0):
            result = result | temp
        else:
            result |= tips & set(component)
            for tip in tips & set(component):
                tips.remove(tip)
    
    for tip in tips:
        GFA.remove_node(tip)

    result = [GFA[node][list(GFA[node].keys())[0]]["strand1"] + node for node in result if(len(GFA[node].keys()) > 0)]
    result.sort(key = lambda node : (GFA.nodes[node[1:]]["coverage"],node[1:]),reverse = True)

    return result,homologous

def test_match(src,strand1,dst,strand2,match):
    src_seq = ""
    dst_seq = ""
    if(strand1 == "+"):
        src_seq = GFA.nodes[src]["seq"][GFA.nodes[src]["length"] - match : ]
    else:
        src_seq = GFA.nodes[src]["seq"].translate(transtable)[::-1][GFA.nodes[src]["length"] - match : ]
    if(strand2 == "+"):
        dst_seq = GFA.nodes[dst]["seq"][ : match]
    else:
        dst_seq = GFA.nodes[dst]["seq"].translate(transtable)[::-1][ : match]

    alignments = pairwise2.align.globalxx(src_seq,dst_seq) # type: ignore
    max_rate = 0
    for alignment in alignments:
        if(alignment.score / match > max_rate):
            max_rate = alignment.score / match
        # print(src,strand1,dst,strand2,match,alignment.score)
    return max_rate

def BFS(queue,G):
    print(queue)
    visited = set()
    prior = {node : 0 for node in GFA.nodes()}
    next_hop = ""
    step = 0
    while(len(queue) != 0):
        strand = queue[0][0]
        node = queue[0][1:]
        if(node == "utg000130l"):
            temp = 1
        # print(node)
        queue.pop(0)
        GFA.nodes[node]["visit_count"] -= 1
        if(GFA.nodes[node]["visit_count"] == 0):
            visited.add(node)
        elif(GFA.nodes[node]["visit_count"] <= -2):
            break
        while("+" + node in queue):
            queue.remove("+" + node)
        while("-" + node in queue):
            queue.remove("-" + node)    
        next_hops = {}
        next_hops_available = {}
        for next_hop in GFA[node].keys():
            if(GFA[node][next_hop]["strand1"] == strand):
                next_hops[next_hop] = GFA[node][next_hop]["strand2"]
                if(next_hop not in visited):
                    next_hops_available[next_hop] = GFA[node][next_hop]["strand2"]
        small_bubble_flag = 0
        for next_hop in sorted(next_hops.keys(),key = lambda next_hop: GFA.nodes[next_hop]["length"],reverse = True):
            if(not loop_next_hop(node,next_hop) and structure_detect(node,next_hop,next_hops,prior)): # small bubble
                small_bubble_flag = 1
                G.add_edge(node,next_hop,strand1 = strand,strand2 = GFA[node][next_hop]["strand2"],match = GFA[node][next_hop]["match"],way = 1,step = step)
                queue.append(GFA[node][next_hop]["strand2"] + next_hop)
                break
        if(small_bubble_flag == 0):
            if(len(next_hops_available) == 1): # only one
                next_hop = list(next_hops_available.keys())[0]
                G.add_edge(node,next_hop,strand1 = strand,strand2 = GFA[node][next_hop]["strand2"],match = GFA[node][next_hop]["match"],way = 2,step = step)
                queue.append(GFA[node][next_hop]["strand2"] + next_hop)
            else:
                next_hops_available_unique = {}
                for next_hop in next_hops_available.keys():
                    if(all(GFA[neighbor][next_hop]["strand2"] != GFA[node][next_hop]["strand2"] for neighbor in GFA[next_hop].keys() if(neighbor != node)) and not loop_next_hop(node,next_hop)): # only node can enter next_hop
                        next_hops_available_unique[next_hop] = GFA[node][next_hop]["strand2"]
                if(len(next_hops_available_unique) == 1): # choose unique path first, in order to use more node, while try not to induce loop
                    next_hop = list(next_hops_available_unique.keys())[0]
                    G.add_edge(node,next_hop,strand1 = strand,strand2 = GFA[node][next_hop]["strand2"],match = GFA[node][next_hop]["match"],way = 3,step = step)
                    queue.append(GFA[node][next_hop]["strand2"] + next_hop)
                else:
                    small_bubble_prior_flag = 0
                    if(any(prior[next_hop] > 0 for next_hop in next_hops_available.keys())): # choose small bubble's path first
                        small_bubble_prior_flag = 1
                        next_hop = max(next_hops_available.keys(),key = lambda next_hop: prior[next_hop])
                        G.add_edge(node,next_hop,strand1 = strand,strand2 = GFA[node][next_hop]["strand2"],match = GFA[node][next_hop]["match"],way = 4,step = step)
                        queue.append(GFA[node][next_hop]["strand2"] + next_hop)
                    if(small_bubble_prior_flag == 0):
                        trustable_flag = 0
                        for next_hop in next_hops_available.keys(): # keep all trustable but not loop-inducing node
                            if(GFA.nodes[next_hop]["coverage"] >= middle_coverage and not loop_next_hop(node,next_hop)):
                                G.add_edge(node,next_hop,strand1 = strand,strand2 = GFA[node][next_hop]["strand2"],match = GFA[node][next_hop]["match"],way = 5,step = step)
                                queue.append(GFA[node][next_hop]["strand2"] + next_hop)
                                trustable_flag = 1
                        if(trustable_flag == 0):
                            out_perform_flag = 0
                            for node1,node2 in itertools.combinations(next_hops_available.keys(),2): # keep outperform node if exist
                                if(abs(GFA.nodes[node1]["coverage"] - GFA.nodes[node2]["coverage"]) >= 3):
                                    next_hop = max([node1,node2],key = lambda node : GFA.nodes[node]["coverage"])
                                    G.add_edge(node,next_hop,strand1 = strand,strand2 = GFA[node][next_hop]["strand2"],match = GFA[node][next_hop]["match"],way = 6,step = step)
                                    queue.append(GFA[node][next_hop]["strand2"] + next_hop)
                                    out_perform_flag = 1
                            if(out_perform_flag == 0):
                                recommend_outperform_flag = 0
                                recommend = {}
                                for next_hop in next_hops_available.keys():
                                    if(node + GFA[node][next_hop]["strand1"] + next_hop + GFA[node][next_hop]["strand2"] in pairs_match):
                                        recommend[next_hop] = len(pairs_match[GFA[node][next_hop]["strand1"] + node + GFA[node][next_hop]["strand2"] + next_hop])
                                    elif(next_hop + GFA[node][next_hop]["strand2"].translate(trans_strand) + node + GFA[node][next_hop]["strand1"].translate(trans_strand) in pairs_match):
                                        recommend[next_hop] = len(pairs_match[GFA[node][next_hop]["strand2"].translate(trans_strand) + next_hop + GFA[node][next_hop]["strand1"].translate(trans_strand) + node])
                                if(len(recommend) > 0):
                                    max_recommend = max(recommend.keys(),key = lambda key : recommend[key])
                                    if(all(recommend[next_hop] < recommend[max_recommend] for next_hop in recommend.keys() if(next_hop != max_recommend))): # recommend outperform
                                        G.add_edge(node,max_recommend,strand1 = strand,strand2 = GFA[node][max_recommend]["strand2"],match = GFA[node][max_recommend]["match"],way = 7,step = step)
                                        queue.append(GFA[node][max_recommend]["strand2"] + max_recommend)
                                        recommend_outperform_flag = 1
                                if(recommend_outperform_flag == 0):
                                    for next_hop in next_hops_available.keys(): # keep all
                                        G.add_edge(node,next_hop,strand1 = strand,strand2 = GFA[node][next_hop]["strand2"],match = GFA[node][next_hop]["match"],way = 8,step = step)
                                        queue.append(GFA[node][next_hop]["strand2"] + next_hop)
        step += 1

def DFS(source,result,stack):
    node = source[1:]
    strand = source[0]

    if(G_global.in_degree[node] > 1 and G_global.out_degree[node] > 1): # branch
        results.append(result)
        if(G_global.nodes[node]["visit_count"] > 1):
            results.append([source])
        result = []
    elif(G_global.in_degree[node] > 1):
        results.append(result)
        result = [source]
    elif(G_global.out_degree[node] > 1):
        # get sequence in result with node
        result.append(source)
        results.append(result)
        result = []
    else:
        result.append(source)

    if(G_global.out_degree[node] == 0):
        results.append(result)
        result = []

    G_global.nodes[node]["visit_count"] -= 1
    if(G_global.nodes[node]["visit_count"] == 0):
        visited.add(node)
        if(node in stack and [source] in results): # loop
            index = results.index([source])
            index1 = index - 1
            index2 = index + 1
            flag1 = 0
            flag2 = 0
            while(index1 >= 0 and results[index1] == []):
                index1 -= 1
            while(index2 < len(results) and results[index2] == []):
                index2 += 1
            result = []
            if((results[index1][-1][1:],node) in G_global.edges() and G_global[results[index1][-1][1:]][node]["strand2"] == strand):
                result += results[index1]
                result += results[index]
                flag1 = 1
            if((node,results[index2][0][1:]) in G_global.edges() and G_global[node][results[index2][0][1:]]["strand1"] == strand):
                result += results[index2]
                flag2 = 1
            if(index2 == len(results) - 1):
                if(flag2):
                    results.pop(index2)
                results.pop(index)
                if(flag1):
                    results.pop(index1)
                result.append(source)
            else:
                if(flag2):
                    results.pop(index2)
                results.pop(index)
                if(flag1):
                    results.pop(index1)
                    results.insert(index1,result)
                else:
                    results.insert(index,result)
                index = len(results) - 1
                while(index >= 0 and results[index] == []):
                    index -= 1
                if((results[index][-1][1:],node) in G_global.edges() and G_global[results[index][-1][1:]][node]["strand2"] == strand):
                    result = results[index]
                    results.pop(index)
                    result.append(source)
                else:
                    result = [source]

    next_hops_available = [next_hop for next_hop in G_global[node].keys() if(next_hop not in visited)]
    if(len(next_hops_available) == 0):
        results.append(result)
        result = []
        return
    else:
        stack.append(node)
        for neighbor in sorted(next_hops_available,key = lambda next_hop: G_global[node][next_hop]["step"]):
            if(neighbor not in visited):
                strand = G_global[node][neighbor]["strand2"]
                DFS(strand + neighbor,result,stack)
                result = []
        stack.pop(-1)

GFA = nx.MultiDiGraph()
with open(gfa_file) as file:
    for line in file:
        line = line.strip().split()
        if(line[0] == "S"):
            GFA.add_node(line[1],length = len(line[2]),seq = line[2],coverage = int(line[4][5:]),visit_count = 1)
        elif(line[0] == "L" and line[1] != line[3]): # no self loop
            GFA.add_edge(line[1],line[3],strand1 = line[2],strand2 = line[4],match = int(line[5][:-1]))
coverage = [GFA.nodes[node]["coverage"] for node in GFA.nodes()]
middle_coverage = sorted(coverage)[len(coverage) // 2]
print(f"middle_coverage = {middle_coverage}")

loop_edge = [key for key,value in Counter(GFA.edges()).items() if value > 1]
temp_loop_edge = []
for edge in loop_edge:
    temp = sorted(list(edge))
    if(temp not in temp_loop_edge):
        temp_loop_edge.append(temp)
loop_edge = temp_loop_edge
for edge in loop_edge: # remove small loop
    src = edge[0]
    dst = edge[1]
    # print(src,dst)
    if((src,dst) in GFA.edges() and len(GFA[src][dst]) > 1):
        up = ""
        down = ""
        if(len(GFA[src].keys()) < len(GFA[dst].keys())):
            up = src
            down = dst
        elif(len(GFA[src].keys()) > len(GFA[dst].keys())):
            up = dst
            down = src
        else:
            continue
        if(GFA.nodes[up]["length"] <= GFA.nodes[down]["length"] or (GFA.nodes[up]["length"] > GFA.nodes[down]["length"] and GFA.nodes[up]["length"] - GFA.nodes[down]["length"] < 100)): # up side shorter,or close length,keep down side
            GFA.remove_node(up)
        else:
            length = GFA.nodes[down]["length"]
            flag = 0
            if(len(GFA[down].keys()) == 2):
                for node in GFA[down].keys():
                    if(node != up):
                        src = node
                        dst = node
            else:
                for node in GFA[down].keys():
                    if(node != up and flag == 0):
                        src = node
                        flag = 1
                        continue
                    if(node != up and flag == 1):
                        dst = node
                        break
            match_src_down = GFA[src][down][0]["match"]
            match_dst_down = GFA[down][dst][0]["match"]
            match_src_down_up = 0
            match_dst_down_up = 0
            strand_up = ""
            if(GFA[src][down][0]["strand2"] == GFA[down][up][0]["strand1"]):
                match_src_down_up = GFA[down][up][0]["match"]
                match_dst_down_up = GFA[down][up][1]["match"]
                strand_up = GFA[down][up][0]["strand2"]
            else:
                match_src_down_up = GFA[down][up][1]["match"]
                match_dst_down_up = GFA[down][up][0]["match"]
                strand_up = GFA[down][up][1]["strand2"]
            match_src_down = match_src_down - (length - match_src_down_up)
            match_dst_down = match_dst_down - (length - match_dst_down_up)
            if(match_src_down > 0 and match_dst_down > 0): # both side changed
                GFA.add_edge(src,up,strand1 = GFA[src][down][0]["strand1"],strand2 = strand_up,match = match_src_down)
                GFA.add_edge(up,src,strand1 = strand_up.translate(trans_strand),strand2 = GFA[src][down][0]["strand1"].translate(trans_strand),match = match_src_down)
                GFA.add_edge(up,dst,strand1 = strand_up,strand2 = GFA[down][dst][0]["strand2"],match = match_dst_down)
                GFA.add_edge(dst,up,strand1 = GFA[down][dst][0]["strand2"].translate(trans_strand),strand2 = strand_up.translate(trans_strand),match = match_dst_down)
                GFA.remove_node(down)
            elif(match_src_down <= 0): # src side unchanged
                strand_src = GFA[src][down][0]["strand1"]
                strand_down = GFA[src][down][0]["strand2"]
                strand_dst = GFA[down][dst][0]["strand2"]
                match_src_down = GFA[src][down][0]["match"]
                length = GFA.nodes[down]["length"]
                seq = GFA.nodes[down]["seq"]
                coverage = GFA.nodes[down]["coverage"]
                visit_count = GFA.nodes[down]["visit_count"]
                GFA.remove_node(down)
                GFA.add_node(down,length = length,seq = seq,coverage = coverage,visit_count = visit_count)
                GFA.add_edge(src,down,strand1 = strand_src,strand2 = strand_down,match = match_src_down)
                GFA.add_edge(down,src,strand1 = strand_down.translate(trans_strand),strand2 = strand_src.translate(trans_strand),match = match_src_down)
                GFA.add_edge(down,up,strand1 = strand_down,strand2 = strand_up,match = match_src_down_up)
                GFA.add_edge(up,down,strand1 = strand_up.translate(trans_strand),strand2 = strand_down.translate(trans_strand),match = match_src_down_up)
                GFA.add_edge(up,dst,strand1 = strand_up,strand2 = strand_dst,match = match_dst_down)
                GFA.add_edge(dst,up,strand1 = strand_dst.translate(trans_strand),strand2 = strand_up.translate(trans_strand),match = match_dst_down)
            elif(match_dst_down <= 0): # dst side unchanged
                strand_src = GFA[src][down][0]["strand1"]
                strand_down = GFA[down][dst][0]["strand1"]
                strand_dst = GFA[down][dst][0]["strand2"]
                match_dst_down = GFA[down][dst][0]["match"]
                length = GFA.nodes[down]["length"]
                seq = GFA.nodes[down]["seq"]
                coverage = GFA.nodes[down]["coverage"]
                visit_count = GFA.nodes[down]["visit_count"]
                GFA.remove_node(down)
                GFA.add_node(down,length = length,seq = seq,coverage = coverage,visit_count = visit_count)
                GFA.add_edge(src,up,strand1 = strand_src,strand2 = strand_up,match = match_src_down)
                GFA.add_edge(up,src,strand1 = strand_up.translate(trans_strand),strand2 = strand_src.translate(trans_strand),match = match_src_down)
                GFA.add_edge(up,down,strand1 = strand_up,strand2 = strand_down,match = match_dst_down_up)
                GFA.add_edge(down,up,strand1 = strand_down.translate(trans_strand),strand2 = strand_up.translate(trans_strand),match = match_dst_down_up)
                GFA.add_edge(down,dst,strand1 = strand_down,strand2 = strand_dst,match = match_dst_down)
                GFA.add_edge(dst,down,strand1 = strand_dst.translate(trans_strand),strand2 = strand_down.translate(trans_strand),match = match_dst_down)
            else: # both side unchanged
                raise Exception("small loop process failed!\n")
GFA = nx.DiGraph(GFA)

PAF = defaultdict(list)
name = ""
with open(paf_file) as file:
    for line in file:
        line = line.strip().split()
        if(line[0] != name and name != ""):
            # process
            for key1,key2 in itertools.combinations(PAF.keys(),2):
                for value1 in PAF[key1]:
                    for value2 in PAF[key2]:
                        pair_process(name,key1,value1,key2,value2,GFA)
            PAF.clear()
        name = line[0]
        PAF[line[5]].append([int(line[1]),int(line[2]),int(line[3]),line[4],int(line[6]),int(line[7]),int(line[8])])
    for key1,key2 in itertools.combinations(PAF.keys(),2):
        for value1 in PAF[key1]:
            for value2 in PAF[key2]:
                pair_process(name,key1,value1,key2,value2,GFA)
    PAF.clear()

pairs_match = defaultdict(list)
for key,values in pairs_matchs.items():
    key1 = key[0:10]
    key2 = key[11:21]
    buck = defaultdict(list)
    for match in values:
        buck[match // 10].append(match)
    tempbuck = copy.copy(buck)
    max_key = max(tempbuck.keys(),key = lambda key : len(buck[key] + buck[key + 1] + buck[key - 1]))
    temp = buck[max_key - 1] + buck[max_key] + buck[max_key + 1]
    temp.sort()
    pairs_match[key] = temp
    print(key,sorted(values))
    # print(key,sorted(values),temp[len(temp) // 2])
raw_pairs_match = defaultdict(list)
for key,values in raw_pairs_matchs.items():
    key1 = key[0:10]
    key2 = key[11:21]
    buck = defaultdict(list)
    for match in values:
        buck[match // 10].append(match)
    tempbuck = copy.copy(buck)
    max_key = max(tempbuck.keys(),key = lambda key : len(buck[key] + buck[key + 1] + buck[key - 1]))
    temp = buck[max_key - 1] + buck[max_key] + buck[max_key + 1]
    temp.sort()
    raw_pairs_match[key] = temp

for node in GFA.nodes(): # resolve small bubble in original GFA to facilitate BFS by alleviating degree effection
    neighbors = list(GFA[node].keys())
    if(len(neighbors) == 2 and (neighbors[0],neighbors[1]) in GFA.edges() and GFA[neighbors[0]][node]["strand2"] == GFA[node][neighbors[1]]["strand1"] and GFA[neighbors[0]][node]["strand1"] == GFA[neighbors[0]][neighbors[1]]["strand1"] and GFA[node][neighbors[1]]["strand2"] == GFA[neighbors[0]][neighbors[1]]["strand2"]):
        GFA.remove_edge(neighbors[0],neighbors[1])
        GFA.remove_edge(neighbors[1],neighbors[0])

G_global = nx.DiGraph()
sources,homologous = get_source_node()
for source in sources:
    for node in GFA.nodes():
        if(len(GFA[node].keys()) // 2 > 1):
            GFA.nodes[node]["visit_count"] = len(GFA[node]) // 2
        else:
            GFA.nodes[node]["visit_count"] = 1
    G = nx.DiGraph()
    BFS([source],G)
    
    same_strand = [edge for edge in G.edges() if((edge[0],edge[1]) in G_global.edges())]
    diff_strand = [edge for edge in G.edges() if((edge[1],edge[0]) in G_global.edges())]
    relative_strand = "+" if(len(same_strand) >= len(diff_strand)) else "-"
    for edge in G.edges():
        src = edge[0]
        dst = edge[1]
        if((src,dst) not in G_global.edges() and (dst,src) not in G_global.edges()):
            if(relative_strand == "+"):
                G_global.add_edge(edge[0],edge[1],strand1 = G[edge[0]][edge[1]]["strand1"],strand2 = G[edge[0]][edge[1]]["strand2"],match = G[edge[0]][edge[1]]["match"],way = G[edge[0]][edge[1]]["way"],step = -1) # G[edge[0]][edge[1]]["step"]
            else:
                G_global.add_edge(edge[1],edge[0],strand1 = G[edge[0]][edge[1]]["strand2"].translate(trans_strand),strand2 = G[edge[0]][edge[1]]["strand1"].translate(trans_strand),match = G[edge[0]][edge[1]]["match"],way = G[edge[0]][edge[1]]["way"],step = -1) # G[edge[0]][edge[1]]["step"]

results = []
visited = set()
sources = [G_global[node][list(G_global[node].keys())[0]]["strand1"] + node for (node,val) in G_global.in_degree if(val == 0)]
for node in G_global.nodes():
    G_global.nodes[node]["visit_count"] = max(min(G_global.in_degree[node],G_global.out_degree[node]),1)
print(sources)
for source in sources:
    stack = []
    DFS(source,[],stack)
while([] in results):
    results.remove([])
print(results)

# pos = nx.spring_layout(G_global)
# node_labels = {node: node[6:9] for node in G_global.nodes()}
# edge_labels = nx.get_edge_attributes(G_global,"way")
# nx.draw(G_global,pos,connectionstyle='arc3, rad = 0.2')
# nx.draw_networkx_labels(G_global,pos,labels = node_labels)
# nx.draw_networkx_edge_labels(G_global,pos,edge_labels = edge_labels)
# plt.show()

sequences = {}
i = 1
for result in results:
    sequence = ""
    match = 0
    for j,source in enumerate(result):
        strand = source[0]
        node = source[1 : ]
        if(j != 0):
            match = G_global[result[j - 1][1 : ]][node]["match"]
        if(strand == "+"):
            sequence += GFA.nodes[node]["seq"][match : ]
        else:
            sequence += GFA.nodes[node]["seq"].translate(transtable)[::-1][match : ]
    if(len(result) == 1):
        sequences[result[0][1 : ]] = [result,sequence]
    else:
        sequences["sfd" + "{:0>6d}".format(i) + "l"] = [result,sequence]
        i += 1

length = [GFA.nodes[node]["length"] for node in GFA.nodes()]
avarage_length = sum(length) / len(length)
for node in GFA.nodes(): # retrieve isolated node
    if(GFA.in_degree[node] == 0 and GFA.out_degree[node] == 0): #  and GFA.nodes[node]["length"] > avarage_length
        sequences[node] = [["+" + str(node)],GFA.nodes[node]["seq"]]

new_GFA = nx.DiGraph()

for key,values in sequences.items():
    new_GFA.add_node(key,seq = values[1],units = values[0])

for sequence1,sequence2 in itertools.combinations(sequences.items(),2):
    name1 = sequence1[0]
    name2 = sequence2[0]
    start1 = sequence1[1][0][0][1:]
    start1_strand = sequence1[1][0][0][0]
    end1 = sequence1[1][0][-1][1:]
    end1_strand = sequence1[1][0][-1][0]
    start2 = sequence2[1][0][0][1:]
    start2_strand = sequence2[1][0][0][0]
    end2 = sequence2[1][0][-1][1:]
    end2_strand = sequence2[1][0][-1][0]
    if((end1,start2) in GFA.edges() and GFA[end1][start2]["strand1"] == end1_strand and GFA[end1][start2]["strand2"] == start2_strand):
        new_GFA.add_edge(name1,name2,strand1 = "+",strand2 = "+",match = GFA[end1][start2]["match"])
        new_GFA.add_edge(name2,name1,strand1 = "-",strand2 = "-",match = GFA[end1][start2]["match"])
    elif((end2,start1) in GFA.edges() and GFA[end2][start1]["strand1"] == end2_strand and GFA[end2][start1]["strand2"] == start1_strand):
        new_GFA.add_edge(name2,name1,strand1 = "+",strand2 = "+",match = GFA[end2][start1]["match"])
        new_GFA.add_edge(name1,name2,strand1 = "-",strand2 = "-",match = GFA[end2][start1]["match"])
    elif((start1,start2) in GFA.edges() and GFA[start1][start2]["strand1"] == start1_strand.translate(trans_strand) and GFA[start1][start2]["strand2"] == start2_strand):
        new_GFA.add_edge(name1,name2,strand1 = "-",strand2 = "+",match = GFA[start1][start2]["match"])
        new_GFA.add_edge(name2,name1,strand1 = "-",strand2 = "+",match = GFA[start1][start2]["match"])
    elif((end1,end2) in GFA.edges() and GFA[end1][end2]["strand1"] == end1_strand.translate(trans_strand) and GFA[end1][end2]["strand2"] == end2_strand):
        new_GFA.add_edge(name1,name2,strand1 = "+",strand2 = "-",match = GFA[end1][end2]["match"])
        new_GFA.add_edge(name2,name1,strand1 = "+",strand2 = "-",match = GFA[end1][end2]["match"])

for node in list(new_GFA.nodes()): # resolve small bubble that generate again by loop
    if(node in new_GFA.nodes()):
        neighbors = list(new_GFA[node].keys())
        if(len(neighbors) == 2 and (neighbors[0],neighbors[1]) in new_GFA.edges() and new_GFA[neighbors[0]][node]["strand2"] == new_GFA[node][neighbors[1]]["strand1"] and new_GFA[neighbors[0]][node]["strand1"] == new_GFA[neighbors[0]][neighbors[1]]["strand1"] and new_GFA[node][neighbors[1]]["strand2"] == new_GFA[neighbors[0]][neighbors[1]]["strand2"]):
            sequence = ""
            units = []
            strand1 = new_GFA[neighbors[0]][node]["strand1"]
            strand = new_GFA[neighbors[0]][node]["strand2"]
            strand2 = new_GFA[node][neighbors[1]]["strand2"]
            if(strand1 == "+"):
                sequence += new_GFA.nodes[neighbors[0]]["seq"]
                units += new_GFA.nodes[neighbors[0]]["units"]
            else:
                sequence += new_GFA.nodes[neighbors[0]]["seq"].translate(transtable)[::-1]
                for unit in new_GFA.nodes[neighbors[0]]["units"][::-1]:
                    units.append(unit.translate(trans_strand))
            if(strand == "+"):
                sequence += new_GFA.nodes[node]["seq"][new_GFA[neighbors[0]][node]["match"] : ]
                units += new_GFA.nodes[node]["units"]
            else:
                sequence += new_GFA.nodes[node]["seq"].translate(transtable)[::-1][new_GFA[neighbors[0]][node]["match"] : ]
                for unit in new_GFA.nodes[node]["units"][::-1]:
                    units.append(unit.translate(trans_strand))
            if(strand2 == "+"):
                sequence += new_GFA.nodes[neighbors[1]]["seq"][new_GFA[node][neighbors[1]]["match"] : ]
                units += new_GFA.nodes[neighbors[1]]["units"]
            else:
                sequence += new_GFA.nodes[neighbors[1]]["seq"].translate(transtable)[::-1][new_GFA[node][neighbors[1]]["match"] : ]
                for unit in new_GFA.nodes[neighbors[1]]["units"][::-1]:
                    units.append(unit.translate(trans_strand))

            key = "sfd" + "{:0>6d}".format(i) + "l"
            new_GFA.add_node(key,seq = sequence,units = units)
            i += 1

            for neighbor in new_GFA[neighbors[0]].keys():
                if(neighbor != node and neighbor != neighbors[1]):
                    new_GFA.add_edge(neighbor,key,strand1 = new_GFA[neighbor][neighbors[0]]["strand1"],strand2 = "+",match = new_GFA[neighbor][neighbors[0]]["match"])
                    new_GFA.add_edge(key,neighbor,strand1 = "-",strand2 = new_GFA[neighbor][neighbors[0]]["strand1"].translate(trans_strand),match = new_GFA[neighbor][neighbors[0]]["match"])
            for neighbor in new_GFA[neighbors[1]].keys():
                if(neighbor != node and neighbor != neighbors[0]):
                    new_GFA.add_edge(key,neighbor,strand1 = "+",strand2 = new_GFA[neighbors[1]][neighbor]["strand2"],match = new_GFA[neighbors[1]][neighbor]["match"])
                    new_GFA.add_edge(neighbor,key,strand1 = new_GFA[neighbors[1]][neighbor]["strand2"].translate(trans_strand),strand2 = "-",match = new_GFA[neighbors[1]][neighbor]["match"])
            new_GFA.remove_node(neighbors[0])
            new_GFA.remove_node(node)
            new_GFA.remove_node(neighbors[1])

for tip in list(new_GFA.nodes()): # remove tips
    if(tip in new_GFA.nodes() and len(new_GFA[tip].keys()) == 1):
        for dst in new_GFA[tip].keys():
            if(len(new_GFA.nodes[tip]["seq"]) < 100000 and tip not in homologous and any(node3 != tip and node3 not in new_GFA[tip].keys() and new_GFA[node3][dst]["strand2"] == new_GFA[tip][dst]["strand2"] for node3 in new_GFA[dst].keys())):
                sequence = ""
                units = []
                if(new_GFA[tip][dst]["strand1"] == "+"):
                    sequence += new_GFA.nodes[tip]["seq"]
                    units +=  new_GFA.nodes[tip]["units"]
                else:
                    sequence += new_GFA.nodes[tip]["seq"].translate(transtable)[::-1]
                    for unit in new_GFA.nodes[tip]["units"][::-1]:
                        units.append(unit.translate(trans_strand))
                if(new_GFA[tip][dst]["strand2"] == "+"):
                    sequence += new_GFA.nodes[dst]["seq"][new_GFA[tip][dst]["match"] : ]
                    units +=  new_GFA.nodes[dst]["units"]
                else:
                    sequence += new_GFA.nodes[dst]["seq"].translate(transtable)[::-1][new_GFA[tip][dst]["match"] : ]
                    for unit in new_GFA.nodes[dst]["units"][::-1]:
                        units.append(unit.translate(trans_strand))
                key = "sfd" + "{:0>6d}".format(i) + "l"
                new_GFA.add_node(key,seq = sequence,units = units)
                i += 1

                for neighbor in new_GFA[dst].keys():
                    if(new_GFA[dst][neighbor]["strand1"] == new_GFA[tip][dst]["strand2"]):
                        new_GFA.add_edge(key,neighbor,strand1 = "+",strand2 = new_GFA[dst][neighbor]["strand2"],match = new_GFA[dst][neighbor]["match"])
                        new_GFA.add_edge(neighbor,key,strand1 = new_GFA[dst][neighbor]["strand1"].translate(trans_strand),strand2 = "-",match = new_GFA[dst][neighbor]["match"])
                
                new_GFA.remove_node(tip)
                new_GFA.remove_node(dst)

for (src,dst) in list(new_GFA.edges()): # connect single link
    if((src,dst) in new_GFA.edges()):
        strand1 = new_GFA[src][dst]["strand1"]
        strand2 = new_GFA[src][dst]["strand2"]
        count1,count2 = 0,0
        neighbor1 = [neighbor for neighbor in new_GFA[src].keys() if(new_GFA[src][neighbor]["strand1"] == strand1)]
        neighbor2 = [neighbor for neighbor in new_GFA[dst].keys() if(new_GFA[neighbor][dst]["strand2"] == strand2)]
        if(len(neighbor1) == 1 and len(neighbor2) == 1):
            sequence = ""
            units = []
            if(new_GFA[src][dst]["strand1"] == "+"):
                sequence += new_GFA.nodes[src]["seq"]
                units +=  new_GFA.nodes[src]["units"]
            else:
                sequence += new_GFA.nodes[src]["seq"].translate(transtable)[::-1]
                for unit in new_GFA.nodes[src]["units"][::-1]:
                    units.append(unit.translate(trans_strand))
            if(new_GFA[src][dst]["strand2"] == "+"):
                sequence += new_GFA.nodes[dst]["seq"][new_GFA[src][dst]["match"] : ]
                units +=  new_GFA.nodes[dst]["units"]
            else:
                sequence += new_GFA.nodes[dst]["seq"].translate(transtable)[::-1][new_GFA[src][dst]["match"] : ]
                for unit in new_GFA.nodes[dst]["units"][::-1]:
                    units.append(unit.translate(trans_strand))
            key = "sfd" + "{:0>6d}".format(i) + "l"
            new_GFA.add_node(key,seq = sequence,units = units)
            i += 1

            for neighbor in new_GFA[dst].keys():
                if(new_GFA[dst][neighbor]["strand1"] == new_GFA[src][dst]["strand2"]):
                    new_GFA.add_edge(key,neighbor,strand1 = "+",strand2 = new_GFA[dst][neighbor]["strand2"],match = new_GFA[dst][neighbor]["match"])
                    new_GFA.add_edge(neighbor,key,strand1 = new_GFA[dst][neighbor]["strand2"].translate(trans_strand),strand2 = "-",match = new_GFA[dst][neighbor]["match"])
            for neighbor in new_GFA[src].keys():
                if(new_GFA[neighbor][src]["strand2"] == new_GFA[src][dst]["strand1"]):
                    new_GFA.add_edge(neighbor,key,strand1 = new_GFA[neighbor][src]["strand1"],strand2 = "+",match = new_GFA[neighbor][src]["match"])
                    new_GFA.add_edge(key,neighbor,strand1 = "-",strand2 = new_GFA[neighbor][src]["strand1"].translate(trans_strand),match = new_GFA[neighbor][src]["match"])
            new_GFA.remove_node(src)
            new_GFA.remove_node(dst)

writer = open(output_file,"w")
for node in new_GFA.nodes():
    print(f"S\t{node}\t" + new_GFA.nodes[node]["seq"],file = writer)
    print(node,new_GFA.nodes[node]["units"],len(new_GFA.nodes[node]["seq"]))
for edge in new_GFA.edges():
    print("\t".join(["L",edge[0],new_GFA[edge[0]][edge[1]]["strand1"],edge[1],new_GFA[edge[0]][edge[1]]["strand2"],str(new_GFA[edge[0]][edge[1]]["match"]) + "M"]),file = writer)
writer.close()
