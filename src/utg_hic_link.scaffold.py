#!/usr/bin/env python
from collections import defaultdict
import copy
import sys
import networkx as nx

trans_strand = str.maketrans('+-',"-+")
fa_file = sys.argv[1]

def getUtgDict(filePath):
    myDict = {}
    preUtg = None
    with open(filePath, mode='r') as f:
        for line in f:
            splits = line.strip().split()
            if len(splits) == 1:
                preUtg = splits[0][1:]
            # elif(splits[2] == "+"):
            else:
                myDict[splits[1]] = [preUtg,splits[0]]
    return myDict


def getUtg(readDict): # amount controll,non-overlapping kmer
    maxutg =  max(readDict.keys(),key = lambda key : len(readDict[key]))
    
    """if(len(readDict[maxutg]) > 1 and len(list(zip(*readDict[maxutg]))[0]) == len(set(list(zip(*readDict[maxutg]))[0])) and int(readDict[maxutg][-1][2],16) - int(readDict[maxutg][0][2],16) >= 31):
        return maxutg"""
    if(len(readDict[maxutg]) / sum(map(len,readDict.values())) >= 0.5 and len(list(zip(*readDict[maxutg]))[0]) == len(set(list(zip(*readDict[maxutg]))[0]))):
        return maxutg
    else:
        return None


def deelPreItem(preRead, readDict, fileType, matchDict, resultDict_same, resultDict_diff):
    utg = getUtg(readDict)
    if(utg != None):
        buck = defaultdict(list)
        for (cordinate,strand,_) in readDict[utg]:
            buck[int(cordinate,16) // 16].append([int(cordinate,16),strand])
        tempbuck = copy.copy(buck)
        max_key = max(tempbuck.keys(),key = lambda key : len(buck[key] + buck[key + 1] + buck[key - 1]))
        temp = buck[max_key - 1] + buck[max_key] + buck[max_key + 1]
        temp.sort(key = lambda item : item[0])
        cordinate = temp[len(temp) // 2][0] # middle
        strand = temp[len(temp) // 2][1]
        # part = 1 if(cordinate <= length[utg] / 2) else 2
        # print(preRead,utg,readDict,cordinate)
        if fileType == 'read1':
            matchDict[preRead] = [utg,cordinate,strand]
        else:
            r1 = preRead.replace('/2', '/1')
            if r1 in matchDict:  # 表示read1已经匹配到了
                lstTmp = [matchDict[r1][0], utg]
                strand1 = matchDict[r1][2]
                strand2 = strand
                writer2.write(preRead + '\t' + utg + '\t' + str(cordinate) + '\n')
                if(lstTmp[0] != lstTmp[1]):
                    lstTmp.sort()
                    # resultDict[lstTmp[0] + "\t" + lstTmp[1]] += 1
                    if(strand1 == strand2):
                        resultDict_same[lstTmp[0] + "\t" + lstTmp[1]] += 1
                    else:
                        resultDict_diff[lstTmp[0] + "\t" + lstTmp[1]] += 1
                """cordinate1 = matchDict[r1][1]
                strand1 = matchDict[r1][3]
                cordinates1 = matchDict[r1][4]
                cordinate2 = cordinate
                strand2 = strand
                cordinates2 = readDict[utg]
                if(lstTmp[0] != lstTmp[1]):
                    print(preRead,lstTmp[0][1:],cordinates1,cordinate1,strand1,lstTmp[1][1:],cordinates2,cordinate2,strand2)
                    if((lstTmp[0][1:],lstTmp[1][1:]) in GFA.edges()):
                        match = GFA[lstTmp[0][1:]][lstTmp[1][1:]]["match"]
                    else:
                        match = 0
                    if(strand1 == "-"):
                        cordinate1 = GFA.nodes[lstTmp[0][1:]]["length"] - cordinate1 + 1
                    if(strand2 == "-"):
                        cordinate2 = GFA.nodes[lstTmp[1][1:]]["length"] - cordinate2 + 1
                    distance = GFA.nodes[lstTmp[0][1:]]["length"] - cordinate1 + cordinate2 - 1 - match
                    if((side1 != side2 and strand1 == strand2) or (side1 == side2 and strand1 != strand2)):
                        print(1)
                        if(strand2.translate(trans_strand) + lstTmp[1][1:] + strand1.translate(trans_strand) + lstTmp[0][1:] in resultDict):
                            resultDict[strand2.translate(trans_strand) + lstTmp[1][1:] + strand1.translate(trans_strand) + lstTmp[0][1:]].append(distance)
                        else:
                            resultDict[strand1 + lstTmp[0][1:] + strand2 + lstTmp[1][1:]].append(distance)"""
                    

def deelReadFile(filePath, matchDict, fileType, resultDict_same=None, resultDict_diff=None):
    readDict = defaultdict(list)
    preRead = None
    with open(filePath, mode='r') as f:
        for line in f:
            splits = line.split()
            if len(splits) == 1:
                if preRead != None and len(readDict) != 0:
                    deelPreItem(preRead, readDict, fileType, matchDict, resultDict_same, resultDict_diff)
                preRead = splits[0][1:]
                readDict.clear()
            else:
                if splits[1] in myDict:
                    theUtg = myDict[splits[1]][0]
                    readDict[theUtg].append((myDict[splits[1]][1],splits[2],splits[0]))
    if len(readDict) != 0:
        deelPreItem(preRead, readDict, fileType, matchDict, resultDict_same, resultDict_diff)
    readDict.clear()


if __name__ == '__main__':
    length = defaultdict(int)
    name = ""
    with open(fa_file) as file:
        for line in file:
            line = line.strip()
            if(len(line) > 0 and line[0] == ">"):
                name = line[1:]
            else:
                length[name] = len(line)
    utgPath, reads1Path, reads2Path, reads1LogPath, reads2LogPath, resultLogPath = 'ref.pos/ref.pos', 'read1.pos/ref.pos', 'read2.pos/ref.pos', 'read1.log', 'read2.log', 'result.log'
    myDict, resultDict_same, resultDict_diff, matchDict = getUtgDict(utgPath), defaultdict(int), defaultdict(int), {}
    writer1 = open(reads1LogPath,"w")
    writer2 = open(reads2LogPath,"w")
    deelReadFile(reads1Path, matchDict, 'read1', resultDict_same, resultDict_diff)
    for k in matchDict:
        writer1.write('{0}\t{1}\t{2}\n'.format(k, matchDict[k][0],matchDict[k][1]))
    deelReadFile(reads2Path, matchDict, 'read2', resultDict_same, resultDict_diff)
    writer1.close()
    writer2.close()
    writer = open(resultLogPath,"w")
    for utgMapping, value in sorted(resultDict_same.items(), key=lambda item: item[1],reverse = True):
        writer.write('{0}\tsame\t{1}\n'.format(utgMapping, value))
    for utgMapping, value in sorted(resultDict_diff.items(), key=lambda item: item[1],reverse = True):
        writer.write('{0}\tdiff\t{1}\n'.format(utgMapping, value))
    # writer.write("mat1_cen1_mat2/pat1_cen2_pat2:\t" + str(resultDict["cen000001l\tmat000001l"] + resultDict["cen000001l\tmat000002l"] + resultDict["cen000002l\tpat000001l"] + resultDict["cen000002l\tpat000002l"]) + "\n")
    # writer.write("mat1_cen2_mat2/pat1_cen1_pat2:\t" + str(resultDict["cen000002l\tmat000001l"] + resultDict["cen000002l\tmat000002l"] + resultDict["cen000001l\tpat000001l"] + resultDict["cen000001l\tpat000002l"]) + "\n")
    # temp = defaultdict(float)
    # for utgMapping, value in sorted(resultDict.items(), key=lambda item: item[1] / (length[item[0].split()[0]] * length[item[0].split()[1]]),reverse = True):
    #    writer.write('{0}\t{1}\n'.format(utgMapping, value / (length[utgMapping.split()[0]] * length[utgMapping.split()[1]])))
    #    temp[utgMapping] = value / (length[utgMapping.split()[0]] * length[utgMapping.split()[1]])
    # writer.write("mat1_cen1_mat2/pat1_cen2_pat2:\t" + str(temp["cen000001l\tmat000001l"] + temp["cen000001l\tmat000002l"] + temp["cen000002l\tpat000001l"] + temp["cen000002l\tpat000002l"]) + "\n")
    # writer.write("mat1_cen2_mat2/pat1_cen1_pat2:\t" + str(temp["cen000002l\tmat000001l"] + temp["cen000002l\tmat000002l"] + temp["cen000001l\tpat000001l"] + temp["cen000001l\tpat000002l"]) + "\n")
    writer.close()
    