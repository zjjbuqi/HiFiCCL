#!/usr/bin/env python3

import argparse
import subprocess
import re
import os
import sys
import multiprocessing
import hificcl_tools
import pysam
import hificcl_ALG
import trans_calls
from collections import Counter


### MAIN PROGRAM ###
def Binningchr(args):
    reffasta, qryfasta, prefix, threads, overwrite, minimapoption1, minimapoption2, max_hang, int_frac, workdir, weight, number= args
    inputdict = hificcl_tools.readFastaAsDict(qryfasta)
    total_len = len(inputdict.keys())

    #start all-vs-all alignment
    #ala_alignment information
    # all_ala_alignment = {}
    # ala_alignment = {}
    print('[Info] Starting all-vs-all alignment in minimap2...')
    alapaf_file = hificcl_tools.minimap(qryfasta, qryfasta, prefix, 'all_vs_all', minimapoption1, overwrite, workdir, 0)
    print('[Info] Starting constructing AlignmentLabelGraph...')
    #set qry-ref alignment information
    ALG = hificcl_ALG.AlignmentLabelGraph(max_hang, int_frac, weight, number)
    # ALG.depth = depth
    ALG.build_ALG(inputdict, alapaf_file)



    #         #only restore the top one in each kind of qryid in paf
    #         if qryid not in all_ala_alignment:
    #             ala_alignment = {'qrylen': qrylen, 'refid': refid, 'reflen': reflen,'match': match,'alignlen': alignlen,'overlap_degree':0,'label':[], 'visit':0}
    #             ala_alignment['overlap_degree'] = min(int(depth)/40,1)*(int(ala_alignment['alignlen'])/int(ala_alignment['qrylen']))+min(int(depth)/40,1)*(int(ala_alignment['alignlen'])/int(ala_alignment['reflen']))\
    #                                           +(int(ala_alignment['match'])/int(ala_alignment['alignlen']))
    #             all_ala_alignment[qryid] = ala_alignment
    #         else:
    #             continue
    #set ref-qry alginment information
    # all_ala_reverse_alignment = {}
    # ala_reverse_alignment = {}
    # with open(alapaf_file, 'r') as f:
    #     for line in f:
    #         qryid, qrylen, qrystart, qryend, strand, refid, reflen, refstart, refend, match, alignlen = line.split()[:11]
    #         #only restore the top one in each kind of qryid in paf
    #         if refid not in all_ala_reverse_alignment:
    #             ala_reverse_alignment = {'qryid': qryid}
    #             all_ala_reverse_alignment[refid] = ala_reverse_alignment
    #         else:
    #             continue
    #
    if ALG.graph == {}:
        print(f'[Error] ALG.graph is NULL.')
        sys.exit(0)
    else:
        print('[Info] The AlignmentLabelGraph is established successfully.')
    # start map_to_reference alignment
    #qryfasta to chr information
    print('[Info] Starting map_to_reference alignment in minimap2...')
    map_file = hificcl_tools.minimap(reffasta, qryfasta, prefix, 'map_to_reference', minimapoption2, overwrite, workdir, 1)
    print('[Info] Starting translocation detection...')
    trans_calling = trans_calls.Trans_calls(map_file)
    trans_calling.run()
    ALG.initial_label(map_file)
    ALG.write_chr_by_chr_reads_intial(workdir, inputdict)
    # ALG.label_filter(trans_calling.removeReadID)
   # ALG.label_set(trans_calling.label_revise_readid)
    #         if rec.query_name not in qry_alignment_chr:
    #             qry_alignment_chr[rec.query_name] = rec.reference_name
    #         else:
    #             continue
    #
    # for key,value in qry_alignment_chr.items():
    #     if key in all_ala_alignment:
    #         all_ala_alignment[key]['label'].append(value)
    # processing the map_to_reference alignment labels' data
    #
    ALG.overlap_degree_sort()
    N = 200
    ALG.label_correct(N)
    # ALG.label_filling(map_file)

    # print('[Info] Starting label correction based on the overlap information between reads...')
    # total_visited = set()
    # corrected = 0
    # for node in ALG.nodelist:
    #     if node in total_visited:
    #         continue
    #     visited, pre_key, pre_overlapdegree = ALG.Findprimary_down(node)
    #     visited_all = ALG.Findprimary_up(pre_key, pre_overlapdegree, visited)
    #     total_visited += visited_all
    #     label_list = []
    #     for node in visited_all:
    #         label_list += node.label
    #     counter = Counter(label_list)
    #     most_common_labels = counter.most_common()
    #     max_label, max_num = most_common_labels[0]
    #     label = max_label
    #
    #     for node in visited_all:
    #          node.label.append(label)
    #
    # for node in ALG.nodelist:
    #     print(f'[Info]{node.node} have {len(node.label)} numbers of label.')
    #     final_label = []
    #     node_counter = Counter(node.label)
    #     print(node_counter)
    #     most_node_labels = node_counter.most_common()
    #     max_node_label, max_node_num = most_node_labels[0]
    #     final_label.append(max_node_label)
    #     if len(node.label)>=100:
    #         for i in range(1, len(most_node_labels)):
    #             alter_node_label, alter_node_num = most_node_labels[i]
    #             if buqi_tools.calculate_digit_difference(max_node_num, alter_node_num):
    #                break
    #             else:
    #                 final_label.append(alter_node_label)
    #     else:
    #         for i in range(1, len(most_node_labels)):
    #             alter_node_label, alter_node_num = most_node_labels[i]
    #             if buqi_tools.calculate_number_difference(max_node_num, alter_node_num):
    #                break
    #             else:
    #                 final_label.append(alter_node_label)
    #     node.label = final_label

    #
    # #using link information
    # label_count = {}
    # corrected_count = 0
    # stack = buqi_tools.Stack()
    # for qryfasta in all_ala_alignment:
    #     print(qryfasta)
    #     label_count = {}
    #     #seek to primary chain.
    #     if all_ala_alignment[qryfasta]['visit']==0:
    #         #Find the precursor node
    #         all_ala_alignment[qryfasta]['visit']=1
    #         if qryfasta in all_ala_reverse_alignment:
    #             current_overlap_degree = all_ala_alignment[qryfasta]['overlap_degree']
    #             pioneerid = all_ala_reverse_alignment[qryfasta]['qryid']
    #             while(all_ala_alignment[pioneerid]['overlap_degree'] >= (1/2) * current_overlap_degree and all_ala_alignment[pioneerid]['visit']==0):
    #                 qryfasta = pioneerid
    #                 all_ala_alignment[qryfasta]['visit']=1
    #                 if qryfasta not in all_ala_reverse_alignment:
    #                     break
    #                 current_overlap_degree = all_ala_alignment[qryfasta]['overlap_degree']
    #                 pioneerid = all_ala_reverse_alignment[qryfasta]['qryid']
    #     else:
    #         continue
    #     if not stack.is_full():
    #         stack.push(qryfasta)
    #     else:
    #         print(f'[Error] Stack is full.')
    #         sys.exit(0)
    #     if all_ala_alignment[qryfasta]['label'][0] not in label_count:
    #         label_count[all_ala_alignment[qryfasta]['label'][0]] = 1
    #     else:
    #         label_count[all_ala_alignment[qryfasta]['label'][0]] += 1
    #     pre_overlap_degree = all_ala_alignment[qryfasta]['overlap_degree']
    #     qryfasta = all_ala_alignment[qryfasta]['refid']
    #     current_overlap_degree = all_ala_alignment[qryfasta]['overlap_degree']
    #     while(current_overlap_degree >= (1/2) * pre_overlap_degree and all_ala_alignment[qryfasta]['visit']==0):
    #         if all_ala_alignment[qryfasta]['label'][0] not in label_count:
    #             label_count[all_ala_alignment[qryfasta]['label'][0]] = 1
    #         else:
    #             label_count[all_ala_alignment[qryfasta]['label'][0]] += 1
    #         if not stack.is_full():
    #             stack.push(qryfasta)
    #         else:
    #             print(f'[Error] Stack is full.')
    #             sys.exit(0)
    #         all_ala_alignment[qryfasta]['visit'] = 1
    #         pre_overlap_degree = all_ala_alignment[qryfasta]['overlap_degree']
    #         qryfasta = all_ala_alignment[qryfasta]['refid']
    #         current_overlap_degree = all_ala_alignment[qryfasta]['overlap_degree']
    #     label = buqi_tools.get_key(label_count,max(label_count.values()))
    #     while(not stack.is_empty()):
    #         qry = stack.peek()
    #         if label not in all_ala_alignment[qry]['label']:
    #             corrected_count += 1
    #             all_ala_alignment[qry]['label'].append (label)
    #         stack.pop()
    # print(f'[Info] There are {corrected} reads which have been corrected.')
    #classify the reads according their 'chr' and write them.
    # chr_classify = {}
    # for node in ALG.nodelist:
    #     for label in node.label:
    #         if label not in chr_classify:
    #             chr_classify[label] = [node.node]
    #         else:
    #             chr_classify[label].append(node.node)
    # #write
    # try:
    #     os.mkdir('chr_by_chr_reads')
    # except:
    #
    # new_dir = os.path.join(workdir, 'chr_by_chr_reads')
    # os.chdir(new_dir)
    # workdir = os.getcwd()
    # file_list = []
    # for label in chr_classify:
    #     filename = workdir +'/'+f"{label}.fasta"
    #     with open(filename, 'w') as chr_bining:
    #         for readid in chr_classify[label]:
    #             seq = inputdict[readid]
    #             chr_bining.write(f'>{readid}\n{seq}\n')
    #     file_list.append(label)
    # print("[Info The reads of each chr have been writted.]")
    file_list = ALG.write_chr_by_chr_reads(workdir, inputdict)
    return file_list

