import collections
import sys

import pysam
from collections import Counter
import hificcl_tools
import os
import random
import copy
import re
#'None' labels corrected
class AlignmentLabelGraph:
    def __init__(self, max_hang, int_frac, weight, number):
        self.node_label = {}
        # self.seg_chr_gfa = {}
        # # self.chr_seg_gfa = {}
        # self.path_gfa = {}
        # self.unkown_chr = {}
        self.max_hang = int(max_hang)
        self.int_frac = float(int_frac)
        self.weight = float(weight)
        self.number = float(number)
        # self.depth = 1
        self.graph = collections.OrderedDict()
        self.graph_reverse = collections.OrderedDict()
        # self.testc = 0

    def contained_data_filter(self):
        pass

    def add_vertex(self, vertex):
        self.graph[vertex] = collections.OrderedDict()
        self.graph_reverse[vertex] = collections.OrderedDict()

    def add_edge(self, source, target, weight):
        if source in self.graph and target in self.graph:
            self.graph[source][target] = weight
            self.graph_reverse[target][source] = weight

    def Findprimary_up(self, start, node_label_temp):
        visited = set()
        end = 0
        pre_key_list = iter(self.graph[start])
        while 1:
            pre_key = next(pre_key_list,None)
            if pre_key != None and node_label_temp[pre_key]['visited'] == 0:
                pre_overlapdegree = self.graph[start][pre_key]
                pre_overlapdegree_initial = self.graph[start][pre_key]
                node_label_temp[pre_key]['visited'] = 1
                visited.add(pre_key)
                break
            elif pre_key == None:

                pre_overlapdegree = 0
                pre_overlapdegree_initial = 0
                pre_key = start
                break
        visited.add(start)
        node_label_temp[start]['visited'] = 1

        while(end!=1):
            # try:
            iter_data = iter(self.graph_reverse[start])
            while(1):
                up_key = next(iter_data, None)
                if (up_key != None and up_key not in visited and node_label_temp[up_key]['visited']==0):
                    if self.weight *pre_overlapdegree  >  self.graph_reverse[start][up_key]:
                        end = 1
                        break
                    visited.add(up_key)
                    node_label_temp[up_key]['visited'] = 1
                    pre_overlapdegree = self.graph_reverse[start][up_key]
                    start = up_key
                    break
                elif up_key == None:
                    end = 1
                    break
            # except StopIteration:
            #     break
        # print(f'[Info] Primary chain {len(visited)}')
        return visited, pre_key, pre_overlapdegree_initial

    def Findprimary_down(self, pre_key, pre_overlapdegree, visited, node_label_temp):

        end = 0
        while(end!=1):
            # try:
            iter_data = iter(self.graph[pre_key])
            while(1):
                next_key = next(iter_data, None)
                if (next_key != None and next_key not in visited and node_label_temp[next_key]['visited']==0):
                    if self.graph[pre_key][next_key] < self.weight * pre_overlapdegree:
                        end = 1
                        break
                    visited.add(next_key)
                    node_label_temp[next_key]['visited'] = 1
                    pre_overlapdegree = self.graph[pre_key][next_key]
                    pre_key = next_key
                    break
                elif next_key == None:
                    end = 1
                    break
            # except StopIteration:
            #     break

        return visited
    # def Findprimary_down(self, start):
    #     visited = set()
    #     end = 0
    #     pre_key_list = iter(self.graph_reverse[start])
    #     # try:
    #     while 1:
    #         pre_key = next(pre_key_list,None)
    #         if pre_key != None and self.node_label[pre_key]['visited'] == 0:
    #             pre_overlapdegree = self.graph[pre_key][start]
    #             pre_overlapdegree_initial = self.graph[pre_key][start]
    #             self.node_label[pre_key]['visited'] = 1
    #             visited.add(pre_key)
    #             break
    #         elif pre_key == None:
    #
    #             pre_overlapdegree = 0
    #             pre_overlapdegree_initial = 0
    #             pre_key = start
    #             break
    #
    #         # pre_overlapdegree = self.graph[pre_key][start]
    #         # pre_overlapdegree_initial = self.graph[pre_key][start]
    #
    #     # except StopIteration:
    #     #     pre_overlapdegree = 0
    #     #     pre_overlapdegree_initial = 0
    #     #     pre_key = start
    #     # if pre_key not in visited:
    #     #     visited.add(pre_key)
    #     if start not in visited:
    #         visited.add(start)
    #     while(end!=1):
    #         # try:
    #         iter_data = iter(self.graph[start])
    #         while(1):
    #             next_key = next(iter_data, None)
    #             if (next_key != None and next_key not in visited and self.node_label[next_key]['visited']==0):
    #                 if self.graph[start][next_key] < 2/3 * pre_overlapdegree:
    #                     end = 1
    #                     break
    #                 visited.add(next_key)
    #                 self.node_label[next_key]['visited'] = 1
    #                 pre_overlapdegree = self.graph[start][next_key]
    #                 start = next_key
    #                 break
    #             elif next_key == None:
    #                 end = 1
    #                 break
    #         # except StopIteration:
    #         #     break
    #
    #     return visited, pre_key, pre_overlapdegree_initial


    # def Findprimary_up(self, next_key, next_overlapdegree, visited):
    #     end = 0
    #     while(end!=1):
    #         # try:
    #         iter_data = iter(self.graph_reverse[next_key])
    #         while(1):
    #             up_key = next(iter_data, None)
    #             if (up_key != None and up_key not in visited and self.node_label[up_key][visited]==0):
    #                 if next_overlapdegree < 2/3 * self.graph_reverse[next_key][up_key]:
    #                     end = 1
    #                     break
    #                 visited.add(up_key)
    #                 self.node_label[up_key][visited] = 1
    #                 next_overlapdegree = self.graph_reverse[next_key][up_key]
    #                 next_key = up_key
    #                 break
    #             elif up_key == None:
    #                 end = 1
    #                 break
    #         # except StopIteration:
    #         #     break
    #     # print(f'[Info] Primary chain {len(visited)}')
    #     return visited


    def add_all_node(self, input_dict):
        for key in input_dict:
            self.add_vertex(key)
            self.node_label[key] = {'label':[], 'visited':0, 'flag':0, 'contained_node_list':[]}

    # def gfa_processed(self, pangfa):
    #     with open(pangfa, 'r') as f:
    #         for line in f:
    #             head_char = line.split()[:7][0]
    #             #adding segment_chr information
    #             if head_char == 'S':
    #                 segmentid= line.split()[:7][1]
    #                 chrIn = line.split()[:7][4].replace('SN:Z:','')
    #                 self.seg_chr_gfa[segmentid] = chrIn
    #             if head_char == 'L':
    #                 seg1 = line.split()[:9][1]
    #                 seg2 =line.split()[:9][3]
    #                 if seg1 not in self.path_gfa:
    #                     self.path_gfa[seg1] = seg2
    #                 else:
    #                     continue
    # def unknown_label_processed(self):
    #     for seg in self.seg_chr_gfa:
    #         if self.seg_chr_gfa[seg].startswith('chr'):
    #             continue
    #         elif self.seg_chr_gfa[seg].startswith('HG'):
    #             print(f'This run {self.seg_chr_gfa[seg]}!')
    #             next_seg = self.path_gfa[seg]
    #             while(self.seg_chr_gfa[next_seg].startswith('HG')):
    #                 next_seg = self.path_gfa[next_seg]
    #             else:
    #                 if self.seg_chr_gfa[seg] not in self.unkown_chr:
    #                     self.unkown_chr[self.seg_chr_gfa[seg]] = self.seg_chr_gfa[next_seg]
    #                     print(f'{self.seg_chr_gfa[seg]} is {self.seg_chr_gfa[next_seg]}')

    # def unkown_label_transform(self, unkown_label):
    #     if unkown_label in self.unkown_chr:
    #         return self.unkown_chr[unkown_label]
    #     else:
    #         print("wrong label_transform!")
    #         sys.exit(1)





    def add_no_contained_edge(self, alapaf_file):
        nodelist = {}
        with open(alapaf_file, 'r') as f:
            for line in f:
                qryid, qrylen, qrystart, qryend, strand, refid, reflen, refstart, refend, match, alignlen = line.split()[:11]

                qrylen = int(qrylen)
                qrystart = int(qrystart)
                qryend = int(qryend)

                reflen = int(reflen)
                refstart = int(refstart)
                refend = int(refend)
                match = int(match)
                alignlen = int(alignlen)
                if qryid not in nodelist:

                    l5 = reflen - refend if strand else refstart
                    l3 = refstart if strand else reflen - refend
                    if qrylen >> 1 > reflen:
                        if l5 > self.max_hang >> 2 or l3 > self.max_hang >> 2 or refend - refstart < reflen * self.int_frac:
                            self.node_label[qryid]['contained_node_list'].append(refid)
                            continue

                    elif qrylen < reflen >> 1:
                        if qrystart > self.max_hang >> 2 or qrylen - qryend > self.max_hang >> 2 or qryend - qrystart < qrylen * self.int_frac:
                            self.node_label[refid]['contained_node_list'].append(qryid)
                            continue

                    # weight = min(int(self.depth) / 40, 1.2) * (int(alignlen) / int(qrylen)) + min(int(self.depth) / 40, 1.2) * (
                    #             int(alignlen) / int(reflen)) + (int(match) / int(alignlen))
                    weight = (int(alignlen) / int(qrylen))*(int(alignlen) / int(reflen)) * (int(match) / int(alignlen)) * 1000
                    self.add_edge(qryid, refid, weight)
                # elif nodelist[qryid] == 1:
                #     nodelist[qryid] += 1
                #     weight = (int(alignlen) / int(qrylen)) * (int(alignlen) / int(reflen)) * (
                #                 int(match) / int(alignlen)) * 1000
                #     self.add_edge(qryid, refid, weight)
                else:
                    if len(self.graph[qryid]) == 0:
                        l5 = reflen - refend if strand else refstart
                        l3 = refstart if strand else reflen - refend
                        if qrylen >> 1 > reflen:
                            if l5 > self.max_hang >> 2 or l3 > self.max_hang >> 2 or refend - refstart < reflen * self.int_frac:
                                self.node_label[qryid]['contained_node_list'].append(refid)
                                continue

                        elif qrylen < reflen >> 1:
                            if qrystart > self.max_hang >> 2 or qrylen - qryend > self.max_hang >> 2 or qryend - qrystart < qrylen * self.int_frac:
                                self.node_label[refid]['contained_node_list'].append(qryid)
                                continue

                        # weight = min(int(self.depth) / 40, 1.2) * (int(alignlen) / int(qrylen)) + min(int(self.depth) / 40, 1.2) * (
                        #             int(alignlen) / int(reflen)) + (int(match) / int(alignlen))
                        weight = (int(alignlen) / int(qrylen) )* (int(alignlen) / int(reflen))*(int(match) / int(alignlen))* 1000
                        self.add_edge(qryid, refid, weight)
                # else:
                #     continue
        print(f'[Info] All edges were processed.')

    def build_ALG(self, input_dict, alapaf_file):
        self.add_all_node(input_dict)
        self.add_no_contained_edge(alapaf_file)

    def overlap_degree_sort(self):
        # for node in self.graph:
        #     sorted_dict = collections.OrderedDict(sorted(self.graph[node].items(), key=lambda item: item[1],reverse=True))
        #     self.graph[node] = sorted_dict
        for node in self.graph_reverse:
            sorted_dict_reverse = collections.OrderedDict(sorted(self.graph_reverse[node].items(), key=lambda item: item[1],reverse=True))
            self.graph_reverse[node] = sorted_dict_reverse

    def initial_label(self, map_file):
        with pysam.AlignmentFile(map_file, 'r') as sam_rd:
            for rec in sam_rd:
                if rec.query_name in self.node_label and rec.reference_name not in self.node_label[rec.query_name]['label']:
                    if rec.reference_name != None:
                        self.node_label[rec.query_name]['label'].append(rec.reference_name)
                    else:
                        self.node_label[rec.query_name]['label'].append("None")
        for node in self.node_label:
            if len(self.node_label[node]['label']) == 0 :
                self.node_label[node]['label'].append("None")
                print(f"{node}'s chr is None!")
        print('[Info] Each reads has been labeled.')

    def initial_label_pan(self, map_file):
        lastid = 0
        count1 = 0
        count2 = 0
        count3 = 0
        # self.unknown_label_processed()
        with open(map_file, 'r') as f:
            for line in f:
                qryid, qrylen, qrystart, qryend, strand, refid, reflen, refstart, refend, match, alignlen = line.split()[:11]
                # if qryid == lastid:
                #     continue
                # else:
                #     lastid = qryid

                chr_id = self.extract_chr(refid)

                if qryid in self.node_label and chr_id not in self.node_label[qryid]['label']:
                    self.node_label[qryid]['label'].append(chr_id)
                    count1 += 1
                    # print(f"{qryid}'s chr is {chr_id}")


                # print(f"{node}'s chr is None!")
        print(f'[Info] There are {count1} sequences assigned labels.')
        # print(f'[Info] There are {self.testc} special!')

    def label_filter(self, removeReadID):
        for id in removeReadID:
            self.node_label[id]['label'] = []

    def label_set(self, label_revise_readid):
        for ind in label_revise_readid:
            for id in label_revise_readid[ind]['readid']:
                if self.node_label[id]['flag'] == 1:
                    continue
                CHR = 'chr' + str(label_revise_readid[ind]['chr'])
                self.node_label[id]['label'].append(CHR)
                self.node_label[id]['flag'] = 1

    def label_correct(self, N):
        global label
        print('[Info] Starting label correction based on the overlap information between reads...')
        # intial_label_dict = self.node_label
        # total_visited = set()
        node_label_total = {}
        corrected_final = 0
        # corrected = 0
        # count = 0
        for i in range(N):
            # for node in self.node_label:
            #     self.node_label[node]['visited'] = 0
            #     node_label_temp[node] = self.node_label[node]['label']
            corrected = 0
            node_label_temp = copy.deepcopy(self.node_label)
            nodes = list(self.node_label.keys())
            random.shuffle(nodes)
            for node in nodes:
                # count += 1
                # if count % 10000 == 0:
                #     print(f'[Info] There were {count} nodes having been processed... ')
                if node_label_temp[node]['visited']==1:
                    continue
                visited, pre_key, pre_overlapdegree = self.Findprimary_up(node, node_label_temp)
                visited_all = self.Findprimary_down(pre_key, pre_overlapdegree, visited, node_label_temp)
                # if len(visited_all)==1:
                #     node = next(iter(visited_all))
                #     new_down_node = next(iter(self.graph[node]), None)
                #     new_up_node = next(iter(self.graph_reverse[node]), None)
                #     if new_up_node != None:
                #         visited_all.add(new_up_node)
                #         node_label_temp[node] = node_label_temp[new_up_node]
                #     elif new_down_node != None:
                #         visited_all.add(new_down_node)
                #         node_label_temp[node] = node_label_temp[new_down_node]
                visited_all_length = len(visited_all)
                for node in visited_all:
                    if len(node_label_temp[node]['contained_node_list']) != 0:
                        visited_all_length = visited_all_length + len(node_label_temp[node]['contained_node_list'])

                # total_visited = total_visited | visited_all
                if visited_all_length < self.number:
                    continue


                    # print(f'[Info] Finding one primary chaining containing {visited_all_length} nodes!')
                label_list = []

                label_unified_flag = 0
                for node in visited_all:
                    if node_label_temp[node]['flag'] == 1:
                        label_unified_flag = 1
                        label = node_label_temp[node]['label'][0]
                        break
                    elif node_label_temp[node]['label'][0] != 'None':
                    # elif len(node_label_temp[node]['label'])!=0 and node_label_temp[node]['label'][0] != 'None':
                        label_list += node_label_temp[node]['label']
                # try:
                if label_unified_flag == 0:
                    if len(label_list) != 0:
                        counter = Counter(label_list)
                    # except IndexError:
                    #     print(f'[Info] {node} is {self.node_label[node]['label']} IndexError')
                    #     continue
                        most_common_labels = counter.most_common()
                        max_label, max_num = most_common_labels[0]
                        label = max_label

                for node in visited_all:
                    if label not in node_label_temp[node]['label']:
                        corrected += 1
                        node_label_temp[node]['label'] = [label]
                    else:
                        node_label_temp[node]['label'] = [label]
                    if len(node_label_temp[node]['contained_node_list'])!=0:
                        for node in node_label_temp[node]['contained_node_list']:
                            if label not in node_label_temp[node]['label']:
                                corrected += 1
                                node_label_temp[node]['label'] = [label]
                            else:
                                node_label_temp[node]['label'] = [label]

            for node in node_label_temp:
                # if len(node_label_temp[node]['label']) != 0 and node_label_temp[node]['label'] != None:
                if node_label_temp[node]['label'][0] != "None":
                    if node not in node_label_total:
                        node_label_total[node] = [node_label_temp[node]['label'][0]]
                    else:
                        node_label_total[node].append(node_label_temp[node]['label'][0])
            # print(f"[Info] In the {i}th round, the {corrected} reads have been corrected.")

        for node in node_label_total:
            # for item in node_label_total[node]:
            #     if isinstance(item, list):
            #         print(f"bug is {node_label_total[node]}")
            counter = Counter(node_label_total[node])
            most_common_labels = counter.most_common()
            max_label, max_num = most_common_labels[0]
            label = max_label

            if self.node_label[node]['label'][0]!=label:
                corrected_final += 1
                self.node_label[node]['label'] = [label]
        print(f"[Info] The {corrected_final} reads were corrected!")
            # print(f"[Info] The {corrected} reads have been corrected.")
        # for node in self.node_label:
        #
        #     final_label = []
        #     try:
        #         for label in intial_label_dict[node]:
        #             self.node_label[node].remove(label)
        #     except ValueError:
        #         print('ValueError')
        #     if len(self.node_label[node])==0:
        #         self.node_label[node]=intial_label_dict[node]
        #     if len(self.node_label[node])==0:
        #         print(f'[Info]{node} have {len(self.node_label[node])} numbers of label.biaoqian 0')
        #         continue
        #     # print(f'[Info]{node} have {len(self.node_label[node])} numbers of label.')

        #     most_node_labels = node_counter.most_common()
        #     max_node_label, max_node_num = most_node_labels[0]
        #     final_label.append(max_node_label)
            # if len(self.node_label[node]) >= 100:
            #     for i in range(1, len(most_node_labels)):
            #         alter_node_label, alter_node_num = most_node_labels[i]
            #         if buqi_tools.calculate_digit_difference(max_node_num, alter_node_num):
            #             break
            #         else:
            #             final_label.append(alter_node_label)
            # elif len(self.node_label[node]) >= 10:
            #     for i in range(1, len(most_node_labels)):
            #         alter_node_label, alter_node_num = most_node_labels[i]
            #         if buqi_tools.calculate_number_difference1(max_node_num, alter_node_num):
            #             break
            #         else:
            #             final_label.append(alter_node_label)
            # else:
            #     for i in range(1, len(most_node_labels)):
            #         alter_node_label, alter_node_num = most_node_labels[i]
            #         if buqi_tools.calculate_number_difference2(max_node_num, alter_node_num):
            #             break
            #         else:
            #             final_label.append(alter_node_label)

            # self.node_label[node] = final_label

    def write_chr_by_chr_reads(self, workdir, inputdict):
        chr_classify = {}
        for node in self.node_label:
            if len(self.node_label[node]['label']) !=0:
                if self.node_label[node]['label'][0] not in chr_classify:
                    chr_classify[self.node_label[node]['label'][0]] = [node]
                else:
                    chr_classify[self.node_label[node]['label'][0]].append(node)
            # for label in self.node_label[node]:
            #     if label not in chr_classify:
            #         chr_classify[label] = [node]
                # else:
                #     chr_classify[label].append(node)
        # write
        os.chdir(workdir)
        try:
            os.mkdir('chr_by_chr_reads')
        except:
            pass
        new_dir = os.path.join(workdir, 'chr_by_chr_reads')
        os.chdir(new_dir)
        workdir = os.getcwd()
        file_list = []
        for label in chr_classify:
            filename = workdir + '/' + f"{label}.fasta"
            with open(filename, 'w') as chr_bining:
                for readid in chr_classify[label]:
                    seq = inputdict[readid]
                    chr_bining.write(f'>{readid}\n{seq}\n')
            file_list.append(label)
        print("[Info] The reads of each chr have been writted.")
        return file_list

    def write_chr_by_chr_reads_intial(self, workdir, inputdict):
        chr_classify = {}
        for node in self.node_label:
            if len(self.node_label[node]['label']) != 0:
                if self.node_label[node]['label'][0] not in chr_classify:
                    chr_classify[self.node_label[node]['label'][0]] = [node]
                else:
                    chr_classify[self.node_label[node]['label'][0]].append(node)
            # for label in self.node_label[node]:
            #     if label not in chr_classify:
            #         chr_classify[label] = [node]
            # else:
            #     chr_classify[label].append(node)
        # write
        os.chdir(workdir)
        try:
            os.mkdir('chr_by_chr_reads_initial')
        except:
            pass
        new_dir = os.path.join(workdir, 'chr_by_chr_reads_initial')
        os.chdir(new_dir)
        workdir = os.getcwd()
        file_list = []
        for label in chr_classify:
            filename = workdir + '/' + f"{label}.fasta"
            with open(filename, 'w') as chr_bining:
                for readid in chr_classify[label]:
                    seq = inputdict[readid]
                    chr_bining.write(f'>{readid}\n{seq}\n')
            file_list.append(label)
        print("[Info] The reads of each chr have been writted.")
        return workdir + '/' + "None.fasta"

    def label_filling(self, map_file):
        with pysam.AlignmentFile(map_file, 'r') as sam_rd:
            for rec in sam_rd:
                if rec.query_name in self.node_label and len(self.node_label[rec.query_name]['label'])==0:
                    self.node_label[rec.query_name]['label'].append(rec.reference_name)
        print("[Info] Each read's label has been filling.")

    # def extract_chr(self, string):
    #     if string.startswith('chr'):
    #         return string
    #     match =re.search(r"(chr[^:]*):", string)
    #     if match:
    #         return match.group(1)
    #     else:
    #         match = re.search(r"(HG[^:]*):", string)
    #         # self.testc += 1
    #         if match:
    #             # print(f'Processing [Info] {match.group(1)}')
    #             # label = self.unkown_label_transform(match.group(1))
    #             # print(f'[Info] Transform is {label}')
    #             return match.group(1)
    #         else:
    #             match = re.search(r"(NA[^:]*):", string)
    #             if match:
    #                 return match.group(1)
    #
    #             else:
    #                 print(f"wrong chr information!!!{string}")
    #                 sys.exit(1)







