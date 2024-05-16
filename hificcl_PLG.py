import collections
import re


class PangenomeLabelGraph:
    def __init__(self):
        self.node_label = {}
        self.graph = collections.OrderedDict()
        self.graph_reverse = collections.OrderedDict()
        self.testc = 0
        self.label_to_chr = {}
        self.label_to_seg = {}
        self.chrid = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                      'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrM','chrX','chrY']

    def add_vertex(self, vertex):
        self.graph[vertex] = collections.OrderedDict()
        self.graph_reverse[vertex] = collections.OrderedDict()

    def add_edge(self, source, target, weight):
        if source in self.graph and target in self.graph:
            self.graph[source][target] = weight
            self.graph_reverse[target][source] = weight


    # def add_all_node(self, input_dict):
    #     for key in input_dict:
    #         self.add_vertex(key)
    #         self.node_label[key] = {'label':[], 'visited':0}

    def gfa_processed(self, pangfa):
        with open(pangfa, 'r') as f:
            for line in f:
                # if len(line.split())<7:
                #     continue
                head_char = line.split()[:7][0]
                # adding segment_chr information
                if head_char == 'S':
                    segmentid = line.split()[:7][1]
                    self.add_vertex(segmentid)
                    chrIn = str(line.split()[:7][4].replace('SN:Z:', ''))
                    if chrIn.startswith("GR"):
                        pattern = r'chr.*'
                        match = re.search(pattern,chrIn)
                        if match:
                            result = match.group()
                            chrIn = result
                        # else:
                        #     chrIn = "None"
                    self.node_label[segmentid] = {'label': chrIn, 'visited': 0}
                    if chrIn.startswith("HG") and chrIn not in self.label_to_seg:
                        self.label_to_seg[chrIn] = segmentid
                    if chrIn.startswith("NA") and chrIn not in self.label_to_seg:
                        self.label_to_seg[chrIn] = segmentid
                    if chrIn.startswith('chr') and chrIn not in self.label_to_seg and chrIn not in self.chrid:
                        self.label_to_seg[chrIn] =segmentid
                if head_char == 'L':
                    seg1 = line.split()[:9][1]
                    seg2 = line.split()[:9][3]
                    self.add_edge(seg1,seg2,0)

    def Findprimary_up(self, string):
        start = self.label_to_seg[string]
        end = 0
        visited = set()
        while (end != 1):
            # try:
            iter_data = iter(self.graph_reverse[start])
            while (1):
                up_key = next(iter_data, None)
                visited.add(up_key)
                if (up_key != None and  self.node_label[up_key]['visited'] == 0):
                    if self.node_label[up_key]['label'].startswith('chr') and self.node_label[up_key]['label'] in self.chrid:
                        for node in visited:
                            self.node_label[node]['visited']=0
                        if string not in self.label_to_chr:
                            self.label_to_chr[string] = self.node_label[up_key]['label']
                        return self.node_label[up_key]['label']
                    self.node_label[up_key]['visited'] = 1
                    start = up_key
                    break
                elif up_key == None:
                    end = 1
                    break
        return visited

    def Findprimary_down(self, string, visited):
        start = self.label_to_seg[string]
        end = 0
        while (end != 1):
            # try:
            iter_data = iter(self.graph[start])
            while (1):
                next_key = next(iter_data, None)
                if (next_key != None and self.node_label[next_key]['visited'] == 0):
                    if self.node_label[next_key]['label'].startswith('chr') and self.node_label[next_key]['label'] in self.chrid:
                        for node in visited:
                            self.node_label[node]['visited'] = 0
                        if string not in self.label_to_chr:
                            self.label_to_chr[string] = self.node_label[next_key]['label']
                        return self.node_label[next_key]['label']
                    self.node_label[next_key]['visited'] = 1
                    start = next_key
                    break
                elif next_key == None:
                    end = 1
                    break

        for node in visited:
            self.node_label[node]['visited'] = 0
        self.testc += 1
        return "None"


    def unknown_label_transform(self, node_label):
        for node in node_label:
            if node_label[node]['label'][0] == "None":
                continue
            if node_label[node]['label'][0].startswith('chr') and node_label[node]['label'][0] in self.chrid:
                continue
            elif node_label[node]['label'][0].startswith('HG') or node_label[node]['label'][0].startswith('NA') or node_label[node]['label'][0].startswith('chr'):
                # print(f"This run {node_label[node]['label'][0]}!")
                if node_label[node]['label'][0] in self.label_to_chr:
                    node_label[node]['label'][0] =  self.label_to_chr[node_label[node]['label'][0]]
                    # print(f"Transform is {node_label[node]['label'][0]}")
                else:
                    flag = self.Findprimary_up(node_label[node]['label'][0])
                    if flag.startswith('chr') and flag in self.chrid:
                        if node_label[node]['label'][0] not in self.label_to_chr:
                            self.label_to_chr[node_label[node]['label'][0]] = flag
                        node_label[node]['label'][0] = flag
                        # print(f"Transform is {node_label[node]['label'][0]}")
                    else:
                        flag2 = self.Findprimary_down(node_label[node]['label'][0], flag)
                        if flag2.startswith('chr') and flag2 in self.chrid:
                            if node_label[node]['label'][0] not in self.label_to_chr:
                                self.label_to_chr[node_label[node]['label'][0]] = flag2
                            node_label[node]['label'][0] = flag2
                            # print(f"Transform is {node_label[node]['label'][0]}")
                        else:
                            node_label[node]['label'][0] = 'None'
                            # print(f"Transform is {node_label[node]['label'][0]}")
        print(f"There are {self.testc} None!")

    def build_PLG(self, pangfa):
        self.gfa_processed(pangfa)








