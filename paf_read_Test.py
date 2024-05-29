

def modify_fasta_headers(input_file, output_file):
    count = 1
    with open(input_file, 'r') as input_f, open(output_file, 'w') as output_f:
        for line in input_f:
            if line.startswith('>'):
                header = '>sequence{}'.format(count)
                output_f.write(header + '\n')
                count += 1
            else:
                output_f.write(line)

input_file = '/data/zjjiang/ninanjie/hifiasm/buqi.fasta'
output_file = '/data/zjjiang/ninanjie/hifiasm/output.fasta'
modify_fasta_headers(input_file, output_file)
# with open('/data/zjjiang/HPRC_PAN/hprc-v1.0-minigraph-chm13.gfa', 'r') as f:
#           for line in f:
#             print(len(line.split()))
#             break
        # head_char = line.split()[:7][0]
        # # adding segment_chr information
        # if head_char == 'S':
        #     segmentid = line.split()[:7][1]
        #     self.add_vertex(segmentid)
        #     chrIn = line.split()[:7][4].replace('SN:Z:', '')
        #     self.node_label[segmentid] = {'label': chrIn, 'visited': 0}
        #     if chrIn.startswith("HG") and chrIn not in self.label_to_seg:
        #         self.label_to_seg[chrIn] = segmentid
        # if head_char == 'L':
        #     seg1 = line.split()[:9][1]
        #     seg2 = line.split()[:9][3]
        #     self.add_edge(seg1, seg2, 0)
