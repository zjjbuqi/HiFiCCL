import pysam
import os
import sys
import re
from collections import Counter

class Trans_calls:
    def __init__(self, file_name, mapq = 20, alignlen = 1000):
        self.filename = file_name
        self.mapq = mapq
        self.alignlen = alignlen
        self.readInformation = {}
        self.removeReadID = []
        self.tran_signal = []
        self.cluster = {}
        self.breakpoint_all = {}
        self.label_revise_readid = {}
        self.ind = 0

        if not os.path.exists(self.filename):
            print(f"File {self.filename} does not exist!")
            sys.exit(0)
        if os.stat(self.filename).st_size == 0:
            print(f"{self.filename} is empty!")
            sys.exit(0)

    def sparse_flag(self,flag):
        hex_num = hex(flag)
        last_digit = hex_num[-1]
        if last_digit == '0':
            return "+"
        else:
            return "-"

    def c_pos(self, cigar, refstart):
        number = ''
        numlist = [str(i) for i in range(10)]
        readstart = False
        readend = False
        refend = False
        readloc = 0
        refloc = refstart
        for c in cigar:
            if (c in numlist):
                number += c
            else:
                number = int(number)
                if (readstart == False and c in ['M', 'I', '=', 'X']):
                    readstart = readloc
                if (readstart != False and c in ['H', 'S']):
                    readend = readloc
                    refend = refloc
                    break

                if (c in ['M', 'I', 'S', '=', 'X']):
                    readloc += number

                if (c in ['M', 'D', 'N', '=', 'X']):
                    refloc += number
                number = ''
        if (readend == False):
            readend = readloc
            refend = refloc

        return refstart, refend, readstart, readend

    def read_from_mapfile(self):
        samfile = pysam.AlignmentFile(self.filename, "r")
        return samfile

    def data_process(self, samfile):
        for read in samfile:
            if read.mapping_quality <= 20:
                self.removeReadID.append(read.query_name)
                continue
            if read.reference_length<=1000:
                self.removeReadID.append(read.query_name)
                continue
            if read.query_name in self.readInformation:
                continue
            if read.has_tag('SA')==True:
                rawsplit = read.get_tag('SA').split(';')
                if len(rawsplit)>7:
                    self.removeReadID.append(read.query_name)
                    continue
            refstart, refend, readstart, readend = self.c_pos(read.cigarstring, read.reference_start)
            strand_read = self.sparse_flag(read.flag)
            self.readInformation[read.query_name]={"chr": read.reference_name, "qrys": readstart,"qrye": readend, "refs": refstart , "refe": refend, "strand": strand_read, "sa": {}}
            if(read.has_tag('SA')==True):
                rawsplit = read.get_tag('SA').split(';')
                for sa in rawsplit:
                    if sa != '':
                        sainfo = sa.split(',')
                        sachr, sarefs, strand, cigar = sainfo[0], int(sainfo[1]), sainfo[2], sainfo[3]
                        if sachr == 'chrM':
                            continue
                        sa_refstart, sa_refend, sa_readstart, sa_readend = self.c_pos(cigar, sarefs)
                        if abs(sa_refstart - refstart) <= 100 and read.reference_name != sachr:
                            self.readInformation[read.query_name]['sa']= {"chr": sachr, "qrys": sa_readstart, "qrye": sa_readend, "refs": sa_refstart, "refe": sa_refend, "strand": strand}
                            break


    def tran_sigal_exact(self):
        pattern = r'\d+'
        for key in self.readInformation:
            if len(self.readInformation[key]['sa'])!=0:
                readid = key
                Chr1 = self.readInformation[key]['chr']
                number1 =int(re.search(pattern, Chr1).group())
                Chr2 = self.readInformation[key]['sa']['chr']
                number2 = int(re.search(pattern, Chr2).group())
                strand1 = self.readInformation[key]['strand']
                strand2 = self.readInformation[key]['sa']['strand']
                ref1s = self.readInformation[key]['refs']
                ref1e = self.readInformation[key]['refe']
                ref2s = self.readInformation[key]['sa']['refs']
                ref2e = self.readInformation[key]['sa']['refe']
                self.tran_sigal_process(number1, number2, strand1, strand2, ref1s, ref1e, ref2s, ref2e, readid)


    def tran_sigal_process(self, Chr1, Chr2, strand1, strand2, ref1s, ref1e, ref2s, ref2e, readid):
        if Chr1 < Chr2 and strand1 == "+" and strand2 == "+":
            self.tran_signal.append((Chr1, ref1e, Chr2, ref2e, readid))
        if Chr2 < Chr1 and strand1 == '+' and strand2 == '+':
            self.tran_signal.append((Chr2, ref2s, Chr1, ref1e, readid))
        if Chr1 < Chr2 and strand1 == '+' and strand2 == '-':
            self.tran_signal.append((Chr1, ref1e, Chr2, ref2e, readid))
        if Chr2 < Chr1 and strand1 == '+' and strand2 == '-':
            self.tran_signal.append((Chr2, ref2e, Chr1, ref1e, readid))
        if Chr1 < Chr2 and strand1 == '-' and strand2 == '+':
            self.tran_signal.append((Chr1, ref1s, Chr2, ref2s, readid))
        if Chr2 < Chr1 and strand1 == '-' and strand2 == '+':
            self.tran_signal.append((Chr2, ref2s, Chr1, ref1s, readid))
        if Chr1 < Chr2 and strand1 == '-' and strand2 == '-':
            self.tran_signal.append((Chr1, ref1s, Chr2, ref2e, readid))
        if Chr2 < Chr1 and strand1 == '-' and strand2 == '-':
            self.tran_signal.append((Chr2, ref2e, Chr1, ref1s, readid))

    def clustering(self):
        latest_signal = {}
        self.tran_signal.sort(key=self.sort_key)
        cluster_i=0

        for signal in self.tran_signal:
            flag = 0
            cluster_number = len(self.cluster)
            if cluster_number == 0:
                self.cluster[cluster_i] = []
                self.cluster[cluster_i].append(signal)
                latest_signal[cluster_i] = signal
                cluster_i = cluster_i + 1
                continue
            for cluster in self.cluster:
                if signal[0]== latest_signal[cluster][0] and signal[2] == latest_signal[cluster][2]:
                    if abs(signal[1]-latest_signal[cluster][1]) <= 200 or abs(signal[3]-latest_signal[cluster][3]) <= 200:
                        flag = 1
                        self.cluster[cluster].append(signal)
                        latest_signal[cluster] = signal
                        break
            if flag == 0:
                self.cluster[cluster_i] = []
                self.cluster[cluster_i].append(signal)
                latest_signal[cluster_i] = signal
                cluster_i = cluster_i + 1

    def breakpoint_call(self):
        for cluster_i in self.cluster:
            if len(self.cluster[cluster_i]) == 1:
                continue
            elif len(self.cluster[cluster_i]) == 2:
                if (abs(self.cluster[cluster_i][0][1] - self.cluster[cluster_i][1][1]) <= 200) and (abs(self.cluster[cluster_i][0][3] - self.cluster[cluster_i][1][3]) <= 200):
                    continue
                if abs(self.cluster[cluster_i][0][1] - self.cluster[cluster_i][1][1]) <= 200:
                    breakpoint_insert_start = (min(self.cluster[cluster_i][0][1], self.cluster[cluster_i][1][1]), max(self.cluster[cluster_i][0][1], self.cluster[cluster_i][1][1]))
                    breakpoint_deletion_start = min(self.cluster[cluster_i][0][3], self.cluster[cluster_i][1][3])
                    breakpoint_deletion_end = max(self.cluster[cluster_i][0][3], self.cluster[cluster_i][1][3])
                    self.breakpoint_all[str(self.ind)] = {}
                    self.breakpoint_all[str(self.ind)]['deletion_chr'] = self.cluster[cluster_i][0][2]
                    self.breakpoint_all[str(self.ind)]['deletion_start'] = breakpoint_deletion_start
                    self.breakpoint_all[str(self.ind)]['deletion_end'] = breakpoint_deletion_end
                    self.breakpoint_all[str(self.ind)]['insertion_chr'] = self.cluster[cluster_i][0][0]
                    self.breakpoint_all[str(self.ind)]['insertion_start'] = breakpoint_insert_start
                    self.label_revise_readid[str(self.ind)] = {}
                    self.label_revise_readid[str(self.ind)]['chr'] = self.cluster[cluster_i][0][0]
                    self.label_revise_readid[str(self.ind)]['readid'] = []
                    self.label_revise_readid[str(self.ind)]['readid'].append(self.cluster[cluster_i][0][4])
                    self.label_revise_readid[str(self.ind)]['readid'].append(self.cluster[cluster_i][1][4])
                    self.ind = self.ind + 1
                elif abs(self.cluster[cluster_i][0][3] - self.cluster[cluster_i][1][3]) <= 200:
                    breakpoint_insert_start = (min(self.cluster[cluster_i][0][3], self.cluster[cluster_i][1][3]),
                                               max(self.cluster[cluster_i][0][3], self.cluster[cluster_i][1][3]))
                    breakpoint_deletion_start = min(self.cluster[cluster_i][0][1], self.cluster[cluster_i][1][1])
                    breakpoint_deletion_end = max(self.cluster[cluster_i][0][1], self.cluster[cluster_i][1][1])
                    self.breakpoint_all[str(self.ind)] = {}
                    self.breakpoint_all[str(self.ind)]['deletion_chr'] = self.cluster[cluster_i][0][0]
                    self.breakpoint_all[str(self.ind)]['deletion_start'] = breakpoint_deletion_start
                    self.breakpoint_all[str(self.ind)]['deletion_end'] = breakpoint_deletion_end
                    self.breakpoint_all[str(self.ind)]['insertion_chr'] = self.cluster[cluster_i][0][2]
                    self.breakpoint_all[str(self.ind)]['insertion_start'] = breakpoint_insert_start
                    self.label_revise_readid[str(self.ind)] = {}
                    self.label_revise_readid[str(self.ind)]['chr'] = self.cluster[cluster_i][0][2]
                    self.label_revise_readid[str(self.ind)]['readid'] = []
                    self.label_revise_readid[str(self.ind)]['readid'].append(self.cluster[cluster_i][0][4])
                    self.label_revise_readid[str(self.ind)]['readid'].append(self.cluster[cluster_i][1][4])
                    self.ind = self.ind + 1
            else:
                self.insertion_point_find_and_filter(cluster_i)



    def insertion_point_find_and_filter(self, cluster_i):
        cluster = {}
        cluster_number = len(cluster)
        latest_signal = []
        ind_ = 0
        for signal in self.cluster[cluster_i]:
            if cluster_number == 0:
                cluster[str(ind_)] = {}
                cluster[str(ind_)]['chr'] = signal[0]
                cluster[str(ind_)]['breakpoint'] = []
                cluster[str(ind_)]['breakpoint'].append(signal[1])
                cluster[str(ind_)]['readid'] = []
                cluster[str(ind_)]['readid'].append(signal[4])
                latest_signal[ind_] = signal[1]
                ind_ = ind_ + 1
                cluster[str(ind_)] = {}
                cluster[str(ind_)]['chr'] = signal[2]
                cluster[str(ind_)]['breakpoint'] = []
                cluster[str(ind_)]['breakpoint'].append(signal[3])
                cluster[str(ind_)]['readid'] = []
                cluster[str(ind_)]['readid'].append(signal[4])
                latest_signal[ind_] = signal[3]
                ind_ = ind_ + 1
                continue
            else:
                for cluster_ind in cluster:
                    if signal[0] == cluster[cluster_ind]['chr']:
                        if abs(signal[1] - latest_signal[int(cluster_ind)]) <= 200:
                            cluster[cluster_ind]['readid'].append(signal[4])
                            cluster[cluster_ind]['breakpoint'].append(signal[1])
                        else:
                            cluster[str(ind_)] = {}
                            cluster[str(ind_)]['chr'] = signal[0]
                            cluster[str(ind_)]['breakpoint'] = []
                            cluster[str(ind_)]['breakpoint'].append(signal[1])
                            cluster[str(ind_)]['readid'] = []
                            cluster[str(ind_)]['readid'].append(signal[4])
                            latest_signal[ind_] = signal[1]
                            ind_ = ind_ + 1
                    if signal[2] == cluster[cluster_ind]['chr']:
                        if abs(signal[3] - latest_signal[int(cluster_ind)]) <= 200:
                            cluster[cluster_ind]['readid'].append(signal[4])
                            cluster[cluster_ind]['breakpoint'].append(signal[3])
                        else:
                            cluster[str(ind_)] = {}
                            cluster[str(ind_)]['chr'] = signal[2]
                            cluster[str(ind_)]['breakpoint'] = []
                            cluster[str(ind_)]['breakpoint'].append(signal[3])
                            cluster[str(ind_)]['readid'] = []
                            cluster[str(ind_)]['readid'].append(signal[4])
                            latest_signal[ind_] = signal[3]
                            ind_ = ind_ + 1
        readid_counts = Counter(k for k,v in cluster.items() for _ in v['readid'])
        top3_cluster_keys = [k for k,v in readid_counts.most_common(3)]
        insertion_chr, deletionpoint1_chr, deletionpoint2_chr = self.find_different_number(cluster[top3_cluster_keys[0]]['chr'],cluster[top3_cluster_keys[1]]['chr'], cluster[top3_cluster_keys[2]]['chr'])
        if insertion_chr == 0:
            return
        if cluster[top3_cluster_keys[0]]['chr'] == insertion_chr:
            insertion_ind = top3_cluster_keys[0]
            if  cluster[top3_cluster_keys[1]]['breakpoint'][0] < cluster[top3_cluster_keys[2]]['breakpoint'][0]:
                deletionpoint1_ind = top3_cluster_keys[1]
                deletionpoint2_ind = top3_cluster_keys[2]
            else:
                deletionpoint1_ind = top3_cluster_keys[2]
                deletionpoint2_ind = top3_cluster_keys[1]
        elif cluster[top3_cluster_keys[1]]['chr'] == insertion_chr:
            insertion_ind = top3_cluster_keys[1]
            if  cluster[top3_cluster_keys[0]]['breakpoint'][0] < cluster[top3_cluster_keys[2]]['breakpoint'][0]:
                deletionpoint1_ind = top3_cluster_keys[0]
                deletionpoint2_ind = top3_cluster_keys[2]
            else:
                deletionpoint1_ind = top3_cluster_keys[2]
                deletionpoint2_ind = top3_cluster_keys[0]
        elif cluster[top3_cluster_keys[2]]['chr'] == insertion_chr:
            insertion_ind = top3_cluster_keys[2]
            if  cluster[top3_cluster_keys[0]]['breakpoint'][0] < cluster[top3_cluster_keys[1]]['breakpoint'][0]:
                deletionpoint1_ind = top3_cluster_keys[0]
                deletionpoint2_ind = top3_cluster_keys[1]
            else:
                deletionpoint1_ind = top3_cluster_keys[1]
                deletionpoint2_ind = top3_cluster_keys[0]

        breakpoint_insert_start = (min(cluster[insertion_ind]['breakpoint']),
                                   max(cluster[insertion_ind]['breakpoint']))
        breakpoint_deletion_start = min(cluster[deletionpoint1_ind]['breakpoint'])
        breakpoint_deletion_end = max(cluster[deletionpoint2_ind]['breakpoint'])
        self.breakpoint_all[str(self.ind)] = {}
        self.breakpoint_all[str(self.ind)]['deletion_chr'] = cluster[deletionpoint1_ind]['chr']
        self.breakpoint_all[str(self.ind)]['deletion_start'] = breakpoint_deletion_start
        self.breakpoint_all[str(self.ind)]['deletion_end'] = breakpoint_deletion_end
        self.breakpoint_all[str(self.ind)]['insertion_chr'] = cluster[insertion_ind]['chr']
        self.breakpoint_all[str(self.ind)]['insertion_start'] = breakpoint_insert_start
        self.label_revise_readid[str(self.ind)] = {}
        self.label_revise_readid[str(self.ind)]['chr'] = insertion_chr
        self.label_revise_readid[str(self.ind)]['readid'] = []
        self.label_revise_readid[str(self.ind)]['readid'] = self.label_revise_readid[str(self.ind)]['readid'] + cluster[deletionpoint1_ind]['readid']
        self.label_revise_readid[str(self.ind)]['readid'] = self.label_revise_readid[str(self.ind)]['readid'] + cluster[deletionpoint2_ind]['readid']
        self.ind = self.ind + 1






    def sort_key(self, item):
        return (item[0], item[1], item[2], item[3])

    def find_different_number(self, a, b, c):
        if a != b and a != c and b==c:
            return a, b, c
        if b != a and b !=c and a==c:
            return b, a, c
        if c != a and c != b and b==c:
            return c, a, b
        if a==b and b==c:
            return 0, 0, 0

    def run(self):
        samfile = self.read_from_mapfile()
        self.data_process(samfile)
        self.tran_sigal_exact()
        self.clustering()
        self.breakpoint_call()


