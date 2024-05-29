import argparse
import time
import sys
import math
import subprocess
import multiprocessing
import os
import re
import hificcl_tools
from hificcl_binningchr_normal import Binningchr
from hificcl_binningchr_pan import Binningchr_pan

from collections import defaultdict

### MAIN PROGRAM ###
def Primaryassemble(args):

    reffasta, pangfa, qryfasta, prefix, threads, overwrite, assembler, hifiasmoption, canuoption, ljaoption, verkkooption, minimapoption1, minimapoption2, max_hang, int_frac, output, model, minigraphoption, \
        weight, number= args
    #calling buqi_binningchr.py to classify raw reads into every chr bin.
    if model == 'n':
        args = reffasta, qryfasta, prefix, threads, overwrite, minimapoption1, minimapoption2, max_hang, int_frac, output, weight, number
        file_list = Binningchr(args)
        file_list = [i for i in file_list if(i!='chrM') ]
        file_list_left, file_list_right = [file_list[i] for i in range(0,int(len(file_list) / 2))], \
                                          [file_list[j] for j in range(int(len(file_list) / 2),len(file_list))]
        #hifiasm using multiprocessing
        if assembler == 'hifiasm':
            pool = multiprocessing.Pool(processes=2)
            pool.starmap(hificcl_tools.hifiasmchr, [(hifiasmoption, prefix,overwrite, output, file_list_left), (hifiasmoption, prefix,overwrite, output, file_list_right)])
            pool.close()
            pool.join()
        if assembler == 'canu':
            pool = multiprocessing.Pool(processes=2)
            pool.starmap(hificcl_tools.canuchr, [(canuoption, prefix, overwrite, output, file_list_left),
                                                 (canuoption, prefix, overwrite, output, file_list_right)])
            pool.close()
            pool.join()
        if assembler == 'lja':
            pool = multiprocessing.Pool(processes=2)
            pool.starmap(hificcl_tools.ljachr, [(ljaoption, overwrite, output, file_list_left),
                                                 (ljaoption, overwrite, output, file_list_right)])
            pool.close()
            pool.join()
        if assembler == 'verkko':
            pool = multiprocessing.Pool(processes=2)
            pool.starmap(hificcl_tools.verkkochr, [(verkkooption, overwrite, output, file_list_left),
                                                   (verkkooption, overwrite, output, file_list_right)])
    elif model == 'p':
        args = reffasta, pangfa, qryfasta, prefix, threads, overwrite, minimapoption1, minimapoption2, max_hang, int_frac, output, minigraphoption, weight, number
        file_list = Binningchr_pan(args)
        file_list = [i for i in file_list if (i != 'chrM')]
        file_list_left, file_list_right = [file_list[i] for i in range(0, int(len(file_list) / 2))], \
                                          [file_list[j] for j in range(int(len(file_list) / 2), len(file_list))]
        # hifiasm using multiprocessing
        if assembler == 'hifiasm':
            pool = multiprocessing.Pool(processes=2)
            pool.starmap(hificcl_tools.hifiasmchr, [(hifiasmoption, prefix, overwrite, output, file_list_left),
                                                    (hifiasmoption, prefix, overwrite, output, file_list_right)])
            pool.close()
            pool.join()
        if assembler == 'canu':
            pool = multiprocessing.Pool(processes=2)
            pool.starmap(hificcl_tools.canuchr, [(canuoption, prefix, overwrite, output, file_list_left),
                                                 (canuoption, prefix, overwrite, output, file_list_right)])
            pool.close()
            pool.join()
    if assembler == 'hifiasm':
        merge_gfa = hificcl_tools.concatenate_gfa_files_hifiasm(prefix,output)
        merge_fasta = f"{output}/{prefix}.fasta"
        hificcl_tools.gfa_to_fasta(merge_gfa, merge_fasta, prefix)
        final_fasta = f"{output}/output.fasta"
        hificcl_tools.modify_fasta_headers(merge_fasta, final_fasta)
        os.remove(merge_fasta)
    if assembler == 'canu':
        merge_fasta = hificcl_tools.concatenate_fasta_files_canu(prefix,output)
        final_fasta = f"{output}/output.fasta"
        hificcl_tools.modify_fasta_headers(merge_fasta, final_fasta)
        os.remove(merge_fasta)
    if assembler == 'lja' or assembler == 'verkko':
        merge_fasta = hificcl_tools.concatenate_fasta_files_lja_verkko(prefix,output)
        final_fasta = f"{output}/output.fasta"
        hificcl_tools.modify_fasta_headers(merge_fasta, final_fasta)
        os.remove(merge_fasta)

### RUN ###
if __name__ == '__main__':
    # Argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', dest='model', help='Select the reference genome, normal or pan-genome. Please enter n or p!')
    parser.add_argument('-r', dest='reference_genome',help='Normal reference genome file, FASTA format.')
    parser.add_argument('-R', dest='pan_reference_genome', help='Pan-reference genome file, GFA format.')
    parser.add_argument('-f', dest='reads', required=True, help='(*Required) raw reads file, FASTA format.')
    parser.add_argument('-p', dest='prefix', default='hificcl',
                        help='The prefix used on generated files, default: hificcl')
    parser.add_argument('-t', dest='threads', default='1', help='Use number of threads, default: 1')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true', default=False,
                        help='Overwrite existing alignment file instead of reuse.')
    parser.add_argument('-a', dest='assembler', choices=['hifiasm', 'canu', 'lja', 'verkko'], default='hifiasm', help='Specify assembler (support hifiasm, hicanu, lja and verkko), default: hifiasm.')
    parser.add_argument('--hifiasmoption', dest='hifiasmoption', default='--primary', help='Pass additional parameters to hifiasm program for primary alignment\t[default --primary].')
    parser.add_argument('--hicanuoption',dest='hicanuoption', default='genomeSize=3.1g -pacbio-hifi',help='Pass additional parameters to hicanu program for primary alignment\t[default genomeSize=3.1g -pacbio-hifi]')
    parser.add_argument('--ljaoption',dest='ljaoption',default='--diploid',help='Pass additional parameters to lja program for primary alignment\t[default --diploid].')
    parser.add_argument('--minimapoption1', dest='minimapoption1', default='-D --dual=no --no-long-join -k19 -w5 -U50,500 --rmq -A1 -B19 -O39,81 -E3,1 -H -e0 -m100 -N20',
                        help='Pass additional parameters to minimap2 program for all-vs-all mapping.')
    parser.add_argument('--minimapoption2', dest='minimapoption2', default='-x map-hifi -a',
                        help='Pass additional parameters to minimap2 program for mapping to reference.')
    parser.add_argument('--max_hang', dest='max_hang', default=1000,
                        help='Maximum overhang length [1000]. An overhang is an unmapped region that should be mapped given a true overlap or true containment. If the overhang is too '
                             'long, the mapping is considered an internal match and will be ignored.')
    parser.add_argument('--int_frac', dest='int_frac', default=0.8,
                        help='Minimal ratio of mapping length to mapping+overhang length for a mapping considered a containment or an overlap [0.8].')
    parser.add_argument('--minigraphoption', dest='minigraphoption',default='-cx lr', help='Pass additional parameters to minigraph program for mapping to pan-reference-genome.')
    parser.add_argument('-o', dest='output', default=os.getcwd(), help='output files path\t[default current directory]')
    parser.add_argument('-W', dest='weight', default=0.75, help="lainning drop weight [0.75].")
    parser.add_argument('-N', dest="number", default=3, help="supporting the numebr of lainning reads [3]")
    # parser.add_argument('-d', dest='depth', required=True, help='(*Required) Sequencing depth')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    # parse input paramater
    model = parser.parse_args().model
    qryfile = hificcl_tools.decompress(parser.parse_args().reads)
    prefix = parser.parse_args().prefix
    threads = parser.parse_args().threads
    overwrite = parser.parse_args().overwrite
    output = parser.parse_args().output
    # depth = parser.parse_args().depth
    assembler = parser.parse_args().assembler
    max_hang = parser.parse_args().max_hang
    int_frac = parser.parse_args().int_frac
    weight = parser.parse_args().weight
    number = parser.parse_args().number
    try:
        os.chdir(output)
    except:
        os.mkdir(output)
        os.chdir(output)
    workdir = os.getcwd()
    hificcl_tools.check_prerequisite(['minimap2'])
    if model == 'p':
        hificcl_tools.check_prerequisite(['minigraph'])
        minigraphoption = parser.parse_args().minigraphoption + f" -t {threads}"
        pangenomefile = parser.parse_args().pan_reference_genome
        refgenomefile = hificcl_tools.decompress(parser.parse_args().reference_genome)
    elif model == 'n':
        minigraphoption = None
        refgenomefile = hificcl_tools.decompress(parser.parse_args().reference_genome)
        pangenomefile = None
    else:
        print("Incorrect mode selection, please choose again.")
        sys.exit(0)
    if assembler == 'hifiasm':
        hificcl_tools.check_prerequisite(['hifiasm'])
        hifiasmoption = parser.parse_args().hifiasmoption + f" -t {threads}"
        hicanuoption = ''
        ljaoption = ''
        verkkooption = ''
    if assembler == 'canu':
        hificcl_tools.check_prerequisite(['canu'])
        hicanuoption = parser.parse_args().hicanuoption
        hifiasmoption = ''
        ljaoption = ''
        verkkooption = ''
    if assembler == 'lja':
        hificcl_tools.check_prerequisite(['lja'])
        ljaoption = parser.parse_args().ljaoption + f" -t {threads}"
        hifiasmoption = ''
        hicanuoption = ''
        verkkooption = ''
    if assembler == 'verkko':
        hificcl_tools.check_prerequisite(['verkko'])
        verkkooption =f" --threads {threads}"
        hifiasmoption = ''
        hicanuoption = ''
        ljaoption = ''
    minimapoption1 = parser.parse_args().minimapoption1 + f" -t {threads}"
    minimapoption2 = parser.parse_args().minimapoption2 + f" -t {threads}"

    # run
    args = [refgenomefile, pangenomefile, qryfile, prefix, threads, overwrite, assembler, hifiasmoption, hicanuoption, ljaoption, verkkooption, minimapoption1, minimapoption2, max_hang, int_frac, workdir, model, minigraphoption, weight, number]
    print(
        f'[Info] Paramater: refgenomefile={refgenomefile}, pangenomefile={pangenomefile}, qryfile={qryfile},  prefix={prefix}, threads={threads}, overwrite={overwrite}, assembler={assembler}, hifiasmoption={hifiasmoption}, hicanuoption={hicanuoption}, ljaoption={ljaoption},'
        f'verkkooption={verkkooption} minimapoption1={minimapoption1}, '
        f'minimapoption2={minimapoption2}, max_hang={max_hang}, int_frac={int_frac},output={workdir}, model={model}, minigraphoption={minigraphoption}, weight={weight}, number={number}')
    hificcl_tools.run(Primaryassemble, args)