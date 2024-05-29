import time
import sys
import math
import subprocess
import os
import re
from collections import defaultdict


def run(func, args):
    starttime = time.time()
    func(args)
    endtime = time.time()
    print('[Info] Complete!')
    timecost = endtime - starttime
    timecostd = math.floor(timecost / (60 * 60 * 24))
    timecosth = math.floor(timecost % (60 * 60 * 24) / (60 * 60))
    timecostm = math.floor(timecost % (60 * 60 * 24) % (60 * 60) / 60)
    timecosts = math.floor(timecost % (60 * 60 * 24) % (60 * 60) % 60)
    print(f'[Info] Time Cost: {timecostd}d, {timecosth}h, {timecostm}m, {timecosts}s')


def check_prerequisite(prerequisitelist: list):
    # print('Checking prerequisites...')
    prerequisitenotfound = []
    for prerequisite in prerequisitelist:
        cmd = subprocess.run(f'which {prerequisite}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if cmd.stdout == b'':
            prerequisitenotfound.append(prerequisite)
        # else:
        #     print(f'{prerequisite} located at: {cmd.stdout.decode("utf-8").strip()}')
    if prerequisitenotfound != []:
        for prerequisite in prerequisitenotfound:
            print(f'[Error] prerequisite not found: {prerequisite}')
        print(
            f'[Error] Please make sure these software have been installed, exported to $PATH, and authorized executable.')
        sys.exit(0)
    # print('All prerequisites found.')


def decompress(file):
    if 'gzip compressed data' in subprocess.run(f'file {file}', stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                                shell=True).stdout.decode():
        subprocess.run(f'gzip -d {file}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        file = '.'.join(file.split('.')[:-1])
    return file


def readFastaAsDict(fastafile):
    fastaDict = {}
    fil = open(fastafile, 'r')
    allline = fil.read()
    fil.close()
    eachidseq = allline.split('>')
    for idseq in eachidseq:
        if idseq != '':
            sidraw, seqraw = idseq.split('\n', 1)
            sid = sidraw.split()[0].strip()
            seq = seqraw.replace('\n', '').upper()
            fastaDict[sid] = seq
    return fastaDict


def reversedseq(seq: str):
    seq = seq[::-1]
    seq = seq.replace('A', 'E')
    seq = seq.replace('T', 'A')
    seq = seq.replace('E', 'T')
    seq = seq.replace('C', 'E')
    seq = seq.replace('G', 'C')
    seq = seq.replace('E', 'G')
    return seq


def changeSuffix(filename, newsuffix):
    namelist = filename.split('.')
    namelist[-1] = newsuffix
    return '.'.join(namelist)


def minimap(reffasta, qryfasta, prefix, suffix, minimapoption, overwrite, output, flag):
    if flag == 0:
        output = output + '/' + f'{prefix}.{suffix}.paf'
    elif flag == 1:
        output = output + '/' + f'{prefix}.{suffix}.sam'
    if not os.path.exists(output) or overwrite == True:
        cmdr = subprocess.run(f'minimap2 {minimapoption} -o {output} {reffasta} {qryfasta}',
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if '[morecore]' in cmdr.stderr.decode("utf-8") or cmdr.returncode < 0:
            print(f'[Error] Memory insufficient.')
            sys.exit(0)
        elif cmdr.returncode != 0:
            print(f'[Error] Unexcepted error occur in minimap2 as follow:')
            print(f'cmd: {cmdr.args}')
            print(f'returncode: {cmdr.returncode}')
            print('stdout:')
            print(cmdr.stdout.decode("utf-8"))
            print('stderr:')
            print(cmdr.stderr.decode("utf-8"))
            sys.exit(1)
    if os.path.getsize(output) == 0:
        print(f'[Error] No alignment found.')
        sys.exit(0)

    return output

def hifiasm(qryfasta, hifiasmoption, prefix,overwrite, output):
    if not os.path.exists(f'{output}.p_ctg.gfa') or overwrite == True:
        cmdr = subprocess.run(f'hifiasm {qryfasta} {hifiasmoption} -o {output}',
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if cmdr.returncode != 0:
            print(f'[Error] Unexcepted error occur in hifiasm as follow:')
            print(f'cmd: {cmdr.args}')
            print(f'returncode: {cmdr.returncode}')
            print('stdout:')
            print(cmdr.stdout.decode("utf-8"))
            print('stderr:')
            print(cmdr.stderr.decode("utf-8"))
            sys.exit(1)
    if os.path.getsize(f'{output}.p_ctg.gfa') == 0:
        print(f'[Error] No p_ctg.gfa found.')
        sys.exit(0)

def hifiasmchr(hifiasmoption, prefix,overwrite, output, file_list):
    print('[Info] Starting excuting hifiasm program for the reads of each chr...')
    for file in file_list:
        qryfasta = output + '/' + 'chr_by_chr_reads' + '/' + file + '.fasta'
        workdir = output +'/'+ file
        try:
            os.mkdir(workdir)
        except:
            pass
        os.chdir(workdir)
        workdir = os.getcwd()
        workdir = workdir + '/' + prefix
        if not os.path.exists(f'{workdir}.p_ctg.gfa') or overwrite == True:
            cmdr = subprocess.run(f'hifiasm {qryfasta} {hifiasmoption} -o {workdir}',
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if cmdr.returncode != 0:
                print(f'[Error] Unexcepted error occur in hifiasm as follow:')
                print(f'cmd: {cmdr.args}')
                print(f'returncode: {cmdr.returncode}')
                print('stdout:')
                print(cmdr.stdout.decode("utf-8"))
                print('stderr:')
                print(cmdr.stderr.decode("utf-8"))
                sys.exit(1)
        # if os.path.getsize(f'{workdir}.p_ctg.gfa') == 0:
        #     print(f'[Error] No p_ctg.gfa found.')
        #     sys.exit(0)
    print('[Info] Primary assembling completed!')



def canuchr(canuoption, prefix,overwrite, output, file_list):
    print('[Info] Starting excuting canu program for the reads of each chr...')
    for file in file_list:
        qryfasta = output + '/' + 'chr_by_chr_reads' + '/' + file + '.fasta'
        workdir = output +'/'+ file

        try:
            os.mkdir(workdir)
        except:
            pass

        os.chdir(workdir)
        workdir = os.getcwd()
        filename = 'hicanu_set'
        filepath = os.path.join(workdir, filename)
        with open(filepath, 'w') as file:
            file.write('minInputCoverage=4 stopOnLowCoverage=4\n')
        workdir = workdir + '/' + prefix

        if not os.path.exists(f'{workdir}.contigs.fasta') or overwrite == True:
            cmdr = subprocess.run(f'canu {canuoption} -p {prefix} -s {filepath} -d {workdir} {qryfasta} ',
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if cmdr.returncode != 0:
                print(f'[Error] Unexcepted error occur in canu as follow:')
                print(f'cmd: {cmdr.args}')
                print(f'returncode: {cmdr.returncode}')
                print('stdout:')
                print(cmdr.stdout.decode("utf-8"))
                print('stderr:')
                print(cmdr.stderr.decode("utf-8"))
                sys.exit(1)
        if os.path.getsize(f'{workdir}.contigs.fasta') == 0:
            print(f'[Error] No contigs.fasta found.')
    print('[Info] Primary assembling completed!')

def ljachr(ljaoption,overwrite, output, file_list):
    print('[Info] Starting excuting lja program for the reads of each chr...')
    for file in file_list:
        qryfasta = output + '/' + 'chr_by_chr_reads' + '/' + file + '.fasta'
        workdir = output +'/'+ file
        try:
            os.mkdir(workdir)
        except:
            pass
        os.chdir(workdir)
        workdir = os.getcwd()
        workdir = workdir + '/'
        if not os.path.exists(f'{workdir}assembly.fasta') or overwrite == True:
            cmdr = subprocess.run(f'lja --reads {qryfasta} {ljaoption} -o {workdir}',
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if cmdr.returncode != 0:
                print(f'[Error] Unexcepted error occur in lja as follow:')
                print(f'cmd: {cmdr.args}')
                print(f'returncode: {cmdr.returncode}')
                print('stdout:')
                print(cmdr.stdout.decode("utf-8"))
                print('stderr:')
                print(cmdr.stderr.decode("utf-8"))
                sys.exit(1)
        # if os.path.getsize(f'{workdir}.p_ctg.gfa') == 0:
        #     print(f'[Error] No p_ctg.gfa found.')
        #     sys.exit(0)
    print('[Info] Primary assembling completed!')

def verkkochr(verkkooption,overwrite, output, file_list):
    print('[Info] Starting excuting verkko program for the reads of each chr...')
    for file in file_list:
        qryfasta = output + '/' + 'chr_by_chr_reads' + '/' + file + '.fasta'
        workdir = output +'/'+ file
        try:
            os.mkdir(workdir)
        except:
            pass
        os.chdir(workdir)
        workdir = os.getcwd()
        workdir = workdir + '/'
        if not os.path.exists(f'{workdir}assembly.fasta') or overwrite == True:
            cmdr = subprocess.run(f'verkko --hifi {qryfasta} {verkkooption} -d {workdir}',
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if cmdr.returncode != 0:
                print(f'[Error] Unexcepted error occur in verkko as follow:')
                print(f'cmd: {cmdr.args}')
                print(f'returncode: {cmdr.returncode}')
                print('stdout:')
                print(cmdr.stdout.decode("utf-8"))
                print('stderr:')
                print(cmdr.stderr.decode("utf-8"))
                sys.exit(1)
        # if os.path.getsize(f'{workdir}.p_ctg.gfa') == 0:
        #     print(f'[Error] No p_ctg.gfa found.')
        #     sys.exit(0)
    print('[Info] Primary assembling completed!')

def minigraph(pangfa, qryfasta, output, overwrite, prefix, minigraphoption):

    output = output + '/' + f'{prefix}.gaf'

    if not os.path.exists(output) or overwrite == True:
        cmdr = subprocess.run(f'minigraph {minigraphoption} {pangfa} {qryfasta} > {output}',
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if cmdr.returncode != 0:
            print(f'[Error] Unexcepted error occur in hifiasm as follow:')
            print(f'cmd: {cmdr.args}')
            print(f'returncode: {cmdr.returncode}')
            print('stdout:')
            print(cmdr.stdout.decode("utf-8"))
            print('stderr:')
            print(cmdr.stderr.decode("utf-8"))
            sys.exit(1)
    if os.path.getsize(output) == 0:
        print(f'[Error] No alignment found.')
        sys.exit(0)

    return output

import sys
import glob
import os

def concatenate_gfa_files_hifiasm(prefix, output_dir):
    """
    Concatenate GFA files with a given prefix in multiple subdirectories under output_dir.

    Parameters:
    prefix (str): The prefix of the GFA files.
    output_dir (str): The directory containing the subdirectories with GFA files.

    Returns:
    str: The path to the concatenated GFA file.
    """
    gfa_files = glob.glob(f"{output_dir}/*/{prefix}.p_ctg.gfa")
    concatenated_gfa_path = f"{output_dir}/{prefix}.gfa"

    with open(concatenated_gfa_path, 'w') as outfile:
        for gfa_file in gfa_files:
            with open(gfa_file, 'r') as infile:
                outfile.write(infile.read())

    return concatenated_gfa_path

def concatenate_fasta_files_canu(prefix, output_dir):
    """
    Concatenate GFA files with a given prefix in multiple subdirectories under output_dir.

    Parameters:
    prefix (str): The prefix of the GFA files.
    output_dir (str): The directory containing the subdirectories with GFA files.

    Returns:
    str: The path to the concatenated GFA file.
    """
    fasta_files = glob.glob(f"{output_dir}/*/asm.contigs.fasta")
    concatenated_fasta_path = f"{output_dir}/{prefix}.fasta"

    with open(concatenated_fasta_path, 'w') as outfile:
        for gfa_file in fasta_files:
            with open(gfa_file, 'r') as infile:
                outfile.write(infile.read())

    return concatenated_fasta_path

def concatenate_fasta_files_lja_verkko(prefix, output_dir):
    """
    Concatenate GFA files with a given prefix in multiple subdirectories under output_dir.

    Parameters:
    prefix (str): The prefix of the GFA files.
    output_dir (str): The directory containing the subdirectories with GFA files.

    Returns:
    str: The path to the concatenated GFA file.
    """
    fasta_files = glob.glob(f"{output_dir}/*/assembly.fasta")
    concatenated_fasta_path = f"{output_dir}/{prefix}.fasta"

    with open(concatenated_fasta_path, 'w') as outfile:
        for gfa_file in fasta_files:
            with open(gfa_file, 'r') as infile:
                outfile.write(infile.read())

    return concatenated_fasta_path

def gfa_to_fasta(gfa_file, fasta_file, prefix):
    """
    Convert a GFA file to a FASTA file.

    Parameters:
    gfa_file (str): The path to the GFA file.
    fasta_file (str): The path to the output FASTA file.
    prefix (str): The prefix to prepend to sequence IDs in the FASTA file.

    Returns:
    None
    """
    with open(gfa_file, 'r') as gfa, open(fasta_file, 'w') as fasta:
        for line in gfa:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                fasta.write(f">{prefix}_{parts[1]}\n{parts[2]}\n")

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