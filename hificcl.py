import argparse
import sys
import multiprocessing
import os
import hificcl_tools
from hificcl_binningchr_normal import Binningchr
from hificcl_binningchr_pan import Binningchr_pan

def Primaryassemble(args):
    reffasta, pangfa, qryfasta, prefix, threads, overwrite, assembler, hifiasmoption, ljaoption, flyeoption, minimapoption1, minimapoption2, max_hang, int_frac, output, model, minigraphoption, \
        weight, number, iterations, process= args
    if model == 'n':
        args = reffasta, qryfasta, prefix, threads, overwrite, minimapoption1, minimapoption2, max_hang, int_frac, output, weight, number, iterations
        file_list, ref_size_information = Binningchr(args)
        file_list = [i for i in file_list if(i!='chrM') ]
        file_list_left, file_list_right = [file_list[i] for i in range(0,int(len(file_list) / 2))], \
                                          [file_list[j] for j in range(int(len(file_list) / 2),len(file_list))]
        if process == 2:
            if assembler == 'hifiasm':
                pool = multiprocessing.Pool(processes=2)
                pool.starmap(hificcl_tools.hifiasmchr, [(hifiasmoption, prefix,overwrite, output, file_list_left), (hifiasmoption, prefix,overwrite, output, file_list_right)])
                pool.close()
                pool.join()
            if assembler == 'lja':
                pool = multiprocessing.Pool(processes=2)
                pool.starmap(hificcl_tools.ljachr, [(ljaoption, overwrite, output, file_list_left),
                                                     (ljaoption, overwrite, output, file_list_right)])
                pool.close()
                pool.join()
            if assembler == 'flye':
                pool = multiprocessing.Pool(processes=2)
                pool.starmap(hificcl_tools.flyechr, [(flyeoption, overwrite, output, file_list_left, ref_size_information),
                                                       (flyeoption, overwrite, output, file_list_right, ref_size_information)])
        else:
            if assembler == 'hifiasm':
                pool = multiprocessing.Pool(processes=1)
                pool.starmap(hificcl_tools.hifiasmchr, [(hifiasmoption, prefix, overwrite, output, file_list)])
                pool.close()
                pool.join()
            if assembler == 'lja':
                pool = multiprocessing.Pool(processes=1)
                pool.starmap(hificcl_tools.ljachr, [(ljaoption, overwrite, output, file_list)])
                pool.close()
                pool.join()
            if assembler == 'flye':
                pool = multiprocessing.Pool(processes=1)
                pool.starmap(hificcl_tools.flyechr,
                             [(flyeoption, overwrite, output, file_list, ref_size_information)])

    elif model == 'p':
        args = reffasta, pangfa, qryfasta, prefix, threads, overwrite, minimapoption1, minimapoption2, max_hang, int_frac, output, minigraphoption, weight, number, iterations
        file_list, ref_size_information = Binningchr_pan(args)
        file_list = [i for i in file_list if (i != 'chrM')]
        file_list_left, file_list_right = [file_list[i] for i in range(0, int(len(file_list) / 2))], \
                                          [file_list[j] for j in range(int(len(file_list) / 2), len(file_list))]
        if process == 2:
            if assembler == 'hifiasm':
                pool = multiprocessing.Pool(processes=2)
                pool.starmap(hificcl_tools.hifiasmchr, [(hifiasmoption, prefix, overwrite, output, file_list_left),
                                                        (hifiasmoption, prefix, overwrite, output, file_list_right)])
                pool.close()
                pool.join()
            if assembler == 'lja':
                pool = multiprocessing.Pool(processes=2)
                pool.starmap(hificcl_tools.ljachr, [(ljaoption, overwrite, output, file_list_left),
                                                    (ljaoption, overwrite, output, file_list_right)])
                pool.close()
                pool.join()
            if assembler == 'flye':
                pool = multiprocessing.Pool(processes=2)
                pool.starmap(hificcl_tools.flyechr,
                             [(flyeoption, overwrite, output, file_list_left, ref_size_information),
                              (flyeoption, overwrite, output, file_list_right, ref_size_information)])
        else:
            if assembler == 'hifiasm':
                pool = multiprocessing.Pool(processes=1)
                pool.starmap(hificcl_tools.hifiasmchr, [(hifiasmoption, prefix, overwrite, output, file_list)])
                pool.close()
                pool.join()
            if assembler == 'lja':
                pool = multiprocessing.Pool(processes=1)
                pool.starmap(hificcl_tools.ljachr, [(ljaoption, overwrite, output, file_list)])
                pool.close()
                pool.join()
            if assembler == 'flye':
                pool = multiprocessing.Pool(processes=1)
                pool.starmap(hificcl_tools.flyechr,
                             [(flyeoption, overwrite, output, file_list, ref_size_information)])
    if assembler == 'hifiasm':
        merge_gfa = hificcl_tools.concatenate_gfa_files_hifiasm(prefix,output)
        merge_fasta = f"{output}/{prefix}.fasta"
        hificcl_tools.gfa_to_fasta(merge_gfa, merge_fasta, prefix)
        final_fasta = f"{output}/output.fasta"
        hificcl_tools.modify_fasta_headers(merge_fasta, final_fasta)
        os.remove(merge_fasta)
    if assembler == 'lja' or assembler == 'flye':
        merge_fasta = hificcl_tools.concatenate_fasta_files_lja_flye(prefix,output)
        final_fasta = f"{output}/output.fasta"
        hificcl_tools.modify_fasta_headers(merge_fasta, final_fasta)
        os.remove(merge_fasta)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(epilog="Example:\n python hificcl.py -r <T2T_Reference.fasta> -f <your_input.fasta> -t <threads> -o <your_dir>")
    parser.add_argument('-m', dest='model', metavar="model [str]", default='n',help='Select the reference genome, normal or pan-genome. Please enter n or p! [n]')
    parser.add_argument('-r', dest='reference_genome [file]', metavar="T2T_reference_genome",help='T2T reference genome file, FASTA format.')
    parser.add_argument('-R', dest='pan_reference_genome [file]', metavar="Pangenome_graph", help='Pan-reference genome file, GFA format. [optional]')
    parser.add_argument('-f', dest='reads', metavar="input_reads [file]", required=True, help='(*Required) raw reads file, FASTA format.')
    parser.add_argument('-p', dest='prefix', metavar="prefix [str]", default='hificcl',
                        help='The prefix used on generated files, default: hificcl')
    parser.add_argument('-t', dest='threads', metavar="threads [str]", default='10', help='Use number of threads, default: 10')
    parser.add_argument('--overwrite', dest='overwrite', action='store_true', default=False,
                        help='Overwrite existing alignment file instead of reuse.')
    parser.add_argument('-a', dest='assembler', metavar="assembler [str]", choices=['hifiasm', 'lja' ,'flye'], default='hifiasm', help='Specify assembler (support hifiasm, lja, and flye), default: hifiasm.')
    parser.add_argument('--hifiasmoption', dest='hifiasmoption', metavar="[str]", default='--primary', help='Pass additional parameters to hifiasm program for primary assembly\t[default --primary].')
    parser.add_argument('--ljaoption',dest='ljaoption', metavar="[str]", default='--diploid',help='Pass additional parameters to lja program for primary assembly\t[default --diploid].')
    parser.add_argument('--flyeoption',dest='flyeoption', metavar="[str]", default='',help='Pass additional parameters to flye program for primary\t')
    parser.add_argument('--minimapoption1', dest='minimapoption1', metavar="[str]", default='-D --dual=no --no-long-join -k19 -w5 -U50,500 --rmq -A1 -B19 -O39,81 -E3,1 -H -e0 -m100 -N20',
                        help='Pass additional parameters to minimap2 program for all-vs-all mapping.')
    parser.add_argument('--minimapoption2', dest='minimapoption2', metavar="[str]", default='-x map-hifi -a',
                        help='Pass additional parameters to minimap2 program for mapping to reference.')
    parser.add_argument('--max_hang', dest='max_hang',metavar="[int]", default=1000,
                        help='Maximum overhang length [1000]. An overhang is an unmapped region that should be mapped given a true overlap or true containment. If the overhang is too '
                             'long, the mapping is considered an internal match and will be ignored.')
    parser.add_argument('--int_frac', dest='int_frac', metavar="[float]", default=0.8,
                        help='Minimal ratio of mapping length to mapping+overhang length for a mapping considered a containment or an overlap [0.8].')
    parser.add_argument('--minigraphoption', dest='minigraphoption', metavar="[float]", default='-cx asm', help='Pass additional parameters to minigraph program for mapping to pan-reference-genome.')
    parser.add_argument('-o', dest='output', default=os.getcwd(), metavar="output_dir [str]", help='Output files path\t[default current directory]')
    parser.add_argument('-W', dest='weight', default=0.75, metavar="weight [int]", help="Lainning drop weight [0.75].")
    parser.add_argument('-N', dest="number", default=3,metavar="number [int]", help="Supporting the numebr of lainning reads [3]")
    parser.add_argument('--iterations', dest="iterations",metavar="[int]", default=200, help="Chromosome label correction rounds [200]")
    parser.add_argument('--process', dest='process',default=1, metavar="[int]", help="Number of processes used, with options of 1 or 2 [1]")
    parser.add_argument('-v', '--version',action='version',version="HiFiCCL-1.0.0", help="The version of HiFiCCL")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

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
    iterations = parser.parse_args().iterations
    process = parser.parse_args().process
    if process > 2:
        process = 2
    elif process <=1:
        process =1
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
        ljaoption = ''
        flyeoption = ''
    if assembler == 'lja':
        hificcl_tools.check_prerequisite(['lja'])
        ljaoption = parser.parse_args().ljaoption + f" -t {threads}"
        hifiasmoption = ''
        flyeoption = ''
    if assembler == 'flye':
        hificcl_tools.check_prerequisite(['flye'])
        flyeoption = f" --threads {threads}"
        hifiasmoption = ''
        ljaoption = ''
    minimapoption1 = parser.parse_args().minimapoption1 + f" -t {threads}"
    minimapoption2 = parser.parse_args().minimapoption2 + f" -t {threads}"

    args = [refgenomefile, pangenomefile, qryfile, prefix, threads, overwrite, assembler, hifiasmoption, ljaoption, flyeoption, minimapoption1, minimapoption2, max_hang, int_frac, workdir, model, minigraphoption, weight, number, iterations, process]
    print(
        f'[Info] Paramater: refgenomefile={refgenomefile}, pangenomefile={pangenomefile}, qryfile={qryfile},  prefix={prefix}, threads={threads}, overwrite={overwrite}, assembler={assembler}, hifiasmoption={hifiasmoption}, ljaoption={ljaoption}'
        f', flyeoption={flyeoption}, minimapoption1={minimapoption1}, '
        f'minimapoption2={minimapoption2}, max_hang={max_hang}, int_frac={int_frac},output={workdir}, model={model}, minigraphoption={minigraphoption}, weight={weight}, number={number}, iterations={iterations}, process={process}')
    hificcl_tools.run(Primaryassemble, args)