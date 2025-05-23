import sys
import hificcl_tools
import hificcl_ALG
import trans_calls

def Binningchr(args):
    reffasta, qryfasta, prefix, threads, overwrite, minimapoption1, minimapoption2, max_hang, int_frac, workdir, weight, number, iterations= args
    inputdict = hificcl_tools.readFastaAsDict(qryfasta)
    total_len = len(inputdict.keys())

    print('[Info] Starting all-vs-all alignment in minimap2...')
    alapaf_file = hificcl_tools.minimap(qryfasta, qryfasta, prefix, 'all_vs_all', minimapoption1, overwrite, workdir, 0)
    print('[Info] Starting constructing AlignmentLabelGraph...')
    ALG = hificcl_ALG.AlignmentLabelGraph(max_hang, int_frac, weight, number)
    ALG.build_ALG(inputdict, alapaf_file)
    ALG.ref_information_process(reffasta)

    if ALG.graph == {}:
        print(f'[Error] ALG.graph is NULL.')
        sys.exit(0)
    else:
        print('[Info] The AlignmentLabelGraph is established successfully.')

    print('[Info] Starting map_to_reference alignment in minimap2...')
    map_file = hificcl_tools.minimap(reffasta, qryfasta, prefix, 'map_to_reference', minimapoption2, overwrite, workdir, 1)
    print('[Info] Starting translocation detection...')
    trans_calling = trans_calls.Trans_calls(map_file)
    trans_calling.run()
    ALG.initial_label(map_file)
    ALG.overlap_degree_sort()
    N = iterations
    ALG.label_correct(N)
    file_list = ALG.write_chr_by_chr_reads(workdir, inputdict)
    return file_list, ALG.ref_size_information

