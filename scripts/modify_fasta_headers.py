def modify_fasta_headers(input_file, output_file):
    count = 1
    with open(input_file, 'r') as input_f, open(output_file, 'w') as output_f:
        for line in input_f:
            if line.startswith('>'):
                header = 'HG000126_{}'.format(count)
                output_f.write(header + '\n')
                count += 1
            else:
                output_f.write(line)

input_file ='/data/zjjiang/HPRC/new_assembly/HG00126/SRR29483148/hifiasm/output.fasta'
output_file = '/data/zjjiang/HPRC/new_assembly/HG00126/SRR29483148/hifiasm/output_name.fasta'
modify_fasta_headers(input_file, output_file)