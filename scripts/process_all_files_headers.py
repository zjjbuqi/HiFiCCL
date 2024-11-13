import os
import glob

def modify_fasta_headers(input_file, output_file, prefix):
    count = 1
    with open(input_file, 'r') as input_f, open(output_file, 'w') as output_f:
        for line in input_f:
            if line.startswith('>'):
                header = '>{0}_sequence{1}'.format(prefix, count)
                output_f.write(header + '\n')
                count += 1
            else:
                output_f.write(line)

def process_all_files(base_directory):

    for folder in os.listdir(base_directory):
        folder_path = os.path.join(base_directory, folder)
        if os.path.isdir(folder_path):
            if folder.startswith('HG0') or folder.startswith('NA') or folder.startswith('GR'):

                subfolders = glob.glob(os.path.join(folder_path, "SR*"))
                for subfolder in subfolders:
                    input_file = os.path.join(subfolder, 'output.fasta')
                    output_file = os.path.join(subfolder, 'output_name.fasta')
                    if os.path.exists(input_file):
                        modify_fasta_headers(input_file, output_file, folder)
                    else:
                        print(f"File not found: {input_file}")

base_directory = '/data/zjjiang/HPRC/newer_assembly/'
process_all_files(base_directory)
