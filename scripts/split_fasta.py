from Bio import SeqIO
def split_fasta(input_file, maternal_file, paternal_file):
    with open(input_file, "r") as infile, \
            open(maternal_file, "w") as mat_outfile, \
            open(paternal_file, "w") as pat_outfile:

        for record in SeqIO.parse(infile, "fasta"):
            if "MATERNAL" in record.id:
                SeqIO.write(record, mat_outfile, "fasta")
            elif "PATERNAL" in record.id:
                SeqIO.write(record, pat_outfile, "fasta")


if __name__ == "__main__":
    input_file = "/data/zjjiang/HG002_T2T/hg002v1.0.1.fasta"
    maternal_file = "/data/zjjiang/HG002_T2T/hg002_maternal.fasta"
    paternal_file = "/data/zjjiang/HG002_T2T/hg002_paternal.fasta"

    split_fasta(input_file, maternal_file, paternal_file)