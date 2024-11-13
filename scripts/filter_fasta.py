
from Bio import SeqIO

def filter_fasta(input_fasta, output_fasta, min_length):

    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if len(record.seq) > min_length:
                SeqIO.write(record, output_handle, "fasta")


if __name__ == "__main__":
    input_fasta = "/data/zjjiang/NA19240_5X/verkko/scaffold/ragtag_output/Verkko.fasta"
    output_fasta = "/data/zjjiang/NA19240_5X/verkko/scaffold/ragtag_output/Verkko_filtered.fasta"
    min_length = 5000000

    filter_fasta(input_fasta, output_fasta, min_length)