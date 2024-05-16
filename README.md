# Getting Started
\# First, install hifiasm and ensure that it is added to the environment variables (requiring g++ and zlib)  
git clone https://github.com/chhylp123/hifiasm  
cd hifiasm && make  
nano ~/.bashrc  
export PATH="$PATH:/home/username/hifiasm"  
source ~/.bashrc  

\# Then, install minimap2 and hificcl (requires Python3 and the pysam package)  
git clone https://github.com/lh3/minimap2  
cd minimap2 && make  
nano ~/.bashrc  
export PATH="$PATH:/home/username/minimap2"  

git clone https://github.com/zjjbuqi/HiFiCCL.git  
pip install pysam  

\# Assembly under the main mode of HiFiCCL with low-coverage HiFi reads   
python hificcl_primaryassemble.py -m n -o ./ -t 30 -f read.fasta -r T2T-reference.fa  

\# Assembly under the optional mode of HiFiCCL with low-coverage HiFi reads  
python hificcl_primaryassemble.py -m p -o ./ -t 30 -f read.fasta -r T2T-reference.fa -R Pan-reference.gfa  
# Introduction
  With the release of more telomere-to-telomere sequences and pangenomic sequences for an increasing number of species, genomics research has entered the era of T2T and pangenomics, significantly advancing population genomics studies. Population genomics often requires a large number of samples and sufficient coverage, accompanied by unsustainable high sequencing costs. However, the assembly of low-coverage resequencing data can recover population-level genetic information while significantly reducing costs. Additionally, the advent of third-generation sequencing with high-fidelity data has greatly improved assembly quality. However, most current assemblers for high-fidelity data perform poorly under low coverage conditions. Therefore, we propose a chromosome-by-chromosome assembly framework for low-coverage high-fidelity resequencing reads, named HiFiCCL. For the first time, this framework utilizes telomere-telomere reference genomes or pangenome graph to guide the binning of reads by chromosome, employing a chromosome-by-chromosome assembly strategy. We compared our method against Hifiasm, LJA, Verkko, and GALA on three human and two non-human low-coverage high-fidelity datasets. The results demonstrate that our approach enhances the assembly performance for low-coverage data, achieving higher contiguity or gene completeness than the second-ranked Hifiasm.  
# User's Guide
