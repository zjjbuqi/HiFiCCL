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
   
# User's Guide
