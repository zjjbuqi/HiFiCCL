# Getting Started
\# First, install hifiasm and ensure that it is added to the environment variables (requiring g++ and zlib)  
git clone https://github.com/chhylp123/hifiasm  
cd hifiasm && make  
nano ~/.bashrc  
export PATH="$PATH:/home/username/hifiasm"  
source ~/.bashrc  

\# Then, install hificcl (requires Python3 and the pysam package)
git clone https://github.com/zjjbuqi/HiFiCCL.git
pip install pysam

\# Assemble huam
