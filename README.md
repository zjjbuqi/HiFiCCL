## <a name="started"></a>Getting Started

```sh
# First, install hifiasm and ensure that it is added to the environment variables (requiring g++ and zlib)  
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make
nano ~/.bashrc  
export PATH="/<your_dir>/hifiasm:$PATH"  
source ~/.bashrc

# Then, install minimap2 and hificcl (requires Python3.10 and the pysam package)
git clone https://github.com/lh3/minimap2  
cd minimap2 && make  
nano ~/.bashrc  
export PATH="/<your_dir>/minimap2:$PATH"
 
git clone https://github.com/zjjbuqi/HiFiCCL.git  
pip install pysam
or
conda create -n hificcl python=3.10.11 pysam=0.21.0 -c conda-forge
conda activate hificcl

#[optional] if you want to use the optional mode of HiFiCCL, you also need to install minigraph and add it to the system path.
git clone https://github.com/lh3/minigraph
cd minigraph && make
nano ~/.bashrc  
export PATH="/<your_dir>/minigraph:$PATH"
source ~/.bashrc

# Assembly under the main mode of HiFiCCL with ultra-low coverage HiFi reads  
python /<your_path>/hificcl.py -o ./ -t 30 -f <absolute_path/Input.fasta> -r <absolute_path/T2T-reference.fasta>

# Assembly under the optional mode of HiFiCCL with ultra-low coverage HiFi reads
python /<your_path>/hificcl.py -m p -o ./ -t 30 -f <absolute_path/Input.fasta> -r <absolute_path/T2T-reference.fasta> -R <absolute_path/Pan-reference.gfa>  

```
## <a name="Dependency"></a>Dependency
1. [minimap2](https://github.com/lh3/minimap2)
2. [hifiasm](https://github.com/chhylp123/hifiasm) or [flye](https://github.com/mikolmogorov/Flye) or [lja](https://github.com/AntonBankevich/LJA)
3. [pysam 0.21.0](https://github.com/pysam-developers/pysam)
4. [python 3.10.11](https://www.python.org/downloads/release/python-31011/)
5. [minigraph](https://github.com/lh3/minigraph) (optional)

## Table of Contents

- [Getting Started](#started)
- [Dependency](#dependency)
- [Introduction](#intro)
- [Why HiFiCCL?](#why)
- [Usage](#use)
  - [HiFiCCL's main mode](#main)
  - [HiFiCCL's optional mode](#optional)
  - [Output files](#output)
- [Results](#results)
- [Tutorial](#help)
  - [Installation](#install)
  - [HiFiCCL's commands](#commands)
  - [Output](#out)
- [Contact](#contact)
- [Citing HiFiCCL](#cite)
- [License](#license)

## <a name="intro"></a>Introduction

Population genomics using short-read resequencing captures single nucleotide polymorphisms and small insertions and deletions but struggles with structural variants (SVs), leading to a loss of heritability in genome-wide association studies. In recent years, long-read sequencing has improved pangenome construction for key eukaryotic species, addressing this issue to some extent. Sufficient-coverage high-fidelity (HiFi) data for population genomics is often prohibitively expensive, limiting its use in large-scale populations and broader eukaryotic species and creating an urgent need for robust ultra-low coverage assemblies. However, current assemblers underperform in such conditions. To address this, we propose HiFiCCL, the first assembly framework specifically designed for ultra-low-coverage high-fidelity reads, using a reference-guided, chromosome-by chromosome assembly approach. We demonstrate that HiFiCCL improves ultra-low coverage assembly performance of existing assemblers and outperforms the state-of-the-art assemblers on human and plant datasets. Tested on 45 human datasets (~5x coverage), HiFiCCL combined with hifiasm reduces the length of misassembled contigs relative to hifiasm by an average of 21.19% and up to 38.58%. These improved assemblies enhance germline structural variant detection, reduce chromosome-level mis-scaffolding, enable more accurate pangenome graph construction, and improve the detection of rare and somatic 
structural variants based on the pangenome graph under ultra-low-coverage conditions.

## <a name="why"></a>Why HiFiCCL?

* HiFiCCL improves the assembly performance of different assemblers, such as Hifiasm, HiFlye, and LJA, under ultra-low coverage conditions.

* HiFiCCL's improvement in assembly performance is also reflected in its enhancement of assembly-based SV detection, particularly in detecting challenging medically relevant SVs.

* Upon scaffolding the HiFiCCL assembly results, it was found that inter-chromosomal mis-scaffolding were significantly reduced compared to the base assemblers.

* HiFiCCL demonstrates exceptional generalizability, as testing on 45 ultra-low coverage human datasets revealed that HiFiCCL statistically achieved better assembly quality than Hifiasm.

* At about 5x coverage human datasets, HiFiCCL runs faster than Hifiasm while using a comparable amount of memory.

## <a name="use"></a>Usage

### <a name="main"></a>HiFiCCL's main mode

A typical HiFiCCL command line looks like:
```sh
python /<your_path>/hificcl.py -f /your_dir/HG002_5x.fasta -r /your_dir/CHM13v2.0.fasta -o <your_dir> -t 32
```
where `-f` specifies the input reads, `-r` specifies the linear reference genome used to guide the assembly, `-t` sets the number of CPUs in use and `-o` specifies the output directory. Finally, the primary contigs are written to `output.fasta`. 

HiFiCCL uses Hifiasm as the default base assembler, but you can specify a different assembler using the `-a` option, such as `-a flye` or `-a lja`, provided that these base assemblers are already installed and added to the system path.

HiFiCCL also generates assembly results for each chromosome. For more details, refer to the [tutorial] section.

### <a name="optional"></a>HiFiCCL' optional mode

HiFiCCL can utilize not only the linear reference genome to guide assembly, but also the pangenome graph simultaneously for assembly guidance.
```sh
python /<your_path>/hificcl.py -m p -f /your_dir/HG002_5x.fasta -r /your_dir/CHM13v2.0.fasta -R /your_dir/hprc-v1.0-minigraph-chm13.gfa -o <your_dir> -t 32
```
In this mode, you need to use the `-m` option to specify the optional mode, use the `-R` option to specify the pangenome graph for assembly guidance.

### <a name="output"></a>Output files

HiFiCCL will generate the alignment information of the reads to the reference genome, which is written to *prefix*.map_to_reference.sam, and the pairwise alignment information between the reads, which is written to *prefix*.all_vs_all.paf. Additionally, it will output the assembly results for different chromosomes, as well as the merged assembly results. For more details, refer to the [tutorial] section.

## <a name="results"></a>Results

The following table shows the statistics of HiFiCCL combined with Hifiasm(0.19.5-r592).

|<sub>Dataset<sub>|<sub>Species size<sub>|<sub>Cov.<sub>|<sub>Asm options<sub>|<sub>Wall time<sub>|<sub>Maximum resident set size<sub>|<sub> NG50<sub>|
|:---------------|-----:|-----:|:---------------------|--------:|----:|----------------:|
|<sub>[HG002][hg002-data]</sub>|<sub>3.1Gb</sub>|<sub>&times;5</sub>|<sub>-t20 --primary</sub>|<sub>5.0h</sub>|<sub>55G</sub>|<sub>199.0Kb</sub>|
|<sub>[NA19240][na19240-data]</sub>|<sub>3.1Gb</sub>|<sub>&times;5</sub>|<sub>-t20 --primary</sub> |<sub>4.1h</sub>|<sub>53G</sub>|<sub>123.4Kb</sub>|
|<sub>[Rice][rice-data]</sub>|<sub>390Mb</sub>|<sub>&times;5</sub>|<sub>-t20 --primary</sub>|<sub>0.67h</sub>|<sub>19.0G</sub>|<sub>61.0Kb</sub>|
|<sub>[Arabidopsisthaliana][Arabidopsisthaliana-data]</sub>|<sub>135Mb</sub>|<sub>&times;5</sub>|<sub>-t20 --primary</sub>|<sub>0.2h</sub>|<sub>17.4G</sub>|<sub>90.4Kb</sub>|

[hg002-data]: https://www.ncbi.nlm.nih.gov/sra/SRR10382244
[na19240-data]: https://www.ncbi.nlm.nih.gov/sra/?term=SRR14611231
[rice-data]: https://www.ncbi.nlm.nih.gov/sra/?term=SRR11606867
[Arabidopsisthaliana-data]: https://ngdc.cncb.ac.cn/gsa/search?searchTerm=CRR573321
## <a name="help"></a>Tutorial
### <a name="install"></a>Installation
HiFiCCL relies on minimap2 for alignment and various assemblers for assembly. If you intend to use HiFiCCL with Hifiasm for assembly, you first need to visit minimap2's GitHub page [https://github.com/lh3/minimap2], install minimap2, and add it to your system path to ensure the tool can be invoked directly by typing minimap2 in the command line. Then, you need to visit Hifiasm's GitHub page [https://github.com/chhylp123/hifiasm], install Hifiasm, and similarly add it to your system path. If you wish to use other base assemblers, follow the same process.
Then, install HiFiCCL directly from GitHub. You will need to install the pysam package. The command is as follows:
```sh
git clone https://github.com/zjjbuqi/HiFiCCL.git  
pip install pysam
```
Similarly, you need to add HiFiCCL to the system path.
### <a name="commands"></a>HiFiCCL's commands
```
usage: hificcl.py [-h] [-m model [str]] [-r linear_reference_genome] [-R Pangenome_graph] -f input_reads [file] [-p prefix [str]] [-t threads [str]]
                  [--overwrite] [-a assembler [str]] [--hifiasmoption [str]] [--ljaoption [str]] [--flyeoption [str]] [--minimapoption1 [str]]
                  [--minimapoption2 [str]] [--max_hang [int]] [--int_frac [float]] [--minigraphoption [float]] [-o output_dir [str]] [-W weight [int]]
                  [-N number [int]] [--iterations [int]] [--process [int]] [-v]

options:
  -h, --help            show this help message and exit
  -m model [str]        Select the reference genome, normal or pan-genome. Please enter n or p! [n]
  -r Linear_reference_genome [file]
                        Linear reference genome file, FASTA format.
  -R Pangenome_graph [file]   Pan-reference genome file, GFA format. [optional]
  -f input_reads [file]
                        (*Required) raw reads file, FASTA format.
  -p prefix [str]       The prefix used on generated files, default: hificcl
  -t threads [str]      Use number of threads, default: 10
  --overwrite           Overwrite existing alignment file instead of reuse.
  -a assembler [str]    Specify assembler (support hifiasm, lja, and flye), default: hifiasm.
  --hifiasmoption [str]
                        Pass additional parameters to hifiasm program for primary assembly [default --primary].
  --ljaoption [str]     Pass additional parameters to lja program for primary assembly [default --diploid].
  --flyeoption [str]    Pass additional parameters to flye program for primary
  --minimapoption1 [str]
                        Pass additional parameters to minimap2 program for all-vs-all mapping.
  --minimapoption2 [str]
                        Pass additional parameters to minimap2 program for mapping to reference.
  --max_hang [int]      Maximum overhang length [1000]. An overhang is an unmapped region that should be mapped given a true overlap or true containment. If
                        the overhang is too long, the mapping is considered an internal match and will be ignored.
  --int_frac [float]    Minimal ratio of mapping length to mapping+overhang length for a mapping considered a containment or an overlap [0.8].
  --minigraphoption [float]
                        Pass additional parameters to minigraph program for mapping to pan-reference-genome.
  -o output_dir [str]   Output files path [default current directory]
  -W weight [int]       Lainning drop weight [0.75].
  -N number [int]       Supporting the numebr of lainning reads [3]
  --iterations [int]    Chromosome label correction rounds [200]
  --process [int]       Number of processes used, with options of 1 or 2 [1]
  -v, --version         The version of HiFiCCL

Example: python /<your_path>/hificcl.py -r <Linear_Reference.fasta> -f <your_input.fasta> -t <threads> -o <your_dir>
```
I hope this tool proves helpful for your research!

### <a name="out"></a>Output
HiFiCCL will generate the alignment information of the reads to the reference genome, which is written to `*prefix*.map_to_reference.sam`, and the pairwise alignment information between the reads, which is written to `*prefix*.all_vs_all.paf`. Additionally, it will output the assembly results for different chromosomes, as well as the merged assembly results. The `chr_by_chr_reads` file contains the results after applying the chromosome binning algorithm. `chr*` files represent the assembly results for different chromosomes. output.fasta is the final assembly result, used for assembly evaluation. The optional mode will also output the `*prefix*.gaf` file, which represents the alignment information of sequences to the pangenome graph.

## <a name="contact"></a>Contact

If you experience any problems or have suggestions please create an issue or a pull request.

## <a name="cite"></a>Citing HiFiCCL

If you use HiFiCCL in your work, please cite:

Jiang Z, Pan W, Gao R, Hu H, Gao W, Zhou M, et al. Reference-guided genome assembly at scale using ultra-low-coverage high-fidelity long-reads with HiFiCCL. bioRxiv. 2025:2025.04.20.649739.


## <a name="license"></a>License

The project is licensed under the MIT License.
