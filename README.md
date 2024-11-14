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
export PATH="$PATH:/<your_dir>/minimap2:$PATH"
 
git clone https://github.com/zjjbuqi/HiFiCCL.git  
pip install pysam
nano ~/.bashrc
export PATH="/<your_dir>/HiFiCCL:$PATH"
source ~/.bashrc

# Assembly under the main mode of HiFiCCL with low-coverage HiFi reads  
python hificcl.py -o ./ -t 30 -f <Input.fasta> -r <T2T-reference.fasta>

# Assembly under the optional mode of HiFiCCL with low-coverage HiFi reads
python hificcl.py -m p -o ./ -t 30 -f <Input.fasta> -r <T2T-reference.fasta> -R <Pan-reference.gfa>  

```
## <a name="Dependency"></a>Dependency
1. [minimap2](https://github.com/lh3/minimap2)
2. [hifiasm](https://github.com/chhylp123/hifiasm) or [flye](https://github.com/mikolmogorov/Flye) or [lja](https://github.com/AntonBankevich/LJA)
3. [pysam](https://github.com/pysam-developers/pysam)
4. [python3.10.11](https://www.python.org/downloads/release/python-31011/)

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
- [Limitations](#limit)
- [Citing Hifiasm](#cite)

## <a name="intro"></a>Introduction

With the increasing release of telomere-to-telomere (T2T) sequences and pan-genome sequences, genomics research has entered the T2T and pan-genome era, substantially advancing the field of population genomics. Sufficient-coverage HiFi data for population genomics is often prohibitively expensive, creating an urgent need for robust ultra-low coverage assemblies. However, current assemblers underperform in such conditions. Therefore, we introduce HiFiCCL, a reference-guided, chromosome-by-chromosome assembly framework for ultra-low coverage HiFi data. Our method effectively enhances ultra-low coverage assembly performance of existing assemblers, as demonstrated on two human datasets. Furthermore, combined with Hifiasm, it performs outstandingly in comparison with other state-of-the-art assemblers. Additionally, we validated HiFiCCL with Hifiasm on two plant datasets, where it also demonstrated strong performance. Meanwhile, it enhances downstream SV detection and significantly reduces inter-chromosomal misscaffoldings. Tested on 45 human datasets, HiFiCCL exhibits strong generalization and shows promising results in population genomics applications.

## <a name="why"></a>Why HiFiCCL?

* HiFiCCL improves the assembly performance of different assemblers, such as Hifiasm, HiFlye, and LJA, under ultra-low coverage conditions.

* HiFiCCL's improvement in assembly performance is also reflected in its enhancement of assembly-based SV detection, particularly in detecting challenging medically relevant SVs.

* Upon scaffolding the HiFiCCL assembly results, it was found that inter-chromosomal mis-scaffoldings were significantly reduced compared to the base assemblers.

* HiFiCCL demonstrates exceptional generalizability, as testing on 45 ultra-low coverage human datasets revealed that HiFiCCL statistically achieved better assembly quality than Hifiasm.

* At about 5x coverage, HiFiCCL runs faster than Hifiasm while using a comparable amount of memory.

## <a name="use"></a>Usage

### <a name="main"></a>HiFiCCL's main mode

A typical HiFiCCL command line looks like:
```sh
python hificcl.py -f HG002_5x.fasta -r CHM13v2.0.fasta -o <your_dir> -t 32
```
where `-f` specifies the input reads, `-r` specifies the T2T reference genome used to guide the assembly, `-t` sets the number of CPUs in use and `-o` specifies the output directory. Finally, the primary contigs are written to `output.fasta`. 

HiFiCCL uses Hifiasm as the default base assembler, but you can specify a different assembler using the `-a` option, such as `-a flye` or `-a lja`, provided that these base assemblers are already installed and added to the system path.

HiFiCCL also generates assembly results for each chromosome. For more details, refer to the [tutorial] section.

### <a name="optional"></a>HiFiCCL' optional mode

HiFiCCL can utilize not only the T2T reference genome to guide assembly, but also the pangenome graph simultaneously for assembly guidance.
```sh
python hificcl.py -m p -f HG002_5x.fasta -r CHM13v2.0.fasta -R hprc-v1.0-minigraph-chm13.gfa -o <your_dir> -t 32
```
In this mode, you need to use the `-m` option to specify the optional mode, use the `-R` option to specify the pangenome graph for assembly guidance.

### <a name="output"></a>Output files

HiFiCCL will generate the alignment information of the reads to the reference genome, which is written to *prefix*.map_to_reference.sam, and the pairwise alignment information between the reads, which is written to *prefix*.all_vs_all.paf. Additionally, it will output the assembly results for different chromosomes, as well as the merged assembly results. For more details, refer to the [tutorial] section.

## <a name="results"></a>Results

The following table shows the statistics of HiFiCCL combined with Hifiasm(0.19.5-r592).

|<sub>Dataset<sub>|<sub>Size<sub>|<sub>Cov.<sub>|<sub>Asm options<sub>|<sub>Wall time<sub>|<sub>Maximum resident set size<sub>|<sub> NG50<sub>|
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
usage: hificcl.py [-h] [-m model [str]] [-r T2T_reference_genome] [-R Pangenome_graph] -f input_reads [file] [-p prefix [str]] [-t threads [str]]
                  [--overwrite] [-a assembler [str]] [--hifiasmoption [str]] [--ljaoption [str]] [--flyeoption [str]] [--minimapoption1 [str]]
                  [--minimapoption2 [str]] [--max_hang [int]] [--int_frac [float]] [--minigraphoption [float]] [-o output_dir [str]] [-W weight [int]]
                  [-N number [int]] [--iterations [int]] [--process [int]] [-v]

options:
  -h, --help            show this help message and exit
  -m model [str]        Select the reference genome, normal or pan-genome. Please enter n or p! [n]
  -r T2T_reference_genome
                        T2T reference genome file, FASTA format.
  -R Pangenome_graph    Pan-reference genome file, GFA format. [optional]
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

Example: python hificcl.py -r <T2T_Reference.fasta> -f <your_input.fasta> -t <threads> -o <your_dir>
```
### <a name="out"></a>Output
HiFiCCL will generate the alignment information of the reads to the reference genome, which is written to `*prefix*.map_to_reference.sam`, and the pairwise alignment information between the reads, which is written to `*prefix*.all_vs_all.paf`. Additionally, it will output the assembly results for different chromosomes, as well as the merged assembly results. The `chr_by_chr_reads_initial` file stores the results of the reads aligned to the reference genome before applying the chromosome binning algorithm, while `chr_by_chr_reads` file contains the results after applying the chromosome binning algorithm. `chr*` files represent the assembly results for different chromosomes. output.fasta is the final assembly result, used for assembly evaluation.

## <a name="limit"></a>Limitations

1. The assembly accuracy of HiFiCCL becomes compromised as HiFi sequencing data coverage increases. This may be due to several factors, such as errors in chromosomal binning due to increased coverage or the impact of alignment errors.
2. HiFiCCL solely utilizes HiFi data; however, incorporating multiple data types, such as Hi-C or nanopore data, may further enhance chromosome-by-chromosome assembly performance. This approach may potentially yield superior results in medium- and high-coverage datasets and enable T2T (telomere-to-telomere) sequence assembly at lower coverage, thus enabling high-precision genetic analyses at reduced costs.
3. The optional mode of HiFiCCL performs comparably to the main mode. With the further development of pangenome construction tools and alignment tools, as well as the refinement of pangenome reference sequences, the optional mode is expected to achieve superior performance and may eventually replace the main mode, which is currently only applicable to the draft human pangenome.

## <a name="cite"></a>Citing HiFiCCL

If you use HiFiCCL in your work, please cite:
