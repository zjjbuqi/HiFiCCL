## <a name="started"></a>Getting Started

```sh
# First, install hifiasm and ensure that it is added to the environment variables (requiring g++ and zlib)  
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make
nano ~/.bashrc  
export PATH="/<your_dir>/hifiasm:$PATH"  
source ~/.bashrc

# Then, install minimap2 and hificcl (requires Python3 and the pysam package)
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
python hificcl.py -m n -o ./ -t 30 -f <Input.fasta> -r <T2T-reference.fasta>

# Assembly under the optional mode of HiFiCCL with low-coverage HiFi reads
python hificcl_primaryassemble.py -m p -o ./ -t 30 -f <Input.fasta> -r <T2T-reference.fasta> -R <Pan-reference.gfa>  

```
## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [Why HiFiCCL?](#why)
- [Usage](#use)
  - [HiFiCCL's main mode](#main)
  - [HiFiCCL's optional mode](#optional)
  - [Output files](#output)
- [Results](#results)
- [Getting Help](#help)
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

HiFiCCL also generates assembly results for each chromosome. For more details, refer to the Getting help section.

### <a name="optional"></a>HiFiCCL' optional mode

HiFiCCL can utilize not only the T2T reference genome to guide assembly, but also the pangenome graph simultaneously for assembly guidance.
```sh
python hificcl.py -m p -f HG002_5x.fasta -r CHM13v2.0.fasta -R hprc-v1.0-minigraph-chm13.gfa -o <your_dir> -t 32
```
In this mode, you need to use the `-m` option to specify the optional mode, use the `-R` option to specify the pangenome graph for assembly guidance.

### <a name="output"></a>Output files

Hifiasm generates different types of assemblies based on the input data. 
It also writes error corrected reads to the *prefix*.ec.bin binary file and
writes overlaps to *prefix*.ovlp.source.bin and *prefix*.ovlp.reverse.bin.
For more details, please see the complete [documentation][tutorial_output].

## <a name="results"></a>Results

The following table shows the statistics of several hifiasm primary assemblies assembled with v0.12:

|<sub>Dataset<sub>|<sub>Size<sub>|<sub>Cov.<sub>|<sub>Asm options<sub>|<sub>CPU time<sub>|<sub>Wall time<sub>|<sub>RAM<sub>|<sub> N50<sub>|
|:---------------|-----:|-----:|:---------------------|-------:|--------:|----:|----------------:|
|<sub>[Mouse (C57/BL6J)][mouse-data]</sub>|<sub>2.6Gb</sub> |<sub>&times;25</sub>|<sub>-t48 -l0</sub> |<sub>172.9h</sub> |<sub>4.8h</sub> |<sub>76G</sub> |<sub>21.1Mb</sub>|
|<sub>[Maize (B73)][maize-data]</sub>     |<sub>2.2Gb</sub> |<sub>&times;22</sub>|<sub>-t48 -l0</sub> |<sub>203.2h</sub> |<sub>5.1h</sub> |<sub>68G</sub> |<sub>36.7Mb</sub>|
|<sub>[Strawberry][strawberry-data]</sub> |<sub>0.8Gb</sub> |<sub>&times;36</sub>|<sub>-t48 -D10</sub>|<sub>152.7h</sub> |<sub>3.7h</sub> |<sub>91G</sub> |<sub>17.8Mb</sub>|
|<sub>[Frog][frog-data]</sub>             |<sub>9.5Gb</sub> |<sub>&times;29</sub>|<sub>-t48</sub>     |<sub>2834.3h</sub>|<sub>69.0h</sub>|<sub>463G</sub>|<sub>9.3Mb</sub>|
|<sub>[Redwood][redwood-data]</sub>       |<sub>35.6Gb</sub>|<sub>&times;28</sub>|<sub>-t80</sub>     |<sub>3890.3h</sub>|<sub>65.5h</sub>|<sub>699G</sub>|<sub>5.4Mb</sub>|
|<sub>[Human (CHM13)][CHM13-data]</sub>   |<sub>3.1Gb</sub> |<sub>&times;32</sub>|<sub>-t48 -l0</sub> |<sub>310.7h</sub> |<sub>8.2h</sub> |<sub>114G</sub>|<sub>88.9Mb</sub>|
|<sub>[Human (HG00733)][HG00733-data]</sub>|<sub>3.1Gb</sub>|<sub>&times;33</sub>|<sub>-t48</sub>     |<sub>269.1h</sub> |<sub>6.9h</sub> |<sub>135G</sub>|<sub>69.9Mb</sub>|
|<sub>[Human (HG002)][NA24385-data]</sub> |<sub>3.1Gb</sub> |<sub>&times;36</sub>|<sub>-t48</sub>     |<sub>305.4h</sub> |<sub>7.7h</sub> |<sub>137G</sub>|<sub>98.7Mb</sub>|

[mouse-data]:      https://www.ncbi.nlm.nih.gov/sra/?term=SRR11606870
[maize-data]:      https://www.ncbi.nlm.nih.gov/sra/?term=SRR11606869
[strawberry-data]: https://www.ncbi.nlm.nih.gov/sra/?term=SRR11606867
[frog-data]:       https://www.ncbi.nlm.nih.gov/sra?term=(SRR11606868)%20OR%20SRR12048570
[redwood-data]:    https://www.ncbi.nlm.nih.gov/sra/?term=SRP251156
[CHM13-data]:      https://www.ncbi.nlm.nih.gov/sra?term=(((SRR11292120)%20OR%20SRR11292121)%20OR%20SRR11292122)%20OR%20SRR11292123

Hifiasm can assemble a 3.1Gb human genome in several hours or a ~30Gb hexaploid
redwood genome in a few days on a single machine. For trio binning assembly:

|<sub>Dataset<sub>|<sub>Cov.<sub>|<sub>CPU time<sub>|<sub>Elapsed time<sub>|<sub>RAM<sub>|<sub> N50<sub>|
|:---------------|-----:|-------:|--------:|----:|----------------:|
|<sub>[HG00733][HG00733-data], [\[father\]][HG00731-data], [\[mother\]][HG00732-data]</sub>|<sub>&times;33</sub>|<sub>269.1h</sub>|<sub>6.9h</sub>|<sub>135G</sub>|<sub>35.1Mb (paternal), 34.9Mb (maternal)</sub>|
|<sub>[HG002][NA24385-data],   [\[father\]][NA24149-data], [\[mother\]][NA24143-data]</sup>|<sub>&times;36</sub>|<sub>305.4h</sub>|<sub>7.7h</sub>|<sub>137G</sub>|<sub>41.0Mb (paternal), 40.8Mb (maternal)</sub>|

<!--
|<sub>[NA12878][NA12878-data], [\[father\]][NA12891-data], [\[mother\]][NA12892-data]</sub>|<sub>&times;30</sub>|<sub>180.8h</sub>|<sub>4.9h</sub>|<sub>123G</sub>|<sub>27.7Mb (paternal), 27.0Mb (maternal)</sub>|
-->

[HG00733-data]: https://www.ebi.ac.uk/ena/data/view/ERX3831682
[HG00731-data]: https://www.ebi.ac.uk/ena/data/view/ERR3241754
[HG00732-data]: https://www.ebi.ac.uk/ena/data/view/ERR3241755
[NA24385-data]: https://www.ncbi.nlm.nih.gov/sra?term=(((SRR10382244)%20OR%20SRR10382245)%20OR%20SRR10382248)%20OR%20SRR10382249
[NA24149-data]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/HG003Run01-13262252/
[NA24143-data]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004Run01-15133132/
[NA12878-data]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/PacBio_SequelII_CCS_11kb/
[NA12891-data]: https://www.ebi.ac.uk/ena/data/view/ERR194160
[NA12892-data]: https://www.ebi.ac.uk/ena/data/view/ERR194161

Human assemblies above can be acquired [from Zenodo][zenodo-human] and
non-human ones are available [here][zenodo-nonh].

[zenodo-human]: https://zenodo.org/record/4393631
[zenodo-nonh]: https://zenodo.org/record/4393750
[unitig]: http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology
[gfa]: https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[yak]: https://github.com/lh3/yak
[tutorial]: https://hifiasm.readthedocs.io/en/latest/index.html
[tutorial_output]: https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#interpreting-output


## <a name="help"></a>Getting Help

For detailed description of options, please see [tutorial][tutorial] or `man ./hifiasm.1`. The `-h`
option of hifiasm also provides brief description of options. If you have
further questions, please raise an issue at the [issue
page](https://github.com/chhylp123/hifiasm/issues).

## <a name="limit"></a>Limitations

1. Purging haplotig duplications may introduce misassemblies.

## <a name="cite"></a>Citating Hifiasm

If you use hifiasm in your work, please cite:

> Cheng, H., Concepcion, G.T., Feng, X., Zhang, H., Li H. (2021)
> Haplotype-resolved de novo assembly using phased assembly graphs with
> hifiasm. *Nat Methods*, **18**:170-175.
> https://doi.org/10.1038/s41592-020-01056-5

> Cheng, H., Jarvis, E.D., Fedrigo, O., Koepfli, K.P., Urban, L., Gemmell, N.J., Li, H. (2022)
> Haplotype-resolved assembly of diploid genomes without parental data. 
> *Nature Biotechnology*, **40**:1332–1335.
> https://doi.org/10.1038/s41587-022-01261-x

> Cheng, H., Asri, M., Lucas, J., Koren, S., Li, H. (2024)
> Scalable telomere-to-telomere assembly for diploid and polyploid genomes with double graph. 
> *Nat Methods*, **21**:967-970.
> https://doi.org/10.1038/s41592-024-02269-8
