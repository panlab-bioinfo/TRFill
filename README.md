# <center>TRFill</center>

## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Others](#others)
- [Citations](#citations)
- [Contact](#contact)

## Introduction
TRFill is a genomic gap-filling tool that relies on fully assembled homologous genomes, HiFi reads, and Hi-C reads. TRFill can fill gaps in complex repetitive regions of T2T assemblies and supports diploid genotype assembly for complex region gaps. In our tests, TRFill has demonstrated excellent performance in filling gaps in the human genome HG002 and several tomato genomes.

## Installation
TRFill consists of submodules in C++ and Python, ultimately assembled into a pipeline using shell scripts. During its operation, it relies on **[hifiasm](https://github.com/chhylp123/hifiasm.git)**, **[winnowmap](https://github.com/marbl/Winnowmap.git)**, and **[jellyfish](https://github.com/gmarcais/Jellyfish.git)**. Despite having multiple functional submodules, TRFill can be easily obtained and installed using the following command.  

```sh
# get software
git clone panlab-bioinfo/TRFill.git
# Compile c++ programs 
cd TRFill/src & make
cd ..
# add the dictionary to environment variable
export PATH=$PATH:you_install_path/TRFill
# pilot run
trfill_main.sh
```

## Usage
### 1. Main Program Usage  

Workflow of TRFill is easy to run as follow:  

```sh
./trfill_main.sh -h
Usage: ./trfill_main.sh [-t THREADS] [-o OUTPUT_PATH] -c CONFIG_FILE

Options:
  -t THREADS        Number of threads to use (default: 32)
  -o OUTPUT_PATH    Path to output directory (default: './')
  -c CONFIG_FILE    Path to the configuration file
  -h                Display this help message
```

This is a usage sample:  
`./trfill_main.sh -t 32 -o trf_result -c config.txt`  


### 2. Input
The key parameters for TRFill are all specified in a configuration file named **config.txt** and use `-c` to input. The sample content of **config.txt** as follow:

```  
# Parameters
phasing=1

# all file must be use absolute path
reference_fa=/data/chm13/T2Tassembly/chm13v2.0.merge.fa


# The current assembly genome that must be Chromesome level genome. 
# if phasing=0, you only provide a haploid genome.
assembly=/data/wenhuaming/data/HG002/assembly/HG002.fa

# if phasing=1, you need provide two haploid chromosomes.
assembly_mat=/data/wenhuaming/data/HG002/assembly/HG002.mat.chr.fa
assembly_pat=/data/wenhuaming/data/HG002/assembly/HG002.pat.chr.fa

# HiFi reads
hifi_reads=/data/HG002/hifi/high/m64015_190922_010918.Q20.fastq 

# HiC reads. The format of these hic reads file is fasta with a single line sequence.
hic_reads1=/data/HG002/hic/high/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R1_001.fasta
hic_reads2=/data/HG002/hic/high/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R2_001.fasta

# Reference chromosomes
chrs=("chr13" "chr14" "chr15" "chr21" "chr22")

# Reference start and end positions of every chrs
# The start and end positions of the gap in "chr13" in reference genome correspond to the first indices of starts and ends respectively. The same applies to the other chromosomes. 
starts=(15511991 10096873 15035508 11002840 12317333)
ends=(17561453 12679517 17652018 11303875 16654016)


# The location of the gap area

# if phasing=0, you need to provide the positions of gaps in current assembly that the positions.
gap_starts=(4041615 7694226 7879462 4222934 5408063)

# if phasing=1, you need to provide pat and mat positions of gaps in current assembly.
# phasing=1, the mat hap
mat_starts=(4041615 7694226 7879462 4222934 5408063)
mat_ends=(4287078 8558319 10395169 5049388 8250405)

# phasing=1, the pat hap
pat_starts=(14028568 6975733 7060302 128989 1846531)
pat_ends=(15471497 8737488 9300453 463125 5603233)
```  

The format of the config.txt file must follow the structure outlined above. The parameter `phasing` is a switch that determines whether the task will run the phasing step. If `phasing=1`, TRFill will enable the phasing step and the `mat_starts, mat_ends, pat_starts, pat_ends` must be provided; otherwise, it will not be executed. 

**All file paths in contig.txt must use absolute paths.**  
For Hi-C reads, two reads files must be **fasta** format as a follow:

```
>r1_1_1
AACCTGCTGCT
>r1_2_1
CGAACGTGCTA
......
```
You can switch the fastq reads file to fasta reads file by the following command:
`sed -n '1~4s/^@/>/p;2~4p' reads.fastq > reads.fa  # For fastq`

**Chromesome of reference and index of gap**  
The indices of the parameters **chrs, starts, and ends** correspond to each other. For example, the first item in chrs is ***chr13***, the first index in starts represents the start position of a gap in ***chr13***, and the first item in ends indicates the end position of that gap in ***chr13***. The same applies to **mat_starts/ends and pat_starts/ends**. It is important to note that the reference itself has no gaps; gaps exist in the current assembly. However, the positions of these gaps in the assembly chromosomes can be mapped to the corresponding positions in the reference. 

As the article of TRFill, the coordinate boundaries (stars/ends) of the currently assembled gap above the reference can be obtained by using [Syri](https://github.com/schneebergerlab/syri.git) collinearity analysis. Users also need to provide the coordinates of the gap position in the current assembly. These coordinates correspond to the start and end coordinates on the reference assembly, derived from the collinearity of the boundaries of the gap on the current assembly. If phasing is not required (phasing=0), users only need to provide the haploid chromosome sequences that fill the gap and their corresponding position coordinates. If phasing assembly is needed (phasing=1), users must provide the coordinates of the gap on two haploid chromosomes separately, as in the example of mat_starts/mat_ends and pat_starts/pat_ends.These positions/index/coordinates correspond to positions of gap on the reference. They are obtained by collinearity analysis ([Syri](https://github.com/schneebergerlab/syri.git)). In practical application, the boundary coordinates should be extended a certain distance from both ends of the gap.

### Output
For haploid samples(phasing=0), the sequence of gap produced by TRFill is in `result/chrN/scaffolding/hifi_paf_link.available.fa` and the filled chromosome N is in `result/chrN/scaffolding/chrN_filled.fasta`  
For diploid samples (phasing=1), the two phasing sequences of gap is in `pat/mat can be found in result/chrN/phasing/to_be_phased_centromere.fa`. These two sequences will be assigned to haplotypes according to the `result/chrN/phasing/phase_centromere/result.log` and the final chromosomes filled by TRFill are in `result/chrN/phasing/phase_centromere/chrN_mat_filled.fasta` and `result/chrN/phasing/phase_centromere/chrN_pat_filled.fasta`.  
The `result.log` sample as follows:  

```sh
#result/chrN/phasing/phase_centromere/result.log
mat000002l	pat000002l	35162
mat000001l	pat000001l	15454
pat000001l	pat000002l	12408
mat000002l	pat000001l	11545
cen000002l	pat000001l	131
mat000001l	pat000002l	92
mat000001l	mat000002l	75
cen000001l	pat000001l	42
cen000002l	pat000002l	5
cen000001l	mat000002l	4
cen000001l	mat000001l	2
cen000002l	mat000002l	2
cen000001l	cen000002l	2
mat1_cen1_mat2/pat1_cen2_pat2:	142
mat1_cen2_mat2/pat1_cen1_pat2:	44

``` 

The last two lines from the file above are the two sequences assigned to the two hap scores.As shown in this example, cen1 and cen2 are two sequences. When cen1 is assigned to mat and cen2 to pat, the score is 142; conversely, when cen1 is assigned to pat and cen2 to mat, the score is 44. This indicates that the correct allocation in this example is that cen1 belongs to the mat genome and cen2 belongs to the pat genome.


## Others
Pending replenishment

## Citations
The TRFill software and correlated algorithm is published in **Jounal(unpublished)**. If you use TRFill, please cite this paper as follows:
****************************************

## Contact
This software is developed by Professor Wei-Hua Pan's team at the Shenzhen Institute of Genome Research, Chinese Academy of Agricultural Sciences. The various functional modules are implemented by Hua-Ming Wen, while the main program integration and this README are completed by Jin-Bao Yang.

If you have any questions or concerns while using the software, please submit an issue in the repository or contact us through the following methods:

### Email:  

**Prof. Pan:**  panweihua@caas.cn  
**Yang Jinbao:**  yangjinbao@caas.cn