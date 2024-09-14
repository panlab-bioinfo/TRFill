# <center>TRFill</center>

## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Others](#others)
- [Citations](#citations)
- [Contact](#contact)

## Introduction
TRFill is a genomic gap-filling tool that relies on fully assembled homologous genomes, HiFi reads, and Hi-C reads. TRFill can repair gaps in complex repetitive regions of T2T assemblies and supports diploid genotype assembly for complex region gaps. In our tests, TRFill has demonstrated excellent performance in filling gaps in the human genome HG002 and several tomato genomes.

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
reference_fa=/data/wenhuaming/data/chm13/T2Tassembly/chm13v2.0.merge.fa
align_paf=/data/wenhuaming/data/HG002/hifi/high/hifitochm13.paf
jellyfish_kmer=/data/wenhuaming/data/chm13/T2Tassembly/rare21mer

# HiFi reads
hifi_reads=/data/wenhuaming/data/HG002/hifi/high/m64015_190922_010918.Q20.fastq 

# HiC reads
hic_reads1=/data/wenhuaming/data/HG002/hic/high/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R1_001.fasta
hic_reads2=/data/wenhuaming/data/HG002/hic/high/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R2_001.fasta

# Reference chromosomes
chrs=("chr13" "chr14" "chr15" "chr21" "chr22")

# Reference start and end positions
starts=(15511991 10096873 15035508 11002840 12317333)
ends=(17561453 12679517 17652018 11303875 16654016)

# if phasing=1, you must be provide following data, otherwise the config is end here.
# phasing=1, the mat hap
mat_starts=(4041615 7694226 7879462 4222934 5408063)
mat_ends=(4287078 8558319 10395169 5049388 8250405)

# phasing=1, the pat hap
pat_starts=(14028568 6975733 7060302 128989 1846531)
pat_ends=(15471497 8737488 9300453 463125 5603233)
```  

The format of the config.txt file must follow the structure outlined above. The parameter `phasing` is a switch that determines whether the task will run the phasing step. If `phasing=1`, TRFill will enable the phasing step and the `mat_starts, mat_ends, pat_starts, pat_ends` must be provided; otherwise, it will not be executed.   

The parameter `align_paf` is result of current HiFi reads mapping to `reference_fa` produced by **[winnowmap](https://github.com/marbl/Winnowmap.git)** (you can click this link to learn about winnowmap). You can use the following command to obtain the `align_paf`:

```sh
# if your hifi reads named "current_hifi.fastq" and reference genome named "reference.fasta"
meryl count k=19 output ref.meryl.k19 reference.fasta
meryl print greater-than distinct=0.9998 ref.meryl.k19 > ref.k19.txt
winnowmap -t 64 -W ref.k19.txt -x map-pb reference.fasta current_hifi.fastq -o hifi2reference.paf
```

The parameter `jellyfish_kmer` is a rare (unique k-mer) k-mer set of reference genome produced by **[jellyfish](https://github.com/gmarcais/Jellyfish.git)**. TRFill use 21 k lengh k-mer to locate and recall the hifi reads. You can use the following command to obtain this input file. 

```
# reference genome named "reference.fasta"
jellyfish count -t $threads -m 21 -s 1G -o reference.jf reference.fasta
jellyfish dump -c -t -U 3 -o reference.rare.21mer reference.jf
```




## Others
<!-- 章节内容 -->

## Citations
<!-- 章节内容 -->

## Contact
<!-- 章节内容 -->