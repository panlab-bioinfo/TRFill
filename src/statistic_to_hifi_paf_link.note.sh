#!/bin/bash

chrs=("chr13" "chr14" "chr15" "chr21" "chr22")
starts=(15511991 10096873 15035508 11002840 12317333) # reference start
ends=(17561453 12679517 17652018 11303875 16654016) # reference end
mat_starts=(4041615 7694226 7879462 4222934 5408063) # target genome gap start 
mat_ends=(4287078 8558319 10395169 5049388 8250405) # target genome gap end
pat_starts=(14028568 6975733 7060302 128989 1846531)
pat_ends=(15471497 8737488 9300453 463125 5603233)

for ((i = 0; i < ${#chrs[@]}; i++))
do
    cd /data/wenhuaming/data/HG002/recall/auto
    mkdir ${chrs[$i]}
    cd ${chrs[$i]}
    echo ${chrs[$i]}

    if [ ! "$i" -eq 0 ]; then
        # 线程、kmer大小、输出目录、参考基因组染色体名、参考基因组染色体起止、置信度p值、参考基因组序列、read比对到参考基因组的paf、jellyfish-reference稀有kmer、hifi read文件
        /data/wenhuaming/data/HG002/kmerpos/statistic_test_combination -t 64 -k 21 -o statistic_combination ${chrs[$i]} ${starts[$i]} ${ends[$i]} 0.05 /data/wenhuaming/data/chm13/T2Tassembly/chm13v2.0.merge.fa /data/wenhuaming/data/HG002/hifi/high/hifitochm13.paf /data/wenhuaming/data/chm13/T2Tassembly/rare21mer /data/wenhuaming/data/HG002/hifi/high/m64012_190920_173625.Q20.fastq /data/wenhuaming/data/HG002/hifi/high/m64012_190921_234837.Q20.fastq /data/wenhuaming/data/HG002/hifi/high/m64015_190920_185703.Q20.fastq /data/wenhuaming/data/HG002/hifi/high/m64015_190922_010918.Q20.fastq > statistic_combination.log
    fi

    mkdir hifiasm
    cd hifiasm
    /data/wenhuaming/software/hifiasm/hifiasm -t 64 -o ${chrs[$i]} ../statistic_combination/*.fasta
    awk '/^S/{print ">"$2;print $3}' ${chrs[$i]}.bp.p_utg.gfa > ${chrs[$i]}.bp.p_utg.fa
    /data/wenhuaming/software/winnowmap/bin/meryl count k=15 output merylDB.utg.k15 ${chrs[$i]}.bp.p_utg.fa
    /data/wenhuaming/software/winnowmap/bin/meryl print greater-than distinct=0.9998 merylDB.utg.k15 > repetitive.utg.k15.txt
    /data/wenhuaming/software/winnowmap/bin/winnowmap -t 64 -W repetitive.utg.k15.txt -x map-pb ${chrs[$i]}.bp.p_utg.fa ../statistic_combination/*.fasta -o hifitoutg.paf # 召回的hifi read比对到utg序列上
    cd ..
    mkdir scaffolding
    cd scaffolding
    # utg序列、utg的gfa、召回的hifi比对到utg上的paf、输出的gfa
    python3 /data/wenhuaming/data/HG002/recall/hifi_paf_link.diptig_auto.py ../hifiasm/${chrs[$i]}.bp.p_utg.fa ../hifiasm/${chrs[$i]}.bp.p_utg.gfa ../hifiasm/hifitoutg.paf hifi_paf_link.gfa > hifi_paf_link.log
    

    awk '/^S/{print ">"$2;print $3}' hifi_paf_link.gfa > hifi_paf_link.fa
    /data/wenhuaming/software/winnowmap/bin/winnowmap -t 64 -W /data/wenhuaming/data/chm13/T2Tassembly/repetitive_2.0_k19.txt -x asm5 /data/wenhuaming/data/chm13/T2Tassembly/chm13v2.0.merge.fa hifi_paf_link.fa -o hifi_paf_linktochm13.paf
    # 参考基因组染色体名以及起止、contig到reference的paf、contig的gfa、contig序列、ploid、输出的可用contig序列、contig的顺序与方向
    python3 /data/wenhuaming/data/HG002/recall/genetic_algorithm.py ${chrs[$i]} ${starts[$i]} ${ends[$i]} hifi_paf_linktochm13.paf hifi_paf_link.gfa hifi_paf_link.fa diploid hifi_paf_link.available.fa goal_combination.log > genetic_algorithm.log

    cd ..
    mkdir phasing
    cd phasing
    # 根据目标基因组gap起点与终点，划分岸边(shores)
    python3 /data/wenhuaming/data/HG002/recall/shores.py ${chrs[$i]} ${mat_starts[$i]} ${mat_ends[$i]} ${pat_starts[$i]} ${pat_ends[$i]}
    /data/wenhuaming/software/winnowmap/bin/meryl count k=15 output merylDB.hifi_paf_link.available.k15 ../scaffolding/hifi_paf_link.available.fa
    /data/wenhuaming/software/winnowmap/bin/meryl print greater-than distinct=0.9998 merylDB.hifi_paf_link.available.k15 > repetitive.hifi_paf_link.available.k15.txt
    # 把全基因组hifi比对到可用的contig序列
    /data/wenhuaming/software/winnowmap/bin/winnowmap -t 64 -W repetitive.hifi_paf_link.available.k15.txt -x map-pb -o hifitohifi_paf_link.available.paf ../scaffolding/hifi_paf_link.available.fa /data/wenhuaming/data/HG002/hifi/high/m64012_190920_173625.Q20.fastq /data/wenhuaming/data/HG002/hifi/high/m64012_190921_234837.Q20.fastq /data/wenhuaming/data/HG002/hifi/high/m64015_190920_185703.Q20.fastq /data/wenhuaming/data/HG002/hifi/high/m64015_190922_010918.Q20.fastq
    # 根据比对结果得到连接 > link.log
    python3 /data/wenhuaming/data/HG002/recall/get_link.py > link.log
    cat mat_shores.fa pat_shores.fa ../scaffolding/hifi_paf_link.available.fa > mat_pat_hifi_paf_link.available.fa
    # 得到岸边序列与可用contig序列的uniquekmer
    jellyfish count -t 64 -m 31 -s 1G -o mat_pat_hifi_paf_link.available.kmer mat_pat_hifi_paf_link.available.fa
    jellyfish dump -c -t -U 1 -o mat_pat_hifi_paf_link.available.uniquekmer mat_pat_hifi_paf_link.available.kmer
    # 用uniquekmer来锚定hic读段
    /data/wenhuaming/data/HG002/kmerpos/kmerpos -t 64 -k 31 -C -o read1.pos mat_pat_hifi_paf_link.available.uniquekmer /data/wenhuaming/data/HG002/hic/high/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R1_001.fasta &
    /data/wenhuaming/data/HG002/kmerpos/kmerpos -t 64 -k 31 -C -o read2.pos mat_pat_hifi_paf_link.available.uniquekmer /data/wenhuaming/data/HG002/hic/high/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R2_001.fasta &
    /data/wenhuaming/data/HG002/kmerpos/kmerpos -t 64 -k 31 -o ref.pos mat_pat_hifi_paf_link.available.uniquekmer mat_pat_hifi_paf_link.available.fa &
    wait
    # 确定锚定的hic数量 > result.log
    python3 /data/wenhuaming/data/HG002/recall/utg_hic_link.scaffold.py mat_pat_hifi_paf_link.available.fa
    # 参考基因组染色体名以及起止、contig的顺序与方向、contig的gfa、可用contig序列、可用的contig序列与岸边序列、hic连接、hifi连接
    python3 /data/wenhuaming/data/HG002/recall/phasing_by_simulated_annealing.py ${chrs[$i]} ${starts[$i]} ${ends[$i]} ../scaffolding/goal_combination.log ../scaffolding/hifi_paf_link.gfa ../scaffolding/hifi_paf_link.available.fa mat_pat_hifi_paf_link.available.fa result.log link.log > phasing_by_simulated_annealing.log
    # 现在连接成了两条序列

    # 把这两条序列用hic连接再分一次
    mkdir phase_centromere
    cd phase_centromere
    ln -s ../to_be_phased_centromere.fa to_be_phased_centromere.fa
    cat ../mat_shores.fa ../pat_shores.fa to_be_phased_centromere.fa > mat_pat_centromere.fa
    jellyfish count -t 64 -m 31 -s 1G -o mat_pat_centromere.kmer mat_pat_centromere.fa
    jellyfish dump -c -t -U 1 -o mat_pat_centromere.uniquekmer mat_pat_centromere.kmer
    # 用uniquekmer来锚定hic读段
    /data/wenhuaming/data/HG002/kmerpos/kmerpos -t 64 -k 31 -C -o read1.pos mat_pat_centromere.uniquekmer /data/wenhuaming/data/HG002/hic/high/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R1_001.fasta &
    /data/wenhuaming/data/HG002/kmerpos/kmerpos -t 64 -k 31 -C -o read2.pos mat_pat_centromere.uniquekmer /data/wenhuaming/data/HG002/hic/high/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R2_001.fasta &
    /data/wenhuaming/data/HG002/kmerpos/kmerpos -t 64 -k 31 -o ref.pos mat_pat_centromere.uniquekmer mat_pat_centromere.fa &
    wait
    # 确定锚定的hic数量
    python3 /data/wenhuaming/data/HG002/recall/utg_hic_link.centromere.py mat_pat_centromere.fa

    # TODO: 分到两个单倍型中
done