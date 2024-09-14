#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH=$PATH:$SCRIPT_DIR/src/

threads=32
kmer_length=21
output_path="./"
config_file=""

# help information
show_help() {
    echo "Usage: $0 [-t THREADS] [-o OUTPUT_PATH] -c CONFIG_FILE"
    echo
    echo "Options:"
    echo "  -t THREADS        Number of threads to use (default: 32)"
    echo "  -o OUTPUT_PATH    Path to output directory (default: './')"
    echo "  -c CONFIG_FILE    Path to the configuration file"
    echo "  -h                Display this help message"
    exit 0
}

# Parse arguments
while getopts ":t:k:o:c:h" opt; do
    case $opt in
        t)
            threads=$OPTARG
            ;;
        k)
            kmer_length=$OPTARG
            ;;
        o)
            output_path=$OPTARG
            ;;
        c)
            config_file=$OPTARG
            ;;
        h)
            show_help
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            show_help
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            show_help
            ;;
    esac
done


if [ -z "$config_file" ]; then
    echo "Error: Configuration file is required."
    show_help
fi

# print arguments
echo "Threads: $threads"
echo "Output Path: $output_path"
echo "Config File: $config_file"

# Check requirement software
software_list=("hifiasm" "winnowmap" "meryl" "jellyfish")
for software in "${software_list[@]}"; do
    if ! command -v "$software" > /dev/null 2>&1; then
        echo "Error: $software is not installed."
        exit 1
    fi
done

echo "All required software are installed."

# read config content
while IFS= read -r line; do
    [[ "$line" =~ ^# || -z "$line" ]] && continue
    eval "$line"
done < $config_file

# print config
echo "Phasing: $phasing"
echo "Reference FASTA: $reference_fa"
echo "Align PAF: $align_paf"
echo "Jellyfish Kmer: $jellyfish_kmer"
echo "HiFi Reads: $hifi_reads"
echo "HiC Reads: $hic_read1"
echo "HiC Reads: $hic_read2"

echo "Chrs: ${chrs[@]}"
echo "Starts: ${starts[@]}"
echo "Ends: ${ends[@]}"
echo "Mat Starts: ${mat_starts[@]}"
echo "Mat Ends: ${mat_ends[@]}"
echo "Pat Starts: ${pat_starts[@]}"
echo "Pat Ends: ${pat_ends[@]}"


# main processing 

# check output
if [ "$output_path" != "./" ]; then
    mkdir -p "$output_path"
fi


for ((i = 0; i < ${#chrs[@]}; i++))
do
    cd $output_path
    mkdir ${chrs[$i]}
    cd ${chrs[$i]}
    echo ${chrs[$i]}

#    if [ ! "$i" -eq 0 ]; then
        # Thread, kmer size, output directory, reference genome chromosome name, reference genome chromosome start and end, confidence P-value, reference genome sequence, read comparison to the reference genome paf, jellyfish Reference rare kmer, hifi read file
    statistic_test_combination -t $threads -k 21 -o statistic_combination ${chrs[$i]} ${starts[$i]} ${ends[$i]} 0.05 \
    $reference_fa \
    $align_paf \
    $jellyfish_kmer \
    $hifi_reads > statistic_combination.log\
#    fi
    echo "hifiasm first"
    mkdir hifiasm
    cd hifiasm
    hifiasm -t 64 -o ${chrs[$i]} ../statistic_combination/*.fasta
    awk '/^S/{print ">"$2;print $3}' ${chrs[$i]}.bp.p_utg.gfa > ${chrs[$i]}.bp.p_utg.fa
    meryl count k=15 output merylDB.utg.k15 ${chrs[$i]}.bp.p_utg.fa
    meryl print greater-than distinct=0.9998 merylDB.utg.k15 > repetitive.utg.k15.txt
    winnowmap -t 64 -W repetitive.utg.k15.txt -x map-pb ${chrs[$i]}.bp.p_utg.fa ../statistic_combination/*.fasta -o hifitoutg.paf 
    
    cd ..
    mkdir scaffolding
    cd scaffolding
    hifi_paf_link.diptig_auto.py ../hifiasm/${chrs[$i]}.bp.p_utg.fa ../hifiasm/${chrs[$i]}.bp.p_utg.gfa ../hifiasm/hifitoutg.paf hifi_paf_link.gfa > hifi_paf_link.log
    

    awk '/^S/{print ">"$2;print $3}' hifi_paf_link.gfa > hifi_paf_link.fa
    meryl count k=19 output ref.meryl.k19 $reference_fa
    meryl print greater-than distinct=0.9998 ref.meryl.k19 > ref.repetitive.k19.txt
    winnowmap -t 64 -W ref.repetitive.k19.txt -x asm5 $reference_fa hifi_paf_link.fa -o hifi_paf_linktochm13.paf
    # reference genome chromosome name, start and end, contig to reference paf, contig gfa, contig sequence, ploid, available contig sequence, sequence and direction of contig
    genetic_algorithm.py ${chrs[$i]} ${starts[$i]} ${ends[$i]} hifi_paf_linktochm13.paf hifi_paf_link.gfa hifi_paf_link.fa diploid hifi_paf_link.available.fa goal_combination.log > genetic_algorithm.log

    cd ..
    if [ "$phasing" -eq 0 ]; then
        continue
    fi
    mkdir phasing
    cd phasing
    # shores are divided according to the target genome gap starting point and ending point.
    shores.py ${chrs[$i]} ${mat_starts[$i]} ${mat_ends[$i]} ${pat_starts[$i]} ${pat_ends[$i]}
    meryl count k=15 output merylDB.hifi_paf_link.available.k15 ../scaffolding/hifi_paf_link.available.fa
    meryl print greater-than distinct=0.9998 merylDB.hifi_paf_link.available.k15 > repetitive.hifi_paf_link.available.k15.txt
    # mapping the hifi reads to reference
    winnowmap -t $threads -W repetitive.hifi_paf_link.available.k15.txt -x map-pb -o hifitohifi_paf_link.available.paf \
    ../scaffolding/hifi_paf_link.available.fa $hifi_reads
    # link the contig > link.log
    get_link.py > link.log
    cat mat_shores.fa pat_shores.fa ../scaffolding/hifi_paf_link.available.fa > mat_pat_hifi_paf_link.available.fa
    # uniquekmer of shore sequence and available contig sequence is obtained
    jellyfish count -t $threads -m 31 -s 1G -o mat_pat_hifi_paf_link.available.kmer mat_pat_hifi_paf_link.available.fa
    jellyfish dump -c -t -U 1 -o mat_pat_hifi_paf_link.available.uniquekmer mat_pat_hifi_paf_link.available.kmer
    # use uniquekmer to locate hic reads
    kmerpos -t 64 -k 31 -C -o read1.pos mat_pat_hifi_paf_link.available.uniquekmer $hic_read1 &
    kmerpos -t 64 -k 31 -C -o read2.pos mat_pat_hifi_paf_link.available.uniquekmer $hic_read2 &
    kmerpos -t 64 -k 31 -o ref.pos mat_pat_hifi_paf_link.available.uniquekmer mat_pat_hifi_paf_link.available.fa &
    wait
    # Determine the number of anchored Hi-c > result.log
    utg_hic_link.scaffold.py mat_pat_hifi_paf_link.available.fa
    # Reference genome chromosome name and start and end, sequence and direction of contig, gfa of contig, available contig sequence, available contig sequence and shore sequence, hic connection, hifi connection
    phasing_by_simulated_annealing.py ${chrs[$i]} ${starts[$i]} ${ends[$i]} ../scaffolding/goal_combination.log ../scaffolding/hifi_paf_link.gfa ../scaffolding/hifi_paf_link.available.fa mat_pat_hifi_paf_link.available.fa result.log link.log > phasing_by_simulated_annealing.log

    # Assign the two sequences to different hap
    mkdir phase_centromere
    cd phase_centromere
    ln -s ../to_be_phased_centromere.fa to_be_phased_centromere.fa
    cat ../mat_shores.fa ../pat_shores.fa to_be_phased_centromere.fa > mat_pat_centromere.fa
    jellyfish count -t 64 -m 31 -s 1G -o mat_pat_centromere.kmer mat_pat_centromere.fa
    jellyfish dump -c -t -U 1 -o mat_pat_centromere.uniquekmer mat_pat_centromere.kmer
    kmerpos -t $threads -k 31 -C -o read1.pos mat_pat_centromere.uniquekmer $hic_reads1 &
    kmerpos -t $threads -k 31 -C -o read2.pos mat_pat_centromere.uniquekmer $hic_reads2 &
    kmerpos -t $threads -k 31 -o ref.pos mat_pat_centromere.uniquekmer mat_pat_centromere.fa &
    wait
    utg_hic_link.centromere.py mat_pat_centromere.fa
    echo "phasing finished!"
    # TODO: assign location manually
done