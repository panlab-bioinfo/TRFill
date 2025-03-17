#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH=$SCRIPT_DIR/src/:$PATH

threads=32
kmer_length=21
output_path="./"
config_file=""
phasing=0
rtype=HiFi
fmt=fq
step=0

# process fastq/fa 
convert_reads() {
    local input_file=$1
    local output_file=$2  # Default output filename

    # Create the output directory if it does not exist
    # mkdir -p "$(dirname "$output_file")"

    # Check the file format and process accordingly
    if [[ "$input_file" == *.fq.gz || "$input_file" == *.fastq.gz ]]; then
        echo "Detected input file as gzipped FASTQ format. Converting to FASTA format..."
        
        # Use zcat to decompress and convert FASTQ to FASTA using awk
        zcat "$input_file" | awk 'NR % 4 == 1 {print ">" substr($0, 2)} NR % 4 == 2 {print}' > "$output_file"
        return 1
    elif [[ "$input_file" == *.fq || "$input_file" == *.fastq ]]; then
        echo "Detected input file as FASTQ format. Converting to FASTA format..."
        # Convert FASTQ to FASTA using awk
        awk 'NR % 4 == 1 {print ">" substr($0, 2)} NR % 4 == 2 {print}' "$input_file" > "$output_file"
        return 1
    elif [[ "$input_file" == *.fa || "$input_file" == *.fasta ]]; then
        echo "Detected input file as FASTA format. No conversion needed."
        cp "$input_file" "$output_file"
        return 0
    elif [[ "$input_file" == *.bam ]]; then
        echo "Detected input file as BAM format. Converting to FASTA format using seqkit..."
        # Use seqkit to convert BAM to FASTA
        seqkit bam2fq "$input_file" -o "$output_file"
        return 1
    else
        echo "Error: Unsupported file format. Please provide a .fq, .fastq, .fq.gz, .fastq.gz, .fa, .fasta, or .bam file."
        return -1
    fi
}


# help information
show_help() {
    echo "Usage: $0 [-t THREADS] [-o OUTPUT_PATH] [-f [HiFi/ONT]] [-p] [-b] [-h] -c CONFIG_FILE"
    echo "Required:"
    echo "  -c CONFIG_FILE    Path to the configuration file"
    echo "Options:"
    echo "  -t THREADS        Number of threads to use (default: 32)"
    echo "  -o OUTPUT_PATH    Path to output directory (default: './')"
    echo "  -f Reads_format   Input format of reads for assembly [HiFi/ONT] (default: HiFi)"
    echo "  -p                When this parameter is enabled, TRFill will conduct phasing assembly for gap regions"
    echo "  -b                When the format of hifi reads is bam, this parameter is required"
    echo "  -h                Display this help message"
    exit 0
}

# Parse arguments
while getopts ":t:k:o:c:f:pbh" opt; do
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
        f)
            rtype=$OPTARG
            ;;
        p)
            phasing=1
            ;;
        b)
            fmt=bam
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

#!/bin/bash

# Function to load configuration based on phasing value
load_config() {
    local config_file=$1
    while IFS= read -r line; do
        [[ "$line" =~ ^# || -z "$line" ]] && continue
        eval "$line"
    done < "$config_file"
}

# Determine which config file to use based on phasing
if [ "$phasing" -eq 0 ]; then
    echo "Loading haploid configuration..."
else
    echo "Loading diploid configuration..."
fi
load_config $config_file

# Print configuration for verification
echo "***********input argvs*************"
echo "Phasing: $phasing"
echo "Reference FASTA: $reference"
if [ "$phasing" -eq 0 ]; then
    echo "Current assembly: $assembly"
else
    echo "Maternal assembly: $assembly_mat"
    echo "Paternal assembly: $assembly_pat"
fi
echo "Reads: $reads"
echo "HiC Reads: $hic_reads1"
echo "HiC Reads: $hic_reads2"

echo "Chrs: ${chrs[@]}"
echo "Starts: ${starts[@]}"
echo "Ends: ${ends[@]}"

if [ "$phasing" -eq 0 ]; then
    echo "Gap Starts: ${gap_starts[@]}"
    echo "Gap Ends: ${gap_ends[@]}"
    ploid="haploid"
else
    echo "Maternal Starts: ${mat_starts[@]}"
    echo "Maternal Ends: ${mat_ends[@]}"
    echo "Paternal Starts: ${pat_starts[@]}"
    echo "Paternal Ends: ${pat_ends[@]}"
    ploid="diploid"
fi

# main processing
config_full=$(readlink -f $config_file)

# check output
if [ "$output_path" != "./" ]; then
    mkdir -p "$output_path"
fi

cd $output_path

# if [ $phasing -eq 1 ]; then
#     mkdir -p hic_reads_fa
#     convert_reads $hic_reads1 hic_reads_fa/hic_R1.fa
#     convert_reads $hic_reads2 hic_reads_fa/hic_R2.fa
#     flag=$?
#     if [ $flag -eq 1 ]; then
#         hic_reads1=$PWD"/hic_reads_fa/hic_R1.fa"
#         hic_reads2=$PWD"/hic_reads_fa/hic_R2.fa"
#     fi
# fi

if [ $fmt = 'bam' ]; then
    convert_reads $reads hifi.fastq
    reads=hifi.fastq
fi

mkdir -p exact_reference 
#align by winnowmap 
if [ "$step" -ne 1 ]; then
    meryl count k=15 output reference2hifi.meryl.k15 $reference
    meryl print greater-than distinct=0.9998 reference2hifi.meryl.k15 > ref.repetitive.k15.txt
    winnowmap -t $threads -W ref.repetitive.k15.txt -x map-pb $reference $reads -o hifi2ref.paf
    exact_ref.py $reference $config_full exact_reference/cut_ref.fa
    # jellyfish
    jellyfish count -t $threads -m 21 -s 1G -o ref.21.jf exact_reference/cut_ref.fa
    jellyfish dump -c -t -U 1 -o ref.rare.21.kmer ref.21.jf
    rm ref.21.jf
fi

for ((i = 0; i < ${#chrs[@]}; i++))
do
    mkdir -p ${chrs[$i]}
    cd ${chrs[$i]}
    echo ${chrs[$i]}

    # if [ ! "$i" -eq 0 ]; then
    # Thread, kmer size, output directory, reference genome chromosome name, reference genome chromosome start and end, confidence P-value, reference genome sequence, read comparison to the reference genome paf, jellyfish Reference rare kmer, hifi read file
    statistic_test_combination -t $threads -k 21 -o statistic_combination ${chrs[$i]} ${starts[$i]} ${ends[$i]} 0.05 \
    $reference \
    ../hifi2ref.paf \
    ../ref.rare.21.kmer \
    $reads > statistic_combination.log
    # fi
    echo "hifiasm first"
    mkdir hifiasm
    cd hifiasm
    if [ "$rtype" = "HiFi" ]; then
        hifiasm -t "$threads" -o "${chrs[$i]}" ../statistic_combination/*.fasta
    elif [ "$rtype" = "ONT" ]; then
        hifiasm -t "$threads" --ont -o "${chrs[$i]}" ../statistic_combination/*.fasta
    else
        echo "Error: Unsupported read type '$rtype'"
        exit 1
    fi
    awk '/^S/{print ">"$2;print $3}' ${chrs[$i]}.bp.p_utg.gfa > ${chrs[$i]}.bp.p_utg.fa
    meryl count k=15 output merylDB.utg.k15 ${chrs[$i]}.bp.p_utg.fa
    meryl print greater-than distinct=0.9998 merylDB.utg.k15 > repetitive.utg.k15.txt
    winnowmap -t 64 -W repetitive.utg.k15.txt -x map-pb ${chrs[$i]}.bp.p_utg.fa ../statistic_combination/*.fasta -o hifitoutg.paf

    cd ..
    mkdir scaffolding
    cd scaffolding
    hifi_paf_link.diptig_auto.py ../hifiasm/${chrs[$i]}.bp.p_utg.fa ../hifiasm/${chrs[$i]}.bp.p_utg.gfa ../hifiasm/hifitoutg.paf hifi_paf_link.gfa > hifi_paf_link.log


    awk '/^S/{print ">"$2;print $3}' hifi_paf_link.gfa > hifi_paf_link.fa
    meryl count k=19 output ref.meryl.k19 $reference
    meryl print greater-than distinct=0.9998 ref.meryl.k19 > ref.repetitive.k19.txt
    winnowmap -t $threads -W ref.repetitive.k19.txt -x asm5 $reference hifi_paf_link.fa -o hifi_paf_linktochm13.paf
    # reference genome chromosome name, start and end, contig to reference paf, contig gfa, contig sequence, ploid, available contig sequence, sequence and direction of contig
    genetic_algorithm.py ${chrs[$i]} ${starts[$i]} ${ends[$i]} hifi_paf_linktochm13.paf hifi_paf_link.gfa hifi_paf_link.fa $ploid hifi_paf_link.available.fa goal_combination.log > genetic_algorithm.log
    
    if [ "$phasing" -eq 0 ]; then
        # argvs: current_assembly_fa gap_sequence gap_Chr_id gap_start gap_end
        declare -i end_dis
        echo $assembly hifi_paf_link.available.fa ${chrs[$i]} ${gap_starts[$i]} ${gap_ends[$i]}
        end_dis=$(insert.py $assembly hifi_paf_link.available.fa ${chrs[$i]} ${gap_starts[$i]} ${gap_ends[$i]})
        for ((j=$((i+1)); j<${#chrs[@]}; j++))
        do
            if [ ${chrs[$j]} == ${chrs[$i]} ]; then
                gap_starts[$j]=$((gap_starts[$j] + end_dis))
                gap_ends[$j]=$((gap_ends[$j] + end_dis))
            else
                break
            fi
        done
        cd ..
    continue
    fi

    cd ..
    mkdir phasing
    cd phasing
    # shores are divided according to the target genome gap starting point and ending point.
    shores.py ${chrs[$i]} ${mat_starts[$i]} ${mat_ends[$i]} ${pat_starts[$i]} ${pat_ends[$i]} $assembly_mat $assembly_pat
    meryl count k=15 output merylDB.hifi_paf_link.available.k15 ../scaffolding/hifi_paf_link.available.fa
    meryl print greater-than distinct=0.9998 merylDB.hifi_paf_link.available.k15 > repetitive.hifi_paf_link.available.k15.txt
    # mapping the hifi reads to reference
    winnowmap -t $threads -W repetitive.hifi_paf_link.available.k15.txt -x map-pb -o hifitohifi_paf_link.available.paf \
    ../scaffolding/hifi_paf_link.available.fa $reads
    # link the contig > link.log
    get_link.py > link.log
    cat mat_shores.fa pat_shores.fa ../scaffolding/hifi_paf_link.available.fa > mat_pat_hifi_paf_link.available.fa
    # uniquekmer of shore sequence and available contig sequence is obtained
    jellyfish count -t $threads -m 31 -s 1G -o mat_pat_hifi_paf_link.available.kmer mat_pat_hifi_paf_link.available.fa
    jellyfish dump -c -t -U 1 -o mat_pat_hifi_paf_link.available.uniquekmer mat_pat_hifi_paf_link.available.kmer
    rm mat_pat_hifi_paf_link.available.kmer
    # use uniquekmer to locate hic reads
    kmerpos -t $threads -k 31 -C -o read1.pos mat_pat_hifi_paf_link.available.uniquekmer $hic_reads1 &
    kmerpos -t $threads -k 31 -C -o read2.pos mat_pat_hifi_paf_link.available.uniquekmer $hic_reads2 &
    kmerpos -t $threads -k 31 -o ref.pos mat_pat_hifi_paf_link.available.uniquekmer mat_pat_hifi_paf_link.available.fa &
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
    jellyfish count -t $threads -m 31 -s 1G -o mat_pat_centromere.kmer mat_pat_centromere.fa
    jellyfish dump -c -t -U 1 -o mat_pat_centromere.uniquekmer mat_pat_centromere.kmer
    rm mat_pat_centromere.kmer
    kmerpos -t $threads -k 31 -C -o read1.pos mat_pat_centromere.uniquekmer $hic_reads1 &
    kmerpos -t $threads -k 31 -C -o read2.pos mat_pat_centromere.uniquekmer $hic_reads2 &
    kmerpos -t $threads -k 31 -o ref.pos mat_pat_centromere.uniquekmer mat_pat_centromere.fa &
    wait
    utg_hic_link.centromere.py mat_pat_centromere.fa


    declare -i mat_dis pat_dis
    read mat_dis pat_dis < <(insert_dip.py $assembly_mat $assembly_pat ${chrs[$i]} ${mat_starts[$i]} ${mat_ends[$i]} ${pat_starts[$i]} ${pat_ends[$i]})
    #insert_dip.py $assembly_mat $assembly_pat ${chrs[$i]} ${mat_starts[$i]} ${mat_ends[$i]} ${pat_starts[$i]} ${pat_ends[$i]}
    for ((j=$((i+1)); j<${#chrs[@]}; j++))
    do
        if [ ${chrs[$j]} == ${chrs[$i]} ]; then
            mat_starts[$j]=$((mat_starts[$j] + mat_dis))
            mat_ends[$j]=$((mat_ends[$j] + mat_dis))
            pat_starts[$j]=$((pat_starts[$j] + pat_dis))
            pat_ends[$j]=$((pat_ends[$j] + pat_dis))
        else
            break
        fi
    done
    echo "${chrs[$i]} phasing finished!"
done

echo "TRFill running finished!"