#!/bin/bash

# Description: Process samples through HMMer search pipeline for Ob1 pHMM analysis

# Function to display usage
usage() {
    echo "Usage: $0 -t <threads> -h <hmm_path> [options]"
    echo ""
    echo "Required arguments:"
    echo "  -t, --threads     Number of threads to use (e.g., 32)"
    echo "  -h, --hmm-path    Path to directory containing Ob1.hmm file"
    echo ""
    echo "Optional arguments:"
    echo "  --help           Show this help message and exit"
    echo ""
    echo "Example:"
    echo "  $0 -t 32 -h /path/to/hmm/directory"
    echo ""
    echo "Note: Run this script from the directory containing sample subdirectories."
    echo "      Each sample directory should contain:"
    echo "      - sample_1.fastq.gz"
    echo "      - sample_2.fastq.gz"
    exit 1
}

# Initialize variables
threads=""
HMMpath=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -h|--hmm-path)
            HMMpath="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        *)
            echo "Error: Unknown argument '$1'"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$threads" ]]; then
    echo "Error: Number of threads (-t) is required"
    usage
fi

if [[ -z "$HMMpath" ]]; then
    echo "Error: HMM path (-h) is required"
    usage
fi

# Validate threads is a positive integer
if ! [[ "$threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Threads must be a positive integer"
    exit 1
fi

# Validate HMM path exists and contains Ob1.hmm
if [[ ! -d "$HMMpath" ]]; then
    echo "Error: HMM directory '$HMMpath' does not exist"
    exit 1
fi

if [[ ! -f "$HMMpath/Ob1.hmm" ]]; then
    echo "Error: Ob1.hmm not found in '$HMMpath'"
    exit 1
fi

# Check for required tools
required_tools=("seqkit" "hmmsearch" "pigz" "awk" "bc")
for tool in "${required_tools[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: Required tool '$tool' is not installed or not in PATH"
        exit 1
    fi
done

# Display configuration
echo "=========================================="
echo "pHMM Analysis Pipeline"
echo "=========================================="
echo "Threads: $threads"
echo "HMM Path: $HMMpath"
echo "Working Directory: $(pwd)"
echo "=========================================="

# Check if we have any sample directories
sample_count=$(find . -maxdepth 1 -type d ! -name . | wc -l)
if [[ "$sample_count" -eq 0 ]]; then
    echo "Warning: No subdirectories found in current directory"
    echo "Make sure you're running this script from the directory containing sample folders"
fi

# Initialize output file
echo -e "sample\tOb1\treads\tOb1_fxn" > pHMM.counts

# Process each sample directory
for sample_dir in */; do
    # Skip if no directories found
    [[ ! -d "$sample_dir" ]] && continue
    
    sample=${sample_dir%/}
    echo '====='
    echo "Starting on $sample"
    
    cd "$sample_dir"
    
    # Check if required input files exist
    if [[ ! -f "${sample}_1.fastq.gz" ]]; then
        echo "Warning: ${sample}_1.fastq.gz not found in $sample_dir, skipping..."
        cd ..
        continue
    fi
    
    if [[ ! -f "${sample}_2.fastq.gz" ]]; then
        echo "Warning: ${sample}_2.fastq.gz not found in $sample_dir, skipping..."
        cd ..
        continue
    fi
    

    
    echo "Translating $sample"
    seqkit replace --threads "$threads" -p "\s.+" "${sample}_1.fastq.gz" | \
        seqkit replace --threads "$threads" -p $ -r _R1 | \
        seqkit translate --threads "$threads" -f 6 -F - > "${sample}.aa"
    
    seqkit replace --threads "$threads" -p "\s.+" "${sample}_2.fastq.gz" | \
        seqkit replace --threads "$threads" -p $ -r _R2 | \
        seqkit translate --threads "$threads" -f 6 -F - >> "${sample}.aa"
    
    echo "Running HMMer on $sample"
    hmmsearch --notextw -o "${sample}.stdout" --tblout "${sample}.tblout" \
        --cpu "$threads" "$HMMpath/Ob1.hmm" "${sample}.aa"
    
    # Compress the amino acid file
    pigz --best "${sample}.aa"
    
    echo "Processing HMMer results for $sample"
    awk '!(/^#/){print $0}' "${sample}.tblout" | \
        awk '($18){split($1,header,/_frame/); split(header[1],name,/_/); print name[1]"\t"$3}' | \
        sort -k1,1 -u > "${sample}.HMMer_hits"
    
    echo "Counting Ob1 pHMM hits in $sample"
    Ob1_hits=$(grep 'COR' "${sample}.HMMer_hits" | wc -l)
    echo "$sample Ob1 pHMM count (at E-val <= 1e-2): $Ob1_hits"
    
    echo "Counting reads in $sample"
    reads=$(seqkit stats "${sample}_1.fastq.gz" | awk '{print $4}' | tail -1 | sed 's/,//g')
    
    # Calculate fraction (handle division by zero)
    if [[ "$reads" -eq 0 ]]; then
        Ob1_fxn="0.00000"
        echo "Warning: No reads found in ${sample}_1.fastq.gz"
    else
        Ob1_fxn=$(echo "scale=5; $Ob1_hits / $reads" | bc)
    fi
    
    echo -e "$sample\t$Ob1_hits\t$reads\t$Ob1_fxn" >> ../pHMM.counts
    
    cd ..
    
    echo "Finished processing $sample"
done

echo "=========================================="
echo "Pipeline completed successfully!"
echo "Results written to: pHMM.counts"
echo "=========================================="