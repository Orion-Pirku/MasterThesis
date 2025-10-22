#!/bin/bash

print_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

This script calculates GC content over genomic windows.

Options:
  -g, --genome-file FILE      Genome size file (required)
  -f, --fasta-file FILE       Genome FASTA file (required)
  -w, --window-size INT       Window size for calculating GC content (required)
  -o, --output-file FILE      Output file name (default: GC_content.bed)
  -h, --help                  Display this help message and exit

Example:
  $(basename "$0") -g genome.size -f genome.fa -w 1000 -o output.bed
EOF
}

# Show help if no arguments are provided
if [[ "$#" -eq 0 ]]; then
    print_help
    exit 1
fi

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in 
        -g|--genome-file)
            GENOME_SIZE="$2"
            if [[ -z "$GENOME_SIZE" || ! -f "$GENOME_SIZE" ]]; then
                echo "Error: Missing or invalid genome size file."
                exit 1
            fi
            shift 2
            ;;
        -f|--fasta-file)
            FASTA_FILE="$2"
            if [[ -z "$FASTA_FILE" || ! -f "$FASTA_FILE" ]]; then
                echo "Error: Missing or invalid genome FASTA file."
                exit 1
            fi
            shift 2
            ;;
        -w|--window-size)
            WINDOW_SIZE="$2"
            if [[ -z "$WINDOW_SIZE" ]]; then
                echo "Error: Window size missing."
                exit 1
            fi
            shift 2
            ;;
        -o|--output-file)
            OUTPUT="$2"
            if [[ -z "$OUTPUT" ]]; then
                echo "Warning: Output file name missing. Using default: GC_content.bed"
                OUTPUT="GC_content.bed"
            fi
            shift 2
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information."
            exit 1
            ;;
    esac
done

bedtools makewindows -g "$GENOME_SIZE" -w "$WINDOW_SIZE" | \
    bedtools nuc -bed - -fi "$FASTA_FILE" | \
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' > "$OUTPUT"


