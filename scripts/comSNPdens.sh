#!/bin/bash

# parse command line arguments
ARGS=$(getopt -o daohv: --long vcf-directory:,min-allele-count:,output:,help,verbose -- "$@")


# check if getopt returned an error
if [ $? -ne 0 ]; then
    echo "Usage: $0 [-d vcf-directory] [-a min-allele-count] [-o output_name] [-h help] [-v verbose] "
    exit 1
fi

eval set -- "$ARGS"

# process command line arguments

while true; do
    case "$1" in
        -d | --vcf-directory)
            vcf_directory="$2"
            shift 2
            ;;
        -a | --min-allele-count)
            min_allele_count="$2"
            shift 2
            ;;
        -o | --output)
            output="$2"
            shift 2
            ;;
        -h | --help)
        echo "Usage: $0 [-d vcf-directory] [-a min-allele-count] [-o output_name] [-h help] [-v verbose] "
        echo "-d, --vcf-directory    Directory path"
            echo "-a, --min-allele-count  Allele count"
            echo "-o, --output        Output file"
            echo "-v, --verbose       Enable verbose mode"
            echo "-h, --help          Show this help message"
            exit 0
            ;;
        -v | --verbose)
        verbose=1
        shift;;
        --)
        shift
        break
        ;;
        *)
        echo "invalid option: $1"
        exit 1
        ;;
    esac
done

# Output parsed values
echo "Directory: $directory"
echo "Allele count: $allele_count"
echo "Output file: $output"
if [ $verbose -eq 1 ]; then
    echo "Verbose mode enabled."
fi
    