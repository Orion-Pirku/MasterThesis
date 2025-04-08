#!/bin/bash
#set -x
# parse command line arguments
ARGS=$(getopt -o d:a:hv: --long vcf-directory:,min-allele-count:,help,verbose -- "$@")


# check if getopt returned an error
if [ $? -ne 0 ]; then
    echo "Usage: $0 [-d vcf-directory] [-a min-allele-count] [-h help] [-v verbose] "
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
            shift
            ;;
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
#echo "Directory: $vcf_directory"
#echo "Allele count: $min_allele_count"

# check if directory exists

if [ -z "$vcf_directory" ]; then
    echo "Error: $vcf_directory does not exist"
    exit 1
fi

# check if vcf_directory is a directory
if [ ! -d "$vcf_directory" ]; then
    echo "Error: $vcf_directory is not a directory"
    exit 1
fi

# check if min_allele_count is a integer 
if [[ ! "$min_allele_count" =~ ^[0-9]+$ ]]; then
    echo "Error: $min_allele_count is not an integer"
    exit 1
fi

for file in "$vcf_directory"/*.vcf.gz; do
    if [ ! -f "$file" ]; then
        echo "Error: $file not found"
        continue
    fi
    output_file="$vcf_directory/sorted_$(basename "$file")"

    bcftools view $file \
        --apply-filter .,PASS \
        --types snps \
        --min-ac "$min_allele_count":minor | \
        tee -a comSNPdens.log | \
        bcftools sort - -O z -o "$output_file"
    
    if [ $? -ne 0 ]; then
        echo "Error processing file: $file"
    else
        echo "Successfully processed $file"
    fi
done
echo "files in $vcf_directory were sorted successfully"
#set +x


