#!/bin/bash

print_help() {
  echo "Usage: $0 [OPTIONS]"
  echo
  echo "Required options:"
  echo "  -f, --vcf-files <FILES...>        One or more VCF files (wildcards allowed)"
  echo "  -o, --output-dir <DIR>            Output directory"
  echo "  -i, --individual-list <FILE>      File with list of individuals to sample"
  echo 
  echo " -w, --window-size-bp               Size of genome window to compute the LD"
  echo
  echo "Optional:"
  echo "  -h, --help                        Show this help message and exit"
  echo
  echo "Example:"
  echo "  $0 -f data/*.vcf -o results/ -i samples.txt -w 10000" 
}

parse_args() {
  if [[ "$#" -eq 0 ]]; then
    print_help
    exit 1
  fi

  while [[ "$#" -gt 0 ]]; do
    case "$1" in
      -f|--vcf-files)
        shift
        while [[ "$#" -gt 0 && ! "$1" == -* ]]; do
          VCF_FILES+=("$1")
          shift
        done
        ;;
      -o|--output-dir)
        OUTPUT_DIR="$2"
        if [[ -z "$OUTPUT_DIR" ]]; then
          echo "Missing argument for --output-dir"
          exit 1
        fi
        shift 2
        ;;
      -i|--individual-list)
        INDV_LIST="$3"
        if [[ -z "$INDV_LIST" ]]; then
          echo "Missing argument for --individual-list"
          exit 1
        fi
        shift 2
        ;;
    -w|--window-size-bp)
        WINDOW_SIZE="$4"
        if [[ -z "$WINDOW_SIZE" ]]; then
            echo "Missing argument for --window-size-bp"
            exit 1
        fi
        shift 2
        ;;
     -h|--help)
        print_help
        exit 0
        ;;
      *)
        echo "Unknown option: $1"
        print_help
        exit 1
        ;;
    esac
  done
}

# Run the parser
parse_args "$@"

# Debug: print parsed values (optional)
echo "VCF files: ${VCF_FILES[*]}"
echo "Output directory: $OUTPUT_DIR"
echo "Individual list: $INDV_LIST"
echo "Window Size: $WINDOW_SIZE"

mkdir -p "${OUTPUT_DIR}"
for VCF in "${VCF_FILES[@]}"; do
    OUTPUT_NAME="${VCF%%.*}"
    bcftools view "$VCF" --samples-file "$INDV_LIST" \
        --min-ac 1:minor | \
        vcftools --vcf - --out "${OUTPUT_DIR}/${OUTPUT_NAME}" --geno-r2 --ld-window-bp "$WINDOW_SIZE"
done

