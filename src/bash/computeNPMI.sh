#!/bin/bash

set -euo pipefail
INPUT_FILE_A=""
INPUT_FILE_B=()
OUTPUT_FILE=""

print_help() {
  echo "Usage: $0 [OPTIONS]"
  echo
  echo "Required options:"
  echo "  -a, --input-file-a <FILE>              One BED input file which contains the regions of interest"
  echo
  echo "  -b, --input-files-b <FILE FILES...>    One or more BED input files to compare the regions of interest"
  echo "  -o, --output-file <FILE>               Output file"
  echo
  echo "Optional:"
  echo "  -h, --help                             Show this help message and exit"
  echo
}

parse_args() {
  if [[ "$#" -eq 0 ]]; then
    print_help
    exit 1
  fi

  while [[ "$#" -gt 0 ]]; do
    case "$1" in
    -a | --input-file-a)
      shift
      if [[ "$#" -eq 0 || "$1" == -* ]]; then
        echo "Error: Missing argument for -a|--input-file-a"
        exit 1
      fi
      INPUT_FILE_A="$1"
      shift
      ;;

    -b | --input-files-b)
      shift
      if [[ "$#" -eq 0 || "$1" == -* ]]; then
        echo "Error: Missing argument(s) for -b|--input-files-b"
        exit 1
      fi
      while [[ "$#" -gt 0 && "$1" != -* ]]; do
        INPUT_FILE_B+=("$1")
        shift
      done
      ;;

    -o | --output-file)
      shift
      if [[ "$#" -eq 0 || "$1" == -* ]]; then
        echo "Error: Missing argument for -o|--output-file"
        exit 1
      fi
      OUTPUT_FILE="$1"
      shift
      ;;

    -h | --help)
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

  if [[ -z "$INPUT_FILE_A" ]]; then
    echo "Error: -a|--input-file-a is required."
    exit 1
  fi

  if [[ ! -f "$INPUT_FILE_A" ]]; then
    echo "Error: input file '$INPUT_FILE_A' does not exist or is not a regular file."
    exit 1
  fi

  if [[ "${#INPUT_FILE_B[@]}" -eq 0 ]]; then
    echo "Error: -b|--input-files-b is required."
    exit 1
  fi

  if [[ -z "$OUTPUT_FILE" ]]; then
    echo "Error: -o|--output-file is required."
    exit 1
  fi
}

parse_args "$@"

printf "Genomic_Feature\tNPMI\tNPMI_Expected\tLower_95_CI\tUpper_95CI\n" >"$OUTPUT_FILE"

for FILE in "${INPUT_FILE_B[@]}"; do
  if [[ ! -f "$FILE" ]]; then
    echo "Error: $FILE does not exist" >&2
  fi
  echo "[$(date +'%H:%M:%S')] NPMI for $(basename "$INPUT_FILE_A") vs $(basename "$FILE")..." >&2
  cobind.py npmi \
    -b 1090000000 \
    --nameA "$(basename "$INPUT_FILE_A" '.bed')" \
    --nameB "$(basename "$FILE" '.bed')" "$INPUT_FILE_A" "$FILE" | awk '
      BEGIN{ OFS = "\t" }

        /B.name/ { feature = $2 }
        /^Coef[[:space:]]/ { npmi=$2 }
        /^Coef\(expected/ { expected=$2 }
        /^Coef\(95%/ {
          gsub(/[\[\]]/, "", $3)
          split($3, arr, ",")
          lower = arr[1]; upper = arr[2]
        }
      END {print feature, npmi, expected, lower, upper}
      ' >>"$OUTPUT_FILE"
done
