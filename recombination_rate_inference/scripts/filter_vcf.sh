#!/bin/bash
 
#SBATCH --job-name=FILTER_VCF_%J
#SBATCH --ntasks=15
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --error=FILTER_VCF_%J.err
#SBATCH --output=FILTER_VCF_%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pirku@evolbio.mpg.de
#SBATCH --partition=fast

module load gcc/6.2.1

if [[ -d "$HOME/MasterThesis" ]]; then
  cd "$HOME/MasterThesis"
else
  echo "$HOME/MasterThesis does not exist"
  exit 1
fi

VCF_FILES=(vcf/chr_*.vcf.gz)
OUTDIR="vcf/filtered_vcf"
mkdir -p "$OUTDIR"

for VCF in "${VCF_FILES[@]}"; do
  if [[ -f "$VCF" ]]; then
    OUTNAME=$(basename "$VCF" .vcf.gz)
    echo "Extracting and filtering 50 rand individuals from $VCF"
    bcftools view "$VCF" \
    --samples-file rand_n50_individuals.txt \
    --min-ac 1:minor | \
    bcftools view \
    --genotype ^miss \
    -O z -o "${OUTDIR}/n50_${OUTNAME}.vcf.gz"
  else
    echo "Error: $VCF is not a file"
  fi
done
