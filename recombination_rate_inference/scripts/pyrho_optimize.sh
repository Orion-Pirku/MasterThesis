#!/bin/bash
 
#SBATCH --job-name=PYRHO_OPTIMIZE_%J
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --error=PYRHO_OPTIMIZE_%J.err
#SBATCH --output=PYRHO_OPTIMIZE_%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pirku@evolbio.mpg.de
#SBATCH --partition=highmemnew

module load gcc/6.2.1

if [[ -x "$HOME/pyrho/bin/activate" ]]; then
  source $HOME/pyrho/bin/activate
else
  echo "Error: Virtual Environment: $HOME/pyrho/bin/activate not found"
  exit 1
fi

if [[ -d "$HOME/MasterThesis" ]]; then
  cd "$HOME/MasterThesis"
else
  echo "$HOME/MasterThesis does not exist"
  exit 1
fi

VCF_FILES=(vcf/filtered_vcf/*)
for VCF in "${VCF_FILES[@]}"; do
  if [[ -f "$VCF" ]]; then
    echo "Computing Recombination Maps for $VCF"
    OUTFILE_PREF=$(basename "$VCF" .phased.vcf.gz)
    pyrho optimize \
      --tablefile recombination_rate_inference/lists/blackcap_n100_N100.hdf \
      --windowsize 50 \
      --blockpenalty 20 \
      --ploidy 2 \
      --outfile "recombination_rate_inference/output/${OUTFILE_PREF}_W50_P20.rmap"
  else
    echo "Error: $VCF is not a file or doesn't exists"
  fi
done
