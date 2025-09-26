#!/bin/bash
#SBATCH --array=0-33%10  
#SBATCH --job-name=PYRHO_OPTIMIZE_A%_a%
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=94G
#SBATCH --error=PYRHO_OPTIMIZE_%A_%a.err
#SBATCH --output=PYRHO_OPTIMIZE_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pirku@evolbio.mpg.de
#SBATCH --partition=highmemnew

module load gcc/6.2.1
module load cmake/3.20.3

source "$HOME/pyrho/bin/activate"

cd "$HOME/MasterThesis" || exit 1
mkdir -p recombination_rate_inference/output/

# 1) Build a sorted array of all your phased‐VCF paths:
VCF_LIST=( vcf/filtered_vcf/*.phased.vcf.gz )

# 2) Pick exactly one VCF for this array index:
VCF="${VCF_LIST[$SLURM_ARRAY_TASK_ID]}"
if [[ ! -f "$VCF" ]]; then
  echo "[ERROR] VCF not found: $VCF" >&2
  exit 1
fi

PREFIX=$(basename "$VCF" .phased.vcf.gz)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Array ID $SLURM_ARRAY_TASK_ID → running Pyrho on $VCF with 5 threads…"

pyrho optimize \
  --vcffile "$VCF" \
  --tablefile recombination_rate_inference/lists/blackcap_n100_N100.hdf \
  --windowsize 100 \
  --blockpenalty 20 \
  --ploidy 2 \
  --numthreads 5 \
  --outfile "recombination_rate_inference/output/${PREFIX}_W100_P20.rmap"

if [[ $? -ne 0 ]]; then
  echo "[ERROR] Pyrho failed on $VCF" >&2
  exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished $VCF → output at recombination_rate_inference/output/${PREFIX}_W100_P20.rmap"
