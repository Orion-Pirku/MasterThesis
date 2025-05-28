#!/bin/bash
 
#SBATCH --job-name=PYRHO_MAKE_TABLE_%J
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --error=PYRHO_MAKE_TABLE_%J.err
#SBATCH --output=PYRHO_MAKE_TABLE_%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pirku@evolbio.mpg.de
#SBATCH --partition=highmemnew

module load gcc/6.2.1
module load python/3.9.12


source $HOME/pyrho/bin/activate
export JOB_ID=$SLURM_JOB_ID

cd $HOME/MastersThesis/recombination_rate_inference/

usage() {
  echo "Usage: $0 <input_filename> <output_filename> <population_size> <moran_population_size>"
  echo
  echo "Arguments:"
  echo "  input_filename   Name of the file in input/msmc2_pop_size/ (no path)"
  echo "  output_filename  Desired output file name in lists/"
  echo "  population_size  Number of haploids"
  echo "  moran_population_size moran population"
  exit 1
}

# Check for required arguments
if [ "$#" -ne 4 ]; then
  echo "Error: Exactly 4 arguments required."
  usage
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"
POP_SIZE="$3"
MPOP_SIZE="$4"

pyrho make_table --msmc_file input/msmc2_pop_size/"$INPUT_FILE" \
	-n "$POP_SIZE" \
	-N "$MPOP_SIZE" \
	--mu 4.6e-9 \
	--outfile lists/"$OUTPUT_FILE" \
	--approx \
	--decimate_rel_tol 0.1 \
	--numthreads 10 \
	--logfile "make_table_${JOB_ID}.log" \
    --verbosity 50
	
