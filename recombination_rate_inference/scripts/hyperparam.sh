#!/bin/bash
 
#SBATCH --job-name=PYRHO_HYPERPARAM_%J
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --error=PYRHO_HYPERPARAM_%J.err
#SBATCH --output=PYRHO_HYPERPARAM_%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pirku@evolbio.mpg.de
#SBATCH --partition=highmemnew

module load gcc/6.2.1
module load python/3.9.12

source $HOME/pyrho/bin/activate
export JOB_ID=$SLURM_JOB_ID
cd $HOME/MasterThesis/recombination_rate_inference/

pyrho hyperparam -n 100 \
  --ploidy 2 \
  --mu 4.6e-9 \
  --blockpenalty 10,20,30,50,100 \
  --windowsize 25,30,50,100 \
  --tablefile lists/blackcap_n100_N100.hdf \
  --msmc_file input/msmc2_pop_size/msmc2.txt \
  --outfile lists/blackcap_n100_hyperparam.txt \
  --num_sim 10 \
  --logfile lists/hyperparam_blackcap_n100_"${JOB_ID}".log 
