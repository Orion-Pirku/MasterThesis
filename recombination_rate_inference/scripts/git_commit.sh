#!/bin/bash
 
#SBATCH --job-name=PYRHO_MAKE_TABLE_COMMIT
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --error=PYRHO_COMMIT_001.%J.err
#SBATCH --output=PYRHO_COMMIT_001.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pirku@evolbio.mpg.de
#SBATCH --partition=standard

module load gcc/6.2.1
module load python/3.6.0

cd $HOME/MasterThesis/
eval "$(ssh-agent -s)"
ssh-add $HOME/.ssh/masterThesis

git add -A
git commit -m "Auto-commit after job $SLURM_JOB_ID run on $(date)"
git push origin main


