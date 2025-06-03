#!/bin/bash
 
#SBATCH --job-name=PYRHO_OPTIMIZE_%J
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --error=PYRHO_OPTIMIZE_%J.err
#SBATCH --output=PYRHO_OPTIMIZE_%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pirku@evolbio.mpg.de
#SBATCH --partition=highmemnew

module load gcc/6.2.1
