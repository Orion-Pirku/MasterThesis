#!/bin/bash
#BATCH --job-name=compute_snp_stats
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --error=compute_snp_stats.err
#SBATCH --output=compute_snp_stats.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pirku@evolbio.mpg.de
#SBATCH --partition=standard

module load gcc/6.2.1
module load python/python-miniconda3

cd $HOME/Masters_Thesis
mkdir -p vcf_stats

for file in vcf/sorted_*; do
	output_name=$(basename "$file" .vcf.gz)
	for i in 10000 100000 1000000; do
		vcftools --gzvcf "$file" --SNPdensity "${i}" --out "vcf_stats/snp_density_${i}_${output_name}"
		vcftools --gzvcf "$file" --window-pi "${i}" --out "vcf_stats/window_pi_${i}_${output_name}"
		vcftools --gzvcf "$file" --TajimaD "${i}" --out "vcf_stats/tajimaD_${i}_${output_name}"
	done
done
echo "done"
