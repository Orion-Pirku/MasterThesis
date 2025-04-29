#!/bin/bash
#SBATCH --job-name=COMPUTE_SNP_STATS
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --error=compute_snp_stats.err
#SBATCH --output=compute_snp_stats.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pirku@evolbio.mpg.de
#SBATCH --partition=fast

# Load necessary modules
module load gcc/6.2.1
module load python/python-miniconda3

# Go to the working directory
cd $HOME/Masters_Thesis || { echo "Failed to change directory to $HOME/Masters_Thesis"; exit 1; }

# Create output directory if it doesn't exist
mkdir -p MasterThesis/vcf_stats

# Change into the MasterThesis directory
cd MasterThesis || { echo "Failed to change directory to MasterThesis"; exit 1; }

# Make sure vcf directory exists
VCF_DIR="$HOME/Masters_Thesis/vcf"
if [ ! -d "$VCF_DIR" ]; then
    echo "Error: Directory $VCF_DIR does not exist!"
    exit 1
fi

# Loop through sorted VCF files
for file in $VCF_DIR/sorted_*; do
    # Check if the file exists and is a regular file
    if [ ! -f "$file" ]; then
        echo "Warning: File $file does not exist or is not a regular file. Skipping..."
        continue
    fi

    # Get the base name of the VCF file
    output_name=$(basename "$file" .vcf.gz)

    # Loop over the different window sizes
    for i in 10000 100000 1000000; do
        # Run vcftools for SNP density
        vcftools --gzvcf "$file" --SNPdensity "${i}" --out "vcf_stats/snp_density_${i}_${output_name}"
        if [ $? -ne 0 ]; then
            echo "Error: vcftools failed for SNP density on $file with window size $i. Exiting..."
            exit 1
        fi

        # Run vcftools for window pi
        vcftools --gzvcf "$file" --window-pi "${i}" --out "vcf_stats/window_pi_${i}_${output_name}"
        if [ $? -ne 0 ]; then
            echo "Error: vcftools failed for window pi on $file with window size $i. Exiting..."
            exit 1
        fi

        # Run vcftools for Tajima's D
        vcftools --gzvcf "$file" --TajimaD "${i}" --out "vcf_stats/tajimaD_${i}_${output_name}"
        if [ $? -ne 0 ]; then
            echo "Error: vcftools failed for TajimaD on $file with window size $i. Exiting..."
            exit 1
        fi
    done
done

# Git operations (ensure the output directory exists before adding)
if [ ! -d "vcf_stats" ]; then
    echo "Error: Directory vcf_stats does not exist. Skipping git commit and push."
    exit 1
fi
# Notify when done
echo "SNP statistics computation completed successfully."

