#!/bin/bash

# Navigate to your project root
cd recombination_rate_inference || exit 1

echo "📁 Creating clean project structure..."

# Create base folders
mkdir -p data/raw data/processed data/interim data/external notebooks src scripts/bash logs results/plots results/tables models

# Move raw input files
mv input/*.gff* input/*.fasta* input/*.txt input/*.tsv data/raw/ 2>/dev/null
mv input/msmc2_pop_size data/raw/ 2>/dev/null

# Move processed data
mv output/processed data/processed/ 2>/dev/null
mv output/raw data/processed/ 2>/dev/null

# Move external results
mv lists/blackcap_n100_* models/ 2>/dev/null
mv input/GCF_*.gz data/external/ 2>/dev/null

# Move logs
mv lists/logs logs/lists_logs 2>/dev/null
mv output/logs logs/output_logs 2>/dev/null
mv scripts/logs logs/script_logs 2>/dev/null
mv *.log *.err *.out logs/ 2>/dev/null

# Move result figures and LaTeX tables
mv output/plots results/plots/ 2>/dev/null
mv output/tables results/tables/ 2>/dev/null

# Move notebooks
mv scripts/*.ipynb notebooks/ 2>/dev/null
mv *.ipynb notebooks/ 2>/dev/null

# Move Python code to src
mv scripts/Hotspotter/*.py src/ 2>/dev/null
mv scripts/merge_msmc2_tables.py src/ 2>/dev/null
mv scripts/RecombinationRate.py src/ 2>/dev/null

# Move bash scripts
mv scripts/bash_scripts/* scripts/bash/ 2>/dev/null

# Optional: move data files like .fasta from Hotspotter
mv scripts/Hotspotter/*.fasta* data/interim/ 2>/dev/null

echo "✅ Project restructuring complete."

