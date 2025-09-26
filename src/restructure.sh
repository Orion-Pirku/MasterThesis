#!/bin/bash
set -e

# Create clean target directories
mkdir -p src notebooks reports configs data/raw data/processed

########################################
# Move scripts into src
########################################
for f in *.sh *.sbatch *.py; do
    if [ -f "$f" ]; then
        echo "⚙️ Moving script $f -> src/"
        mv "$f" src/
    fi
done

########################################
# Move notebooks into notebooks/
########################################
for f in *.ipynb; do
    if [ -f "$f" ]; then
        echo "📓 Moving notebook $f -> notebooks/"
        mv "$f" notebooks/
    fi
done

########################################
# Move configs into configs/
########################################
for f in *.yml *.yaml *.json; do
    if [ -f "$f" ] && [[ "$f" != "omics_requirements.yml" ]]; then
        echo "⚙️ Moving config $f -> configs/"
        mv "$f" configs/
    fi
done

########################################
# Move analysis outputs into reports/
########################################
for d in *_stats *_inference *results* *output*; do
    if [ -d "$d" ]; then
        echo "📦 Moving analysis folder $d -> reports/"
        mv "$d" reports/
    fi
done

for f in *.pdf *.png *.jpeg *.jpg *.tiff *.html; do
    if [ -f "$f" ]; then
        echo "📊 Moving report $f -> reports/"
        mv "$f" reports/
    fi
done

########################################
# Move data files into raw/processed
########################################
# Raw data
for f in *.rmap *.vcf *.sam *.gfasta *.fasta *.fa *.fastq; do
    if [ -f "$f" ]; then
        echo "🧬 Moving raw data $f -> data/raw/"
        mv "$f" data/raw/
    fi
done

# Processed data
for f in *.tsv *.csv *.txt; do
    if [ -f "$f" ]; then
        echo "📂 Moving processed data $f -> data/processed/"
        mv "$f" data/processed/
    fi
done

########################################
# Normalize requirements at root
########################################
[ -f "data_analytics_requirements.txt" ] && mv data_analytics_requirements.txt requirements.txt
[ -f "omics_requirements.yml" ] && mv omics_requirements.yml environment.yml

########################################
# Gitignore setup
########################################
[ ! -f ".gitignore" ] && touch .gitignore

# Ensure data/ is ignored
grep -q "^data/" .gitignore || echo -e "\n# Ignore data directory\ndata/" >> .gitignore

# Add common ignores once
if ! grep -q "Ignore large data files" .gitignore; then
cat <<EOL >> .gitignore

# Ignore raw and processed data
*.csv
*.tsv
*.vcf
*.bam
*.sam
*.fasta
*.fa
*.gfasta
*.fastq
*.rmap
*.txt

# Logs and checkpoints
*.out
*.err
*.log
.ipynb_checkpoints/

# System files
.DS_Store
EOL
fi

echo "✅ Directory fully restructured into src/, notebooks/, reports/, configs/, data/raw, data/processed."

