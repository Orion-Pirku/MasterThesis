#!/bin/bash
set -euo pipefail

ROOT="$HOME/MasterThesis-main"
cd "$ROOT"
REGIONS_FILES="$ROOT/data/processed/blackcap_genome_features"
computeMatrix reference-point \
  --regionsFileName \
  "$REGIONS_FILES/blackcap_genes.bed" \
  "$REGIONS_FILES/blackcap_5_UTR.bed" \
  "$REGIONS_FILES/blackcap_3_UTR.bed" \
  --scoreFileName "$ROOT/results/blackcap_hotspots/bed/blackcap_hotspots.bw" \
  --referencePoint center \
  --outFileName "$ROOT/blackcap_hotspots_features_matrix" \
  --beforeRegionStartLength 10000 \
  --afterRegionStartLength 10000 \
  --skipZeros \
  --samplesLabel "European Blackcap" \
  -p 10
