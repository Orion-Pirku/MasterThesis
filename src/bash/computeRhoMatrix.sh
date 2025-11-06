#!/bin/bash
set -euo pipefail

ROOT="$HOME/MasterThesis-main/"
cd "$ROOT"
REGIONS_FILES="$ROOT/data/processed/blackcap_genome_features"
computeMatrix reference-point \
  --regionsFileName \
  "$REGIONS_FILES/blackcap1_genes.bed" \
  "$REGIONS_FILES/blackcap_5_UTR.bed" \
  "$REGIONS_FILES/blackcap_3_UTR.bed" \
  "$REGIONS_FILES/blackcap_cds.bed" \
  "$REGIONS_FILES/blackcap_exon.bed" \
  --scoreFileName "$ROOT/data/raw/n100_blackcap_rmaps.bw" \
  --referencePoint center \
  --outFileName "$ROOT/blackcap_rho_features_matrix" \
  --beforeRegionStartLength 20000 \
  --afterRegionStartLength 20000 \
  --binSize 100 \
  --skipZeros \
  --samplesLabel "European Blackcap" \
  -p 10
