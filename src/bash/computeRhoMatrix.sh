#!/bin/bash
set -euo pipefail

ROOT="$HOME/MasterThesis-main"
cd "$ROOT"
REGIONS_FILES="$ROOT/data/processed/blackcap_genome_features"
computeMatrix reference-point \
  --regionsFileName \
  "$REGIONS_FILES/blackcap_gene.bed" \
  "$REGIONS_FILES/blackcap_5_UTR.bed" \
  "$REGIONS_FILES/blackcap_3_UTR.bed" \
  "$REGIONS_FILES/blackcap_cds.bed" \
  "$REGIONS_FILES/blackcap_exon.bed" \
  "$REGIONS_FILES/blackcap_intergenic.bed" \
  --scoreFileName "$ROOT/data/raw/n100_blackcap_rmaps.bw" \
  --referencePoint center \
  --outFileName "$ROOT/blackcap_rho_all_features_matrix" \
  --beforeRegionStartLength 20000 \
  --afterRegionStartLength 20000 \
  --binSize 50 \
  --skipZeros \
  --samplesLabel "European Blackcap" \
  -p 10
