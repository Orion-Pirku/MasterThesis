#!/bin/bash

cd $HOME/MasterThesis-main

plotHeatmap \
  --heatmapHeight 12 \
  --heatmapWidth 6 \
  --matrixFile blackcap_rho_features_matrix \
  --outFileName blackcap_rho_features_heatmap.png \
  --plotFileFormat png \
  --plotType se \
  --refPointLabel center \
  --averageType mean \
  --yAxisLabel "Recombination Rate $\rho$" \
  --samplesLabel "European Blackcap" \
  --regionsLabel "Genes" "5' UTR" "3' UTR" "CDS" "EXONS" \
  --colorMap "afmhot" "afmhot" "afmhot" "afmhot" "afmhot" \
  --boxAroundHeatmaps no \
  --perGroup
