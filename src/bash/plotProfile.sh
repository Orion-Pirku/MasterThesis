#!/bin/bash

cd $HOME/MasterThesis-main

plotProfile \
  --plotHeight 8 \
  --plotWidth 12 \
  --matrixFile blackcap_rho_features_matrix \
  --outFileName blackcap_rho_features_profile.png \
  --plotFileFormat png \
  --plotType se \
  --refPointLabel center \
  --averageType mean \
  --yAxisLabel "Recombination Rate $\rho$" \
  --samplesLabel "European Blackcap" \
  --regionsLabel "Genes" "5' UTR" "3' UTR" \
  --dpi 200
