# ecological-niche-modeling
# Ecological Niche & Spatial Invasion Risk Modeling of *Aedes albopictus* in Turkey

This repository contains the end-to-end computational pipeline developed for my undergraduate graduation thesis in Molecular Biology and Genetics at Istanbul Technical University (ITU).

## Overview
The study models the fundamental and realized ecological niche of the invasive arboviral vector *Aedes albopictus* across Turkey and the Mediterranean basin under current and future climate projections (CMIP6).

## Methodology
- **Algorithm:** Maximum Entropy modeling (`maxnet` in R).
- **Validation:** 5-Fold Cross-Validation, Permutation Importance, and ROC-AUC evaluation.
- **Environmental Covariates:** Bioclimatic variables from WorldClim 2.1 (BIO11, BIO12, BIO14, BIO15).
- **Anthropogenic Variables:** OpenStreetMap road network proximity (15-km sigmoidal decay buffer) and gridded population density (WorldPop).
- **Future Projections:** CMIP6 global climate models (MPI-ESM1-2-HR) under SSP2-4.5 and SSP5-8.5 scenarios for mid-century (2041–2060) and late-century (2061–2080).

## Requirements
- R >= 4.0
- Libraries: `terra`, `sf`, `maxnet`, `ggplot2`, `tidyterra`, `geodata`, `pROC`, `corrplot`
