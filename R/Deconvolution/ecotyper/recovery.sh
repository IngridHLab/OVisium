#!/usr/bin/env bash

#' 4. Manually copy filtered feature bc matrices and spatial tissue position list and high resolution image
#' A bash script can be applied for large number of samples

#' 5. In terminal run ecotyper Visium method for each sample:
#' Edit the yml file:
#' Input Visium directory : "OVisium_data/OVisium_Carcinoma/SampleID"
cd ecotyper
Rscript EcoTyper_recovery_visium.R -c config_recovery_OVisium.yml
