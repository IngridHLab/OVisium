#!/usr/bin/env bash

#' 6. In terminal run ecotyper scRNAseq method on individual sample:
#' Use loupe browser to export matrices and cluster annotation files 
Rscript EcoTyper_recovery_visium.R --discovery 'carcinoma' --matrix OVisium_data/OVisium_Carcinoma/OvCa707_2
