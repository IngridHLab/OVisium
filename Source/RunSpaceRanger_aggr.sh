#!/usr/bin/env bash

cd /home/minerva/VisiumST/Data/Spaceranger

export PATH=/home/minerva/Spaceranger/spaceranger-2.0.1:$PATH

## Perform sitecheck
# spaceranger sitecheck > sitecheck.txt

## Perform aggregation of multiple slices. Specifiy 'none' if no seq depth normalization
spaceranger aggr --id=AGG_BRCA_20_20230511 --csv=agg_BRCA_18.csv --normalize=none 


