#!/usr/bin/env bash

home=$HOME
path=$home/OVisium/Spaceranger
input=$path/input.csv
cd $path/Outs
current_date_time=$(date +'%Y%m%d')

#' Perform sitecheck
export PATH=$path/spaceranger-2.0.1:$PATH
spaceranger sitecheck > sitecheck.txt

#' Create a txt file for all samples
aggTable=$path/aggTable_$current_date_time.csv
printf "library_id,molecule_h5,cloupe_file,spatial_folder,patient,mutation,variant,age,sample,origin,slide,area,run,fastq\n" > $aggTable 

while IFS=, read -r patient mutation variant age sample origin library slide area image	run fastq json
do

molecule_h5=$path/Outs/$library/outs/molecule_info.h5
cloupe_file=$path/Outs/$library/outs/cloupe.cloupe
spatial_folder=$path/Outs/$library/outs/spatial

printf "$library,$molecule_h5,$cloupe_file,$spatial_folder,$patient,$mutation,$variant,$age,$sample,$origin,$slide,$area,$run,$fastq\n" >> $aggTable 

done < <(tail -n +2 $input)

#' Perform aggregation of multiple slices
#' A sample sheet should be prepared accordingly as the example csv file
#' Specifiy 'none' if no seq depth normalization during the aggregation or 'mapped' to have seq depth normalization
spaceranger aggr --id=AGG_BRCA_18_$current_date_time --csv=$aggTable --normalize=none


