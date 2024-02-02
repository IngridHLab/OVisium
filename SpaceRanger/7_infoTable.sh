#!/usr/bin/env bash
home=$HOME
out=$home/OVisium/Spaceranger/Outs
input=$home/OVisium/Spaceranger/input.csv
current_date_time=$(date +'%Y%m%d')

#' Create a txt file for Rascal solutions of all samples
infoTable=$home/OVisium/SCT_RDS/infoTable_$current_date_time.csv
printf "samples,spotfiles,imgs,json,patient,mutation,variant,age,sample,origin,library,slide,area,run,fastq\n" > $infoTable 

while IFS=, read -r patient mutation variant age sample origin library slide 
area image	run fastq json
do
samples=$out/$library/outs/filtered_feature_bc_matrix.h5
spotfiles=$out/$library/outs/spatial/tissue_positions.csv
imgs=$out/$library/outs/spatial/tissue_hires_image.png
json=$out/$library/outs/spatial/scalefactors_json.json

printf "$samples,$spotfiles,$imgs,$json,$patient,$mutation,$variant,$age,$sample,$origin,$library,$slide,$area,$run,$fastq\n" >> $infoTable 

done < <(tail -n +2 $input)