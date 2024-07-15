#!/usr/bin/env bash
#' connect server through ssh Kepler78@130.235.5.221 -p 1220
home=$HOME
path=$home/OVisium/Spaceranger
input=$path/input_rerun.csv
cd $path/Outs

#' Perform sitecheck
export PATH=$path/spaceranger-2.0.1:$PATH
spaceranger sitecheck > sitecheck.txt

#' A sample sheet contains 13 columns 
while IFS=, read -r patient mutation variant age sample origin library slide area image run fastq json
do

if [[ "$json" == "NA" ]]
then
#' Extact count on autodetected area   
#' Jason file is generated automatically based on the provided image   
#' Image should be in either TIFF, QPTIFF or JPEG format
#' The following step is done by the CTG facility after the sequencing usinbg 
#' version 1.2.2               
spaceranger count --id=$library --transcriptome=$path/refdata-gex-GRCh38-2020-A --fastqs=$path/$run/fastq/$fastq --sample=$fastq --image=$path/$run/image/$image --slide=$slide --area=$area --localcores=8 --localmem=500 --no-bam --r1-length=28

else
#' Extract count on manually selected area
#' A jason file is generated in Loupe Browser on the selected area
#' Image should be in either TIFF, QPTIFF or JPEG format  

spaceranger count --id=$library --transcriptome=$path/refdata-gex-GRCh38-2020-A --fastqs=$path/$run/fastq/$fastq --sample=$fastq --image=$path/$run/image/$image --slide=$slide --area=$area --loupe-alignment=$path/$run/image/$json --localcores=8 --localmem=500 --no-bam

fi 

done < <(tail -n +2 $input)
