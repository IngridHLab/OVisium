# OVisium
## Analysis Pipeline for Spatial Transcription Profiling of Fallopian Tube Tissues from Germline BRCA1/2 Mutation Carriers

### OVisium data analysis workflow overview:
<p align ="center">
<img width="500" alt="image" src="https://github.com/NyKepler/OVisium/assets/111468388/7a712847-113a-4790-b8d5-3f6fe2333d29">
</p>

### Analysis of 10X Visium Spatial data

### 1. SpaceRanger and Loupe Browser

`spaceRanger count` default output includes feature-barcode matrices (gene x UMI counts) in tsv.gz or h5 format, a cloup file which can be opened in the *Loupe Browser* and visualize the analysis result ex. graph-based clustering and differential gene expression spatially, and a web_summary.html file to overview the QC matrices of the analysis. Manually annotatioin based on morphology can be exported and integrated to the *Seurat* analysis.

<p align ="center">
<img height="300" alt="lb_2" src="https://github.com/user-attachments/assets/23a10fa6-91ce-4764-9296-876e251c86fe"> <img height="300" alt="lb_3" src="https://github.com/user-attachments/assets/91e2ab63-300a-4b87-9762-55428ba3c1d8"> 
</p>

*SpaceRanger* also provide a pipeline called `spaceranger aggr` to allow user to aggregate multiple capture areas. The command takes a CSV file specifying a list of output files from individual samples and additional categories information for *Loupe Browser* to view or subset.

| library_id |        molecule_h5        |      cloupe_file       |  spatial_folder  |  sample_info   | slide_info |   seq_info   |
|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|
|    PO20    | /path/to/molecule_info.h5 | /path/to/cloupe.cloupe | /path/to/spatial | BRCA_Fimbrial  | V10T06-110 | CTG_2021_075 |
|   PO20_2   | /path/to/molecule_info.h5 | /path/to/cloupe.cloupe | /path/to/spatial | BRCA_Proximal  | V10T06-110 | CTG_2021_075 |
|   PO20a    | /path/to/molecule_info.h5 | /path/to/cloupe.cloupe | /path/to/spatial | BRCA_Otherside | V10S21-050 | CTG_2021_099 |

In the `spaceranger aggr` function, one parameter can be adjusted is either perform a sequencing depth normalization or not before merging, which is recommended by *10XGenomics* to avoid batch effect introduced by sequencing depth. The normalization can be turned off to maximize sensitivity (keep the original amount of reads of each sample) and perform batch correction in the downstream analysis (ex. Harmony). The output is basically same as the output for a single capture area except for the feature-barcode matrix which contains tissue related barcode nucleotide sequence plus number id to distinguish same spots from different tissues/capture areas. For example, `AAACAACGAATAGTTC-1` and `AAACAACGAATAGTTC-2.` The numbering order reflects the sample order in the 'Aggregation CSV' above. 

``` bash
#' Run spaceranger on individual sample
./OVisium/SpaceRanger-2.0.1/SR_count.sh

#' Run spaceranger to aggregate all samples
./OVisium/SpaceRanger-2.0.1/SR_aggregation.sh

#' Generate a samplesheet of the spaceranger outputs for downstream analysis
./OVisium/SpaceRanger-2.0.1/infoTable.sh
```
We also have two samples from two different patients and the capture areas were manually selected in the *Loupe Browser* using the `Visium manual alignment` -> `Fiducial Alignment` and obtained individual jason files and re-run the `spaceRanger count` on fastQ analysis, and lastly generated two separated outputs instead of one.

