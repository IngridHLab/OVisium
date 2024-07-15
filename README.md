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

### 2. Visium Data Analysis with Seurat

Based on the default graph-based clustering results from *Spaceranger* and *Loupe Browser*, we will need to run some batch correctionor sample integration to bring different sequencing libraries as well as the FTEs clusters into alignment. The workflow below is adapted to the *10XGenomics* [batch correction method on Visium data](https://www.10xgenomics.com/resources/analysis-guides/correcting-batch-effects-in-visium-data), [Seurat](https://satijalab.org/seurat/) R package for [analysis of spatial datasets](https://satijalab.org/seurat/articles/spatial_vignette.html), [Harmony](https://www.nature.com/articles/s41592-019-0619-0) and other R
packages.

#### 2.1. Workflow details

The newer version of Seurat introduces the [SCTransform](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) workflow, which is an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow. This procedure skips heuristic steps such as pseudocount addition or log-transformation and improves downstream analysis such as variable gene selection, dimensional reduction and differential expression.

*STutility* is an R-package with the goal of providing an easy-to-use visualization and analysis tool kit for spatial transcriptomics data. The package was developed by Ludvig Larsson and Joseph Bergenstråhle in the Genomics research group [spatialresearch.org](https://www.spatialresearch.org/) at the Royal Institute of Technology (KTH).

```{r}
#' Manually Merge individual spacerange outputs instead of using the `spacerange agg`
#' Use the `infoTable` as mentioned above  
./OVisium/R/Data Preprocessing/Merge_Visium_Samples.R
```
#### 2.2. Imaging and QC
Visualize the QC matrices of individuals in the merged object including "nFeature_RNA", "nCount_RNA", "Mito.percent", "Ribo.percent" and "Hb.percent" by spatial plots and violin plots as well as overall density histogram of the merged object splited by tissue origins (fimbrial, proximal)
Most of the spots have good number of features and expression level.
However, based on the statistic below one could filter away those spots
which have ex.

-   minGenesPerSpot : sets a threshold for the minimum allowed number of
    unique genes in a spot (A). Removing spots which have low number of
    RNA transcripts/features. ex. \>=100

-   minUMICountsPerSpot : sets a threshold for the minimum allowed UMI
    counts in a spot (B). Removing spots which have low RNA expression.
    ex. \>=500

-   minUMICountsPerGene: sets a threshold for the minimum allowed UMI
    counts of a gene across the whole dataset (C). Removing genes low
    expression. ex. \>=100

-   minSpotsPerGene : sets a threshold for the minimum allowed number of
    spots where a gene is detected cross the whole dataset (D). Removing
    low abundant gene. ex. \>=5

-   topN : subset the expression matrix to include only the topN most
    expressed genes.
```{r}
./OVisium/R/Data Preprocessing/Visualization_Data_Quality_Distribution.R

Filter, SCTransform individual samples and merge them through *STutility*. The variable features identified by the SCTransform from individuals (3000) were united and resulted in 6189 variable features:

```{r}
#' According to our data distribution, filter the spots as below before SCTransform:
#' nCount_RNA > 500 &
#' nFeature_RNA > 500 &
#' log10GenesPerUMI > 0.83 &
#' Mito.percent < 15 &
#' Ribo.percent >5 &
#' Ribo.percent < 50 
./OVisium/R/Data Preprocessing/Subset_SCT_Merge.R
```
Alternative, filter, SCTransform individual samples and integrate them through *Seurat*. The 3000 integrated variable features were selected afterward. 
```{r}
#' According to our data distribution, filter the spots as below before SCTransform:
#' nCount_RNA > 500 &
#' nFeature_RNA > 500 &
#' log10GenesPerUMI > 0.83 &
#' Mito.percent < 15 &
#' Ribo.percent >5 &
#' Ribo.percent < 50 
./OVisium/R/Data Preprocessing/Subset_SCT_Integrate_Seurat.R
```



