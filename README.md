# OVisium
## Analysis Pipeline for Spatial Transcription Profiling of Fallopian Tube Tissues from Germline BRCA1/2 Mutation Carriers

### OVisium data analysis workflow overview:
<p align ="center">
<img width="800" alt="image" src="https://github.com/NyKepler/OVisium/assets/111468388/7a712847-113a-4790-b8d5-3f6fe2333d29">
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

#### 2.1. Workflow Details

The newer version of Seurat introduces the [SCTransform](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) workflow, which is an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow. This procedure skips heuristic steps such as pseudocount addition or log-transformation and improves downstream analysis such as variable gene selection, dimensional reduction and differential expression.

*STutility* is an R-package with the goal of providing an easy-to-use visualization and analysis tool kit for spatial transcriptomics data. The package was developed by Ludvig Larsson and Joseph Bergenstråhle in the Genomics research group [spatialresearch.org](https://www.spatialresearch.org/) at the Royal Institute of Technology (KTH).

```{r}
#' Manually Merge individual spacerange outputs instead of using the `spacerange agg`
#' Use the `infoTable` as mentioned above to merge the individuals in the same order  
./OVisium/R/Data_Preprocessing/Merge_Visium_Samples.R
```
#### 2.2. Imaging and QC
Visualize the QC matrices of individuals in the merged object including "nFeature_RNA", "nCount_RNA", "Mito.percent", "Ribo.percent" and "Hb.percent" by spatial plots and violin plots as well as overall density histogram of the merged object splited by tissue origins (fimbrial, proximal and other fimbrial). "log10FeaturesPerUMI", the overall complexity of the gene expression (novelty score) by visualizing the number of genes detected per UMI. Additionally, we also checked which genes contribute the most reads by plotting the percentage of counts per gene per spot.
```{r}
#' This script generate QC plots in the supplementary figure 3 & 4
Rscript ./OVisium/R/Data_Preprocessing/Visualization_Data_Quality_Distribution.R
```

#### 2.3. Subset, Filter, SCTransform and Select Variable Features
Filter, SCTransform individual samples and merge them through *STutility*. The variable features will be identified by the SCTransform from the individuals (3000) and then combined:
```{r}
#' According to our data QC distribution, filter the spots as below before SCTransform:
#' nCount_RNA > 500 &
#' nFeature_RNA > 500 &
#' log10GenesPerUMI > 0.83 &
#' Mito.percent < 15 &
#' Ribo.percent >5 &
#' Ribo.percent < 50 
Rscript ./OVisium/R/Data_Preprocessing/Subset_SCT_Merge.R

#' Output file: OVisium_SCT_merged.rds

#' Get statistic summary of the merged data
Rscript ./OVisium/R/Help_functions/Statistic_Result.R
```

Alternative, filter, SCTransform individual samples and integrate or merge them through *Seurat*. Gene symbols will be first converted to Ensembl ids. Individuals will be integrated and 3000 variable features will be selected after integration. Individuals can be also merged and individual variable features will be combined instead. The output files can be used in the deconvolution later.     
```{r}
#' According to our data distribution, filter the spots as below before SCTransform:
#' nCount_RNA > 500 &
#' nFeature_RNA > 500 &
#' log10GenesPerUMI > 0.83 &
#' Mito.percent < 15 &
#' Ribo.percent >5 &
#' Ribo.percent < 50 
Rscript ./OVisium/R/Data_Preprocessing/Subset_SCT_Integrate_Seurat.R

#' Output file 1: OVisium_ens_SCT_integrated.rds and
#' Output file 2: OVisium_ens_SCT_merged.rds
```

#### 2.4. PCA Dimensional Reduction and Harmony Integration (Standard *Seurat*)

Using the standard seurat analysis pipeline to process the data. The PCA is performed on the scale.data from the variable features. The scale.data of the SCTransform data is obtained through `GetResidual` function. The optimal number of principle components is determined by  quantitative elbowplot method. The PCA step is applied in the previous step 2.3 right after merging. 

The PCs threshold in the quantitative elbowplot:

-  Cutoff_1: The point where the PCs contribute less than 5% of standard deviation and the PCs cumulatively contribute over 90% of the standard deviation. 

-  Cutoff_2: The point where the percent change in variation between the consecutive PCs is less than 0.1%. 

<p align ="center">
<img width="400" alt="image" src="https://github.com/user-attachments/assets/7458b577-7437-4c4a-8abb-9002b82d0345">
</p>

```{r}
#' This script is sourced in the step 2.3
#' Generate supplementary figure 8A
Rscript ./OVisium/R/Data_Preprocessing/Choose_Num_PCs.R
```
Next, *[Harmony](https://www.nature.com/articles/s41592-019-0619-0)* R package will be used to integrates individual Visium samples (sample ids as covarient) by adjusting PCA embedding value of the features accordingly instead of the actual expression level. Further dimensional reduction ex. TSNE, UMAP and [UWOT](https://github.com/satijalab/seurat/issues/2025) will be performed on both the PCA and the harmony embedding values so the data can be visualized on 2-dimensional plots. More details of how to use *Harmony* in *Seurat* can be found in this [vignette](https://cran.r-project.org/web/packages/harmony/vignettes/Seurat.html).

<p align ="center">
<img width="1000" alt="image" src="https://github.com/user-attachments/assets/f29d4d53-e9a0-46fd-85d6-52245189f793">
</p>

```{r}
#' Covarient: group.by.vars = "sample"
#' In the covergence plot, Harmony objective function value should descend over each harmony_idx
#' Generate plots in the supplementary figure 8B-D
Rscript ./OVisium/R/Data_Preprocessing/Harmony_Integration_DimReduction.R

#' Output file: OVisium_SCT_merged_dims.rds
```

#### 2.5. Determine the Number of Clusters
The optimal number of clusters for the Visium data can be determined by the *clustree* R package, which utilizes [C7 (Gene730)](https://www.ncbi.nlm.nih.gov/gene/730) house keeping gene and [SC3 stability](https://pubmed.ncbi.nlm.nih.gov/28346451/) to visualize the clusterings, together with the morphological appearance of the tissue. The resolution of clustering is tested between 0.1 and 1 with 0.1 increasement in the *Seurat* `FindClusters` function.

```{r}
#' Perform clustering on PCA and Harmony data
#' clustree plots for the supplementary figure 9
Rscript ./OVisium/R/Data_Preprocessing/Choose_Number_Clusters.R

#' Output file: OVisium_SCT_merged_clust.rds

#' Export clusters annotation for Loupe Browser
Rscript ./OVisium/R/Help_functions/Export_Clusters_LoupeBrowser.R

#' Subset and rename clusters based on the morphology
Rscript ./OVisium/R/Help_functions/Subset_Rename_Clusters.R

#' Plot clusters annotation spatially
#' The optimal resolution based on morphology and preliminary DEG is 0.6 for the OVisium data 
Rscript ./OVisium/R/Help_functions/Plot_Clusters_Spatial.R
```

#### 2.6 Differential Expression Analaysis
Before running the DEA using the standard *Seurat* `FindMarkers` function, we will perform log2 transformation for better fold-change calculation and harmony batch correction on the expression matrix (SCTransformed count data) to remove batch effects from the Visium preparation. 

The detail steps for the OVisium data:
1. Use all spots from 12 clusters in the OVisium_SCT_merged_clust.rds.
2. Subset SCTcounts from 6189 most variable genes united from 18 samples.
3. Calculate log2(SCTcounts+1) on spot x gene matrix.
4. Perform HarmonyMatrix using sample id as covariant.

$$
HarmonyLog2Data = HarmonyMatrix(log2(SCT\ Counts_{variable\ features} + 1))
$$

6. Subset away spots from cluster #4
7. Filter variable genes based on their distribution.
     -  Total counts per gene >500​ &
     -  Max counts > 4 & 
     -  Percentage of spots > 1% &​
     -  98 Percentile of log2(Data+1) > log2(2.5) &
     -  Remove genes contain: "^MT-", "^RP[SL]",  "^MTRNR", "^LINC"​

<p align ="center">
<img width="500" alt="image" src="https://github.com/user-attachments/assets/b02ddba9-82b2-4a03-8a6a-7217baf80513">
</p>

```{r}
#' Input file: OVisium_SCT_merged_clust.rds
Rscript ./OVisium/R/Differential_Expression_Analysis/HarmonyMatrix_log2DataPlus1_QC_Filter_VariableGenes.R

#' Output file: Variable_features_filt_SCT_log2counts+1_harmony.RData
#' Additional output files for the manuscript:
#' 11clusters_filt_vfeatureas_SCTcounts_log2+1_harmony_samples.csv
#' Filtered_spots_sample_cluster_annotation.csv
```

#### 2.6.1 `FindMarkers` between Clusters using MAST 
As a default, *Seurat* performs differential expression based on the non-parameteric Wilcoxon rank sum test. Here, we choose “MAST”, GLM-framework that treates cellular detection rate as a covariate ([Finak et al, Genome Biology, 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5)).

The following steps are used for the OVisium data:
1. Scale and center *HarmonyLog2Data*.​
2. Average difference of *HarmonyLog2Data* between the two groups (`ident.1`, `ident.2`, 11x10).
3. Calculate the `rank`: -log10(p_val_adj+2.225074e-308)*avg_log2FC​
4. Filter with `min.diff.pct` > 0.1 or 0.25, abs(`avg_log2FC`) > log2(1.5 or 2.5) and abs(`rank`) > -log10(10E-5)* log2(1.5 or 2.5). Depends on `ident.1` & `ident.2` are from the same tissue compartment or not.
5. Combine the results with the same `ident.1`.​
6. Order they by `avg_log2FC` descending.​
7. Remove duplicates with smaller `avg_log2FC`.​
8. Top20 of each: group by `ident.2 and slice max by `rank`.​
9. Merge all results from different `ident.1`.​
10. Order they by `avg_log2FC` descending.​
11. Remove duplicates with smaller `avg_log2FC`.
12. Top20 of all: Filter `avg_log2FC` > log2(2), then group by `ident.1` and slice max by `rank`. 

```{r}
#' DEA between 11 clusters
#' Generate complexheatmap Figure 3A in the manuscript
Rscript ./OVisium/R/Differential_Expression_Analysis/FindMarkers_ComplexHeatmap.R

#' Output folder name: 2024-04-20_harmony_sample_log2SCTcounts1
```
The results data frame has the following columns :
- `p_val` : p_val (unadjusted)
- `avg_log2FC` : log2 fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group.
- `pct.1` : The percentage of cells where the feature is detected in the first group `ident.1`
- `pct.2` : The percentage of cells where the feature is detected in the second group `ident.2`
- `p_val_adj` : Adjusted p-value, based on bonferroni correction using all features in the dataset.

  
#### 2.6.2 `FindMarkers` between Clusters within the Same Tissue Compartment
In the OVisium data, we have 2 clusters from the epithelial compartment and 8 clusters from the stroma compartment. One mix cluster located on conjunction of epithelial and stroma is excluded in this analysis.

##### FTE Clusters:
1. Subset FTE clusters (0_FTE and 1_FTE)
2. Scale and center *HarmonyLog2Data*.​
3. Average difference of *HarmonyLog2Data* between the two groups (`ident.1`, `ident.2`, 2x1).
4. Calculate the `rank`: -log10(p_val_adj+2.225074e-308)*avg_log2FC​
5. Filter with `min.diff.pct` > 0.1, abs(`avg_log2FC`) > 0.25 and abs(`rank`) > -log10(10E-5)* 0.25.
6. Combine the results with the same `ident.1`.​
7. Order they by `avg_log2FC` descending.​
8. Remove duplicates with smaller `avg_log2FC`.​
9. Top20 of each: group by `ident.2` and slice max by `rank`.​
10. Merge all results from different `ident.1`.​
11. Order they by `avg_log2FC` descending.​
12. Remove duplicates with smaller `avg_log2FC`.
13. Top20 of all: Filter `avg_log2FC` > log2(1.5), then group by `ident.1` and slice max by `rank`. 

```{r}
#' DEA between 2 clusters
#' Generate complexheatmap
Rscript ./OVisium/R/Differential_Expression_Analysis/FindMarkers_FTE_ComplexHeatmap.R

#' Output fold name: 2024-05-16_harmony_sample_log2SCTcounts1_FTE_only  
```
##### Stroma Clusters:
1. Subset Stroma clusters (2, 5-11_Stroma)
2. Scale and center *HarmonyLog2Data*.​
3. Average difference of *HarmonyLog2Data* between the two groups (`ident.1`, `ident.2`, 8x7).
4. Calculate the `rank`: -log10(p_val_adj+2.225074e-308)*avg_log2FC​
5. Filter with `min.diff.pct` > 0.1, abs(`avg_log2FC`) > log2(1.5) and abs(`rank`) > -log10(10E-5)* log2(1.5).
6. Combine the results with the same `ident.1`.​
7. Order they by `avg_log2FC` descending.​
8. Remove duplicates with smaller `avg_log2FC`.​
9. Top20 of each: group by `ident.2` and slice max by `rank`.​
10. Merge all results from different `ident.1`.​
11. Order they by `avg_log2FC` descending.​
12. Remove duplicates with smaller `avg_log2FC`.
13. Top20 of all: Filter `avg_log2FC` > log2(2), then group by `ident.1` and slice max by `rank`. 

```{r}
#' DEA between 8 clusters
#' Generate complexheatmap
Rscript ./OVisium/R/Differential_Expression_Analysis/FindMarkers_Stroma_ComplexHeatmap.R

#' Output fold name: 2024-05-04_harmony_sample_log2SCTcounts1_Stroma  
```
