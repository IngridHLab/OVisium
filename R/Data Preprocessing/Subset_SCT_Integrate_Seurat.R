#' Merge sample using standard seurat functions
#' This workflow is for preparing query in deconvolution
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/10_choose_pca_dims.R", sep = "/"), 
       print.eval = TRUE, local = TRUE)

#' Extract features from aggregated Visium matrix
features <- 
  readr::read_tsv(
    paste(sr.dir, "Outs/AGG_BRCA_18_20231018/outs/filtered_feature_bc_matrix/features.tsv.gz", 
          sep = "/"), col_names = FALSE) %>% 
  dplyr::rename("ENSEMBL" = 1, "SYMBOL" = 2, "TYPE" = 3)

mito.features <- features %>%
  dplyr::filter((grepl("^MT-", SYMBOL)))
ribo.features <- features %>%
  dplyr::filter((grepl("^RP[SL]", SYMBOL)))

setwd(rds.dir)
sct.data.list <- list()
inputTable <- readr::read_csv("inputTable_20230804.csv")

for (i in seq_along(inputTable$sample)) { 
  
  set.seed(1220)
  
  ## Create Seurat object with spatial image 
  image <- Read10X_Image(image.dir = inputTable$imgs[[i]], 
                         image.name = "tissue_lowres_image.png",
                         filter.matrix = TRUE)
  
  data <- Load10X_Spatial(data.dir = inputTable$samples[[i]], 
                          filename = "filtered_feature_bc_matrix.h5", 
                          assay = "Spatial", 
                          filter.matrix = TRUE, 
                          image = image, 
                          use.names = FALSE)
  
  #' Get mitochondrial percentage and gene/UMI 
  data[["Mito.percent"]] <- 
    PercentageFeatureSet(data, features = mito.features$ENSEMBL)
  data[["Ribo.percent"]] <- 
    PercentageFeatureSet(data, features = ribo.features$ENSEMBL)
  data$log10GenesPerUMI <- 
    log10(data$nFeature_Spatial) / log10(data$nCount_Spatial)
  data <- subset(data, subset = Mito.percent < 15 & 
                                Ribo.percent > 5 &
                                Ribo.percent < 50 &
                            nFeature_Spatial > 500 &
                              nCount_Spatial > 500 & 
                            log10GenesPerUMI > 0.83)

  #' SCTransform and regress out difference caused by technical batch.
  data <- SCTransform(data, assay = "Spatial", verbose = TRUE, 
                      return.only.var.genes = T)
  
  #' Add sample order as prefix to the barcodes and remove extension/suffix "-1"
  data <- RenameCells(data, new.names = paste(colnames(data), i, sep = "_"))
  
  #' Add metadata
  metadata <- inputTable[i, 5:15][rep(1, ncol(data)), ] %>% as.data.frame() 
  row.names(metadata) <- Cells(data)
  data <- AddMetaData(data, metadata = metadata)

  #' Save under list for later-on merging    
  sct.data.list[[i]] <- data
}
saveRDS(sct.data.list, file = "OVisium_ens_SCT_data_list.rds")


#' perform integration:
set.seed(1220)
ifeatures <- SelectIntegrationFeatures(sct.data.list, nfeatures = 3000) 
sct.data.list <- 
  PrepSCTIntegration(sct.data.list, anchor.features = ifeatures)
anchors <- 
  FindIntegrationAnchors(sct.data.list, normalization.method = "SCT",
                         anchor.features = ifeatures) 
merged <- IntegrateData(anchors, normalization.method = "SCT")
saveRDS(anchors, file = paste(rds.dir, "Query_ens_SCT_integrated_anchors.rds", 
                                  sep = "/"))


#' npcs is decided by 'choose.pcs' function 
threshold.integ <- choose.pcs(data=merged, assay = "integrated")
file.rename("pca_elbow.png", "pca_elbow_Ovisium_ens_SCT_integrated.png")
file.rename("pca_quant_elbow.png", 
            "pca_quant_elbow_Ovisium_ens_SCT_integrated.png")
file.rename("pca_heatmap.png", "pca_heatmap_Ovisium_ens_SCT_integrated.png")

set.seed(1220)
merged <- RunPCA(merged, assay = "integrated", npcs = threshold.integ, 
                 reduction.name = "pca_integrated")

#' save data
saveRDS(merged, file = paste(rds.dir, "OVisium_ens_SCT_integrated.rds", 
                             sep = "/"))


#' merge all 18 SCT data
merged <- merge(x = sct.data.list[[1]],
                      y = sct.data.list[2:length(sct.data.list)],
                      project = "Visium_BRCA", merge.data = TRUE)
vfeatures <- rownames(merged[["SCT"]]@scale.data)

#' set variable features and residues for the SCT assay. 
DefaultAssay(merged) <- "SCT"
VariableFeatures(merged) <- vfeatures
merged <- GetResidual(merged, 
                      assay = "SCT", 
                      features = rownames(merged[["SCT"]]@data),
                      umi.assay = "Spatial",
                      replace.value = TRUE)

#' npcs is decided by 'choose.pcs' function 
npcs.sct <- choose.pcs(data=merged, assay = "SCT")
file.rename("pca_elbow.png", "pca_elbow_Ovisium_ens_SCT_merged.png")
file.rename("pca_quant_elbow.png", 
            "pca_quant_elbow_Ovisium_ens_SCT_merged.png")
file.rename("pca_heatmap.png", "pca_heatmap_Ovisium_ens_SCT_merged.png")

set.seed(1220)
merged <- RunPCA(merged, assay = "SCT", npcs = npcs.sct, 
                 reduction.name = "pca_SCT")
saveRDS(merged, file = paste(rds.dir, "OVisium_ens_SCT_merged.rds", 
                             sep = "/"))
