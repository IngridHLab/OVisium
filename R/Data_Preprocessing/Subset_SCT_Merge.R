#' Subset individual sample from the merged data file
#' SCTransform on individual samples 
#' Perform spatial correlation analysis 
#' Import morphology annotation 
#' Merge all samples again for downstream analysis
#' https://github.com/satijalab/seurat/issues/5183
#' https://github.com/satijalab/seurat/issues/6185

library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Data_Preprocessing/Choose_Num_PCs.R", sep = "/"),
       print.eval = TRUE, local = TRUE)
       
merged <- readRDS(paste(rds.dir, "OVisium_merged_qc.rds", sep = "/"))
dir.create(paste(work.dir, "Spatgenes", "spat_corr_list", sep = "/"), 
           recursive = T)
setwd(rds.dir)
Idents(merged) <- "sample"
sct.data.list <- list()
spatgene.list <- list()
VariableFeaturePlot.list <- list()

for (i in levels(Idents(merged))) { 
  
  set.seed(1220)
  #' subset each library 
  data <- SubsetSTData(merged, idents = i)
  id <- levels(Idents(data))
  
  
  #' According to our data distribution
  #' nCount_RNA > 500
  #' nFeature_RNA > 500
  #' log10GenesPerUMI > 0.83
  #' Mito.percent < 15
  #' Ribo.percent < 50 to remove extreme outliner from sample #4
  data <- SubsetSTData(data, expression = Mito.percent < 15 &
                                          Ribo.percent > 5 &
                                          Ribo.percent < 50 &
                                          nFeature_RNA > 500 & 
                                            nCount_RNA > 500 & 
                                      log10FeaturesPerUMI > 0.83)
  
  #' SCTransform and regress out difference caused by technical batch.
  data <- SCTransform(data, verbose = TRUE, return.only.var.genes = T)

  #' run spatial correlation analysis
  #' https://rdrr.io/github/jbergenstrahle/STUtility/man/CorSpatialGenes.html
  spatgenes <- CorSpatialGenes(data, assay = "SCT", slot = "scale.data")
  write.table(spatgenes, 
              file=paste0(work.dir, "/Spatgenes/spat_corr_list/", 
                          paste("Sample_", id, ".csv", sep ="")), 
              sep =",", quote = F, row.names = F, col.names = T)
  
  #' Save under list for later-on merging    
  sct.data.list[[i]] <- data
  spatgene.list[[i]] <- spatgenes
  
  #' plot variable features
  vfeatures <- data %>% 
    VariableFeatures() %>% 
    as.data.frame() %>% 
    dplyr::rename("Feature" = 1)
  p <- VariableFeaturePlot(data) %>% 
    LabelPoints(plot = . , points = vfeatures$Feature[1:20], repel = TRUE) +
    labs(title = paste0("Sample_", i),
         subtitle = data@meta.data[["origin"]][1],
         caption  = "") +
    theme(plot.caption = element_text(hjust = 0, face= "italic"), 
          plot.title.position = "plot",
          plot.caption.position =  "plot",
          legend.position="bottom") 
  VariableFeaturePlot.list[[i]] <- p
}
saveRDS(sct.data.list, file = "OVisium_SCT_data_list.rds")

#' Obtain all variable features after merging
merged <- MergeSTData(x = sct.data.list[[1]],
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
file.rename("pca_elbow.png", "pca_elbow_Ovisium_SCT_merged.png")
file.rename("pca_quant_elbow.png", 
            "pca_quant_elbow_Ovisium_SCT_merged.png")
file.rename("pca_heatmap.png", "pca_heatmap_Ovisium_SCT_merged.png")

set.seed(1220)
merged <- RunPCA(merged, assay = "SCT", npcs = npcs.sct, 
                 reduction.name = "pca_SCT")
saveRDS(merged, file = paste(rds.dir, "OVisium_SCT_merged.rds", 
                             sep = "/"))

#' variable feature plots
setwd(qc.dir)
png(file = "VariableFeaturePlot_individual.png", 
    width = 7000, height = 9000, res = 350) 
print(gridExtra::grid.arrange(grobs=VariableFeaturePlot.list, ncol=4, nrow=5)
)
dev.off()
saveRDS(VariableFeaturePlot.list, file = "VariableFeaturePlot_list.rds")
