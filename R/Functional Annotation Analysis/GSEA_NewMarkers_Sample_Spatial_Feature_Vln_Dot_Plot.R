#' Individual sample
#' Visualization of individual DEGs that related to TNFa and AGR2/3
#' UMAP Feature plot, Violin plot split by cluster, H&E image and spatial feature plot of individual sample and feature
#' dotplot of all features split by clusters

library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"

#' use the harmonized data with all features 
load(paste(rds.dir, "Variable_features_filt_SCT_log2counts+1_harmony.RData",sep = "/"))

#' load TNFa markers identify from FTE_DE or Stroma_DE
hallmark_FTE_Stroma_genelist <- readRDS("~/OVisium/DE_functional_analysis/OVisium_SCT_merged/harmony_SCT_res_0.6/Variable_features_filt/MAST/FindMarkers/Combined_filtered_markers_v1/2024-05-04_harmony_sample_log2SCTcounts1_Stroma/GSEA/h.all.v2023.2/h.all.0_FTE&10_Stroma_genelist.rds") 
hallmark_FTE_Stroma_genelist[["TNFA_SIGNALING (FTE)"]] <- 
  c("AGR2", "AGR3", "TNFAIP2","TNFAIP3","CCL2","CLCF1", "CD83","IRF1","ZC3H12A",
    "SAT1", "GADD45A", "BTG2", "PTGS2", "NAMPT", "KLF4", "NFKBIA", "KDM6B", "CCNL1")
markers <- hallmark_FTE_Stroma_genelist[c(16,4)] %>% unlist() %>% unique()


#' plot individual sample reference gene expression
cscale <- c("darkblue", "blue", "lightblue", "white", 
                      "lightgray", "mistyrose", "red", "darkred", "black")
                      
#' SDC1, MS4A1, KLRC1, TPSB2 and SOX2 doesn't express in all individual sample
Idents(data.sub.filt) <- "sample"
id <- levels(data.sub.filt)
for (i in seq_along(id)) {
  out.dir <- paste(deg.dir, file.name, cluster.ident, id[i], sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  
  #' plot reference gene spatially
  for (ftr in seq_along(markers)) {
    
    #' Feature expression in UWOP plot
    cells <- WhichCells(data.sub.filt, idents = id[i])
    p1 <- FeaturePlot(data.sub.filt, 
                      features = markers[ftr], 
                      reduction = "uwot_harmony_SCT",
                      cells = cells, 
                      slot = "data", 
                      order = T, 
                      label = T, 
                      label.size = 2, 
                      label.color = "black", 
                      repel = T
    ) + 
      theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8)) + 
      xlab("Harmony:UWOT_1") + 
      ylab("Harmony:UWOT_2") +
      labs(color = "Data")
    
    #' Violin plot show the expression in each cluster
    p2 <- VlnPlot(data.sub.filt, 
                  split.by = "Visium_clusters",
                  assay = "SCT", 
                  idents = id[i], 
                  features = markers[ftr]) + 
      theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8)) + 
      xlab("") + ylab("")
    
    #' the data slot is log1p of re-corrected counts after Harmony
    p3 <- FeatureOverlay(data.sub.filt,
                         features = markers[ftr], 
                         sampleids = i,
                         pt.size = 1,
                         pt.alpha = 0.1,
                         pt.border = F,
                         show.sb = F) +
      theme(plot.title = element_blank(), plot.subtitle = element_blank()) +
      labs(fill = "Data")
    
    #' the scale.data is the residual from the SCT models
    p4 <- ST.FeaturePlot(data.sub.filt,
                         features = markers[ftr], 
                         indices = i,
                         pt.size = 1,
                         pt.border = F,
                         slot = "scale",
                         center.zero = T,
                         cols = cscale,
                         show.sb = F) +
      theme(plot.title = element_blank(), plot.subtitle = element_blank()) +
      labs(fill = "ScaleData")
    
    png(file = paste0(paste("deg", id[i], ftr, markers[ftr], sep = "_"), 
                      ".png"), width = 4500, height = 4000, res = 500)
    print(    
      plot_grid(p1, p2, p3, p4, ncol = 2, labels = c("A", "B", "C", "D"))
    )
    dev.off() 
    
    png(file = paste0(paste("dotplot", id[i], "ref_by_cluster", sep = "_"), 
                      ".png"), width = 5000, height = 1200, res = 300)
    print(
      DotPlot(data.sub.filt, assay = "SCT", features = markers, idents = id[i], 
              group.by = "Visium_clusters", cluster.idents = F, cols="RdYlBu") + 
        theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8)) + 
        xlab("") + ylab("")
    )
    dev.off()
  }
  
}
