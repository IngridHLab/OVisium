#' Spatial Visualization of DEGs that related to TNFa as well as AGR2/3
#' dotplot of markers in individual sample
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"

#' use the harmonized data with all features 
load(paste(rds.dir, "Variable_features_filt_SCT_log2counts+1_harmony.RData",sep = "/"))
data.sub.filt@assays[["SCT"]]@scale.data <- as.matrix(data.sub.filt@assays[["SCT"]]@data)
data.sub.filt <- ScaleData(data.sub.filt)

#' load TNFa markers identify from FTE_DE or Stroma_DE
hallmark_FTE_Stroma_genelist <- readRDS("~/OVisium/DE_functional_analysis/OVisium_SCT_merged/harmony_SCT_res_0.6/Variable_features_filt/MAST/FindMarkers/Combined_filtered_markers_v1/2024-05-04_harmony_sample_log2SCTcounts1_Stroma/GSEA/h.all.v2023.2/h.all.0_FTE&10_Stroma_genelist.rds") 
# hallmark_FTE_Stroma_genelist[["TNFA_SIGNALING (FTE)"]] <- c("TNFAIP2","TNFAIP3","CCL2","CLCF1","CD83","IRF1","ZC3H12A","SAT1", "GADD45A", "BTG2", "PTGS2", "NAMPT", "KLF4", "NFKBIA", "KDM6B", "CCNL1")
hallmark_FTE_Stroma_genelist[["TNFA_SIGNALING (FTE)"]] <- hallmark_FTE_Stroma_genelist[[1]] 
hallmark_FTE_Stroma_genelist[["TNFA_SIGNALING (Stroma)"]] <- hallmark_FTE_Stroma_genelist[[4]] 
hallmark_FTE_Stroma_genelist[["P53_PATHWAY (FTE)"]] <- hallmark_FTE_Stroma_genelist[[3]] 
hallmark_FTE_Stroma_genelist[["P53_PATHWAY (Stroma)"]] <- hallmark_FTE_Stroma_genelist[[7]] 
hallmark_FTE_Stroma_genelist[["New_Markers"]] <- c("AGR2", "AGR3", "TXNIP")
# markers <- hallmark_FTE_Stroma_genelist[c(16:20)] %>% unlist() %>% unique()
marker.list <- hallmark_FTE_Stroma_genelist[c(16:20)]
markers <- marker.list %>% unlist()
un.marker.list <- Map(`[`, marker.list, relist(!duplicated(markers), 
                                                skeleton = marker.list)) 
names(un.marker.list)[1] <- "TNFa (FTE)"
names(un.marker.list)[2] <- "TNFa (Stroma)"
names(un.marker.list)[3] <- "P53 (FTE)"
names(un.marker.list)[4] <- "P53 (Stroma)"
names(un.marker.list)[5] <- "HGSC"
un.markers <- un.marker.list %>% unlist()

## plot individual sample reference gene expression
cscale <- c("darkblue", "blue", "lightblue", "white", 
                      "lightgray", "mistyrose", "red", "darkred", "black")
                      
#' SDC1, MS4A1, KLRC1, TPSB2 and SOX2 doesn't express in all individual sample
Idents(data.sub.filt) <- "sample"
id <- levels(data.sub.filt)
for (i in seq_along(id)) {
  out.dir <- paste(deg.dir, file.name, cluster.ident, 
                   "Spatial_feature_plot/NewMarkers/By_sample", 
                   id[i], sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  data.sample <- SubsetSTData(data.sub.filt, idents = id[i])
  data.sample <- ScaleData(data.sample)
  
  #' plot reference gene spatially
  for (ftr in seq_along(un.markers)) {
    
    #' Feature expression in UWOP plot
    Idents(data.sample) <- "Visium_clusters"
    p1 <- FeaturePlot(data.sample, 
                      features = un.markers[[ftr]], 
                      reduction = "umap_harmony_SCT",
                      order = F, 
                      label = T, 
                      label.size = 3, 
                      label.color = "black", 
                      repel = T
    ) + 
      theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8)) + 
      xlab("UMAP_1") + 
      ylab("UMAP_2") + 
      theme(plot.title = element_blank(), plot.subtitle = element_blank())
    
    #' Violin plot show the expression (log1p of re-corrected counts) in each cluster
    p2 <- VlnPlot(data.sample, 
                  assay = "SCT", 
                  features = un.markers[[ftr]],
                  slot = "data") + 
      theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8)) + 
      xlab("") + 
      ylab("Log2 Normalized Expression") + 
      NoLegend() +
      theme(plot.title = element_blank(), plot.subtitle = element_blank())
    
    #' H&E image
    p3 <- FeatureOverlay(data.sample,
                         features = un.markers[[ftr]], 
                         pt.size = 1,
                         pt.alpha = 0.1,
                         pt.border = F,
                         show.sb = F) +
      theme(plot.title = element_blank(), plot.subtitle = element_blank()) + 
      NoLegend()
    
    #' The scaled data of the harmonized expression count
    p4 <- ST.FeaturePlot(data.sample,
                         features = un.markers[[ftr]], 
                         pt.size = 1,
                         pt.border = F,
                         slot = "scale.data",
                         center.zero = T,
                         cols = cscale,
                         show.sb = F) +
      theme(plot.title = element_blank(), plot.subtitle = element_blank(), 
            strip.text = element_blank()) +
      labs(fill = "ScaleData")
    
    title <- ggdraw() + 
      draw_label(
        paste0("Sample_", id[[i]], ": ", 
               names(un.markers[ftr]), "_", un.markers[[ftr]]),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
    
    p <- plot_grid(p1, p2, p3, p4, ncol = 2, labels = c("A", "B", "C", "D")) 
      
    png(file = paste0(paste("Sample", id[[i]], ftr, un.markers[[ftr]], sep = "_"), 
                      ".png"), width = 4500, height = 4500, res = 500)
    print(    
     plot_grid(title, p, ncol = 1, rel_heights = c(0.08, 1)) 
    )
    dev.off() 
  }
  png(file = paste0(paste("Sample", id[i], "dotplot_by_cluster", sep = "_"), 
                    ".png"), width = 6000, height = 1500, res = 300)
  print(
    DotPlot(data.sample, assay = "SCT", features = un.marker.list,
            group.by = "Visium_clusters", cluster.idents = F, cols="RdYlBu") + 
      theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8), 
            strip.text.x = element_text(angle = 45, size = 11)) + 
      guides(size = guide_legend(title = "Fraction of spots", 
                                 override.aes = list(shape = 21, 
                                                     colour = "grey", 
                                                     fill = "grey"))) +
      xlab("") + ylab("")
  )
  dev.off()
}

## Plot all samples for each gene/marker
Idents(data.sub.filt) <- "Visium_clusters"
for (ftr in seq_along(un.markers)) {
  out.dir <- paste(deg.dir, file.name, cluster.ident, 
                   "Spatial_feature_plot/NewMarkers/By_marker", 
                   paste(ftr, un.markers[[ftr]], sep = "_"), sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  
  #' UMAP feature plots for all samples
  png(file = paste0(paste("AllSamples_UMAPplot", ftr, un.markers[[ftr]], sep = "_"), 
                    ".png"), width = 8000, height = 7000, res = 500)
  print(
    FeaturePlot_scCustom(data.sub.filt, 
                         features = un.markers[[ftr]], 
                         reduction = "umap_harmony_SCT",
                         split.by = "sample",
                         num_columns = 5, 
                         label = T, 
                         label.size = 2, 
                         label.color = "black") & 
      xlab("UMAP_1") & 
      ylab("UMAP_2") 
  )
  dev.off()
  
  #' Vlnplot 
  png(file = paste0(paste("AllSamples_Vlnplot", ftr, un.markers[[ftr]], sep = "_"), 
                    ".png"), width = 9000, height = 2000, res = 500)
  print(
    VlnPlot_scCustom(data.sub.filt, 
                     assay = "SCT", 
                     slot = "data",
                     features = un.markers[[ftr]],
                     group.by = "sample",
                     split.by = "Visium_clusters",
                     pt.size = 0.1) +
      xlab("Sample") +
      ylab("Log2 Normalized Expression") + 
      theme(axis.text.x=element_text(angle = 0))
  )
  dev.off()
  
  png(file = paste0(paste("AllSamples_cVlnplot", ftr, un.markers[[ftr]], sep = "_"), 
                    ".png"), width = 9000, height = 2000, res = 500)
  print(
    VlnPlot_scCustom(data.sub.filt, 
                     assay = "SCT", 
                     slot = "counts",
                     features = un.markers[[ftr]],
                     group.by = "sample",
                     split.by = "Visium_clusters",
                     pt.size = 0.1) +
      xlab("Sample") +
      ylab("Expression") + 
      theme(axis.text.x=element_text(angle = 0))
  )
  dev.off()
  
  #' Spatial feature plot with scaled data
  png(file = paste0(paste("AllSamples_scaleData_spatial", ftr, un.markers[[ftr]], sep = "_"), 
                    ".png"), width = 8500, height = 7000, res = 500)
  print(
    ST.FeaturePlot(data.sub.filt, 
                   features = un.markers[[ftr]], 
                   ncol = 5, 
                   show.sb = F, 
                   pt.size = 0.9,
                   slot = "scale.data",
                   center.zero = T,
                   cols = c("purple", "black", "yellow"), 
                   dark.theme = T) &
      theme(plot.title = element_blank()) &
      theme(text=element_text(size=25)) &
      theme(legend.key.size = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'),
            legend.key.height = unit(2, 'cm')) &
      labs(fill = paste("ScaleData", un.markers[[ftr]], sep = "\n"))
  )
  dev.off()
  
  #' Spatial feature plot with H&E background
  png(file = paste0(paste("AllSamples_log2NormData_spatial", ftr, un.markers[[ftr]], sep = "_"), 
                    ".png"), width = 9000, height = 6000, res = 500)
  print(
    FeatureOverlay(data.sub.filt,
                   features = un.markers[[ftr]],
                   sampleids = 1:18,
                   ncols = 5,
                   pt.size = 1,
                   pt.alpha = 0.5,
                   pt.border = F,
                   show.sb = F, 
                   cols = c("black", "yellow"),
                   dark.theme = T,
                   value.scale = "all") &
      theme(plot.title = element_blank(), plot.subtitle = element_blank()) & 
      labs(fill = paste("Log2Norm", un.markers[[ftr]], sep = "\n")) 
  )
  dev.off()
}  

#' dotplot of all 72 selected markers across 11 OVisium clusters
png(file = "AllSample_dotplot_by_cluster.png", 
    width = 6000, height = 1500, res = 300)
print(
  DotPlot(data.sub.filt, assay = "SCT", features = un.marker.list,
          group.by = "Visium_clusters", cluster.idents = F, cols="RdYlBu") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8), 
          strip.text.x = element_text(angle = 45, size = 11)) + 
    guides(size = guide_legend(title = "Fraction of spots", 
                               override.aes = list(shape = 21, 
                                                   colour = "grey", 
                                                   fill = "grey"))) +
    xlab("") + ylab("")
)
dev.off()

#' dotplot of all 72 selected markers across 18 OVisium samples
png(file = "AllSample_dotplot_by_sample.png", 
    width = 6000, height = 2000, res = 300)
print(
  DotPlot(data.sub.filt, assay = "SCT", features = un.marker.list,
          group.by = "sample", cluster.idents = F, cols="RdYlBu") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8), 
          strip.text.x = element_text(angle = 45, size = 11)) + 
    guides(size = guide_legend(title = "Fraction of spots", 
                               override.aes = list(shape = 21, 
                                                   colour = "grey", 
                                                   fill = "grey"))) +
    xlab("") + ylab("")
)
dev.off()

#' kmean clustered dotplot of all 72 selected markers and default unsupervised clustered (distance = "euclidean" and method = "complete") on 11 OVisium clusters
png(file = "AllSample_clusterd_dotplot_by_marker.png", 
    width = 4000, height = 6000, res = 500)
Clustered_DotPlot(data.sub.filt, assay = "SCT", features = un.markers, k=5) 
dev.off()
