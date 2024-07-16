#' Visualization of differential gene expression on reference cell markers
library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

#' Reference cell markers from single cell data
ref.markers<-read.csv(paste(deg.dir, "single_cell_markers.csv", sep = "/"))
sc.markers <- ref.markers$Hammound_2022 %>% data.frame() %>% dplyr::rename(., "gene"=1)
sc.markers.filt <- sc.markers$gene[c(1:25,27,30:34,36:40,42:45)]
features = list("Epithelial" = sc.markers$gene[1:6], 
                "Stromal" = sc.markers$gene[7:25],
                "Immune" = sc.markers$gene[c(27,30:34,36:40)], #27,30:34,36:40
                "OC" = sc.markers$gene[c(42:45)]) #41,43:45

file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"
#' use the harmonized data with all features 
load(paste(rds.dir, "Variable_features_filt_SCT_log2counts+1_harmony.RData",sep = "/"))
#data <- readRDS(paste(rds.dir, paste0(file.name, "_11_clust_harmony_data.rds"), sep = "/"))

#' Use all the most variable genes as some of the reference genes are not 
#' full filled the filtration.
data <- data.sub.filt
data@assays[["SCT"]]@data <- log1p(data@assays[["SCT"]]@counts)
data <- ScaleData(data, features = row.names(data))
#' reordering only for better dotplot but not for vlnplot or featureplot
Idents(data) <- factor(Idents(data),
                       levels = c("11_Stroma", "10_Stroma", "9_Stroma", 
                                  "8_Stroma", "7_Stroma", "6_Stroma", 
                                  "5_Stroma","2_Stroma","3_Mix","1_FTE", 
                                  "0_FTE"))
data[["Visium_clusters_rev"]] <- Idents(data)

#' Dotplot of reference markers of merged data
out.dir <- paste(deg.dir, file.name, cluster.ident, sep = "/")
dir.create(out.dir, recursive = T)
setwd(out.dir)
png(file = paste("dotplot_ref_by_cluster_origin.png", sep = "."), 
    width = 5000, height = 3000, res = 350)
DotPlot(data, assay = "SCT", 
        features = features, 
        group.by = "Visium_clusters_rev", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 7,
        split.by = "Tissue_origin") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(angle = 0, size = 14)) +
  guides(size=guide_legend(title="Fraction of spots", 
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("") 
dev.off()

png(file = paste("dotplot_ref_by_cluster_ny.png", sep = "."), 
    width = 5000, height = 2000, res = 400)
DotPlot(data, assay = "SCT", 
        features = features, 
        group.by = "Visium_clusters_rev", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 7) + 
  theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8),
        strip.text.x = element_text(angle = 0, size = 14),
        axis.text.y = element_text(angle = 0, size = 14, hjust = 0.5, vjust = 0.5)) +
  guides(size=guide_legend(title="Fraction of spots", 
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("") +
  theme(legend.position = "bottom")
dev.off()

png(file = paste("dotplot_ref_by_sample_ny.png", sep = "."), 
    width = 5000, height = 2400, res = 350)
DotPlot(data, assay = "SCT", 
        features = features, 
        group.by = "sample", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 7,
        split.by = "Tissue_origin") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(angle = 0, size = 14)) +
  guides(size=guide_legend(title="Fraction of spots",
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("") 
dev.off()


png(file = paste("dotplot_ref_by_tissue.png", sep = "."), 
    width = 5000, height = 1200, res = 300)
DotPlot(data, assay = "SCT", 
        features = features, 
        group.by = "Tissue_origin", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 7) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(angle = 0, size = 14)) +
  guides(size=guide_legend(title="Fraction of spots",
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("") 
dev.off()

png(file = paste("dotplot_ref_by_cluster_mutation.png", sep = "."), 
    width = 5000, height = 3000, res = 350)
DotPlot(data, assay = "SCT", 
        features = features, 
        group.by = "Visium_clusters_rev", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 7,
        split.by = "mutation") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(angle = 0, size = 14)) +
  guides(size=guide_legend(title="Fraction of spots", 
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("")
dev.off()

png(file = paste("dotplot_ref_by_variant.png", sep = "."), 
    width = 5000, height = 1500, res = 350)
DotPlot(data, assay = "SCT", 
        features = features, 
        group.by = "variant", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 7) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(angle = 25, size = 14)) +
  guides(size=guide_legend(title="Fraction of spots",
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("")
dev.off()

#' Vlnplot of reference markers of merged data
png(file = paste("vlnplot_ref_by_cluster_ny.png", sep = "."), 
    width = 8200, height = 5000, res = 300)
VlnPlot(data, assay = "SCT", group.by = "Visium_clusters", 
        features = sc.markers$gene, ncol = 8, pt.size = 0) + 
  theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
  xlab("") + 
  ylab("")
dev.off()

png(file = paste("vlnplot_ref_by_sample_ny.png", sep = "."), 
    width = 12000, height = 5000, res = 300)
VlnPlot(data, assay = "SCT", group.by = "sample", 
        features = sc.markers$gene, ncol = 8, pt.size = 0) + 
  theme(axis.text.x=element_text(angle = 0)) + 
  xlab("") + 
  ylab("")
dev.off()

png(file = paste("vlnplot_ref_by_tissue.png", sep = "."), 
    width = 8000, height = 5000, res = 300)
VlnPlot(data, assay = "SCT", group.by = "Tissue_origin", 
        features = sc.markers$gene, ncol = 8, pt.size = 0) + 
  theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
  xlab("") + 
  ylab("")
dev.off()

png(file = paste("vlnplot_ref_by_mutation.png", sep = "."), 
    width = 8000, height = 5000, res = 300)
VlnPlot(data, assay = "SCT", group.by = "mutation", 
        features = sc.markers$gene, ncol = 8, pt.size = 0) + 
  theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
  xlab("") + 
  ylab("")
dev.off()

png(file = paste("vlnplot_ref_by_variant.png", sep = "."), 
    width = 8000, height = 6000, res = 300)
VlnPlot(data, assay = "SCT", group.by = "variant", 
        features = sc.markers$gene, ncol = 8, pt.size = 0) + 
  theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
  xlab("") + 
  ylab("")
dev.off()


Idents(data) <- "Visium_clusters"
#' FeaturePlot of reference markers of merged data
png(file = paste("featureplot_ref_umap.png", sep = "."), 
width = 9000, height = 5000, res = 300)
FeaturePlot(data, 
            features = sc.markers$gene, 
            reduction = "umap_harmony_SCT", 
            slot = "data", 
            ncol= 8, 
            order = T, 
            label = T, 
            label.size = 2, 
            label.color = "black", 
            repel = T) & 
  theme(axis.text.x=element_text(angle = -45, hjust = 0)) & 
  xlab("Harmony:UMAP_1") & 
  ylab("Harmony:UMAP_2") &
  labs(color = "Data")
dev.off()

png(file = paste("featureplot_ref_uwot.png", sep = "."), 
    width = 9000, height = 5000, res = 300)
FeaturePlot(data, 
            features = sc.markers$gene, 
            reduction = "uwot_harmony_SCT", 
            slot = "data", 
            ncol= 8, 
            order = T, 
            label = T, 
            label.size = 2, 
            label.color = "black", 
            repel = T) & 
  theme(axis.text.x=element_text(angle = -45, hjust = 0)) &
  xlab("Harmony:UWOT_1") & 
  ylab("Harmony:UWOT_2") &
  labs(color = "Data")
dev.off()

#' plot individual sample reference gene expression
cscale <- c("darkblue", "blue", "lightblue", "white", 
                      "lightgray", "mistyrose", "red", "darkred", "black")
                    
#' SDC1, MS4A1, KLRC1, TPSB2 and SOX2 doesn't express in all individual sample
Idents(data) <- "sample"
#data <- ScaleData(data, assay = "SCT", features = ref.markers)
id <- levels(data)
for (i in seq_along(id)) {
  out.dir <- paste(deg.dir, file.name, cluster.ident, id[i], sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  
#' plot reference gene spatially
for (ftr in seq_along(sc.markers$gene)) {
  #' the data slot is log1p of re-corrected counts after PreSCTFindMarkers
   p1 <- FeatureOverlay(data,
                        features = sc.markers$gene[ftr], 
                        sampleids = i,
                        pt.size = 1,
                        pt.alpha = 0.5,
                        pt.border = F,
                        show.sb = F) +
     theme(plot.title = element_blank(), plot.subtitle = element_blank()) +
     labs(fill = "Data")
  
  #' the scale.data is the residual from the SCT models
  p2 <- ST.FeaturePlot(data,
                       features = sc.markers$gene[ftr], 
                       indices = i,
                       pt.size = 1,
                       pt.border = F,
                       slot = "scale",
                       center.zero = T,
                       cols = cscale,
                       show.sb = F) +
    theme(plot.title = element_blank(), plot.subtitle = element_blank()) +
    labs(fill = "ScaleData")
    
 
  p3 <- VlnPlot(data, assay = "SCT", idents = id[i], 
                features = sc.markers$gene[ftr]) + 
    theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
    xlab("") + ylab("") + NoLegend()
  
  cells <- WhichCells(data, idents = id[i])
  p4 <- FeaturePlot(data, 
                    features = sc.markers$gene[ftr], 
                    reduction = "uwot_harmony_SCT",
                    cells = cells, 
                    slot = "data", 
                    order = T, 
                    label = T, 
                    label.size = 2, 
                    label.color = "black", 
                    repel = T
                    ) + 
    theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
    xlab("Harmony:UWOT_1") + 
    ylab("Harmony:UWOT_2") +
    labs(color = "Data")

  png(file = paste0(paste("deg", id[i], ftr, sc.markers$gene[ftr], sep = "_"), 
                    ".png"), width = 4500, height = 4000, res = 500)
  print(    
    plot_grid(p3, p2, p4, p1, ncol = 2, labels = c("A", "C", "B", "D"))
  )
  dev.off() 
}
  png(file = paste0(paste("dotplot", id[i], "ref_by_cluster", sep = "_"), 
                    ".png"), width = 3500, height = 1200, res = 300)
  print(
    DotPlot(data, assay = "SCT", features = sc.markers$gene, idents = id[i], 
            group.by = "Visium_clusters", cluster.idents = F) + 
      theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
      xlab("") + ylab("")
  )
  dev.off()
}
