#' Utilize Scillus a Seurat wrapper for enhanced processing and visualization of
#' scRNA-seq data https://scillus.netlify.app/
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/17_geneAnnotations.R", sep = "/"))


file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"
data <- readRDS(paste(rds.dir, paste0(file.name, "_deg.rds"), sep = "/"))

#' Prepare differential gene expression data 
#' only plot positive and top 100 genes for better heatmap
deg.list <- 
  c("Variable_features","Variable_a500UMI_features",
    "Scaled_features","Scaled_a500UMI_features","All_features")
test <- c("wilcox", "MAST")

for (i in seq_along(deg.list)) {
  for (t in seq_along(test)) {
out.dir<-paste(deg.dir, file.name, cluster.ident, deg.list[i], test[t], 
               sep = "/")
setwd(out.dir)
all.markers.filt <- readr::read_csv("Clusters.markers.filt.csv")
all.markers.filt$cluster<-
  factor(all.markers.filt$cluster, 
         levels = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                    "5_Stroma", "6_Stroma", "7_Stroma",
                    "8_Stroma", "9_Stroma", "10_Stroma",
                    "11_Stroma"))

pos.markers.filt <- all.markers.filt %>% 
  dplyr::filter(avg_log2FC > 1 )

neg.markers.filt <- all.markers.filt %>% 
  dplyr::filter(avg_log2FC <= 1 )

top <- pos.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 100)

write.table(top, 
            file = "heatmap_pos_avgFC_top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)

#' Heatmap on all clusters
clust.col <- 
  c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#F0027F","#BF5B17",
             "#666666","#1B9E77","#D95F02","#386CB0","black")
png("heatmap_pos_top100.png", width = 8000 , height = 9000, res = 400) 
print(
  p<-DoHeatmap(data, assay = "SCT", features = top$gene, 
            group.by = "Visium_clusters", 
            group.colors = clust.col, 
            size = 8, 
            angle = 70) +
    theme(text = element_text(size = 18)) +
    theme(axis.text = element_text(size = 4)) + 
    guides(color="none") +
    labs(fill = "SCT:Residual")
)
dev.off()

#' Heatmap on the same gene set but split fimbrial from proximal
Idents(data) <- "Tissue_origin"
cells.1 <- WhichCells(data, idents = "Fimbrial")
cells.2 <- WhichCells(data,  idents = "Proximal")

png("heatmap_pos_top100_fimbrial.png", width = 8000, height = 9000, res = 400) 
print(
  DoHeatmap(data, assay = "SCT", features=top$gene, 
            cells = cells.1, 
            group.by = "Visium_clusters", 
            group.colors = clust.col, 
            size = 8, 
            angle = 70) +
    theme(text = element_text(size = 18)) +
    theme(axis.text = element_text(size = 4)) + 
    guides(color="none") +
    labs(fill = "SCT:Residual")
)
dev.off()

png("heatmap_pos_top100_proximal.png", width = 8000, height = 9000, res = 400) 
print(
  DoHeatmap(data, assay = "SCT", features=top$gene, 
            cells = cells.2,
            group.by = "Visium_clusters", 
            group.colors = clust.col,
            size = 9, 
            angle = 70) +
    theme(text = element_text(size = 18)) +
    theme(axis.text = element_text(size = 4)) + 
    guides(color="none") +
    labs(fill="SCT:Residual")
)
dev.off()

#' The following features were omitted as they were not found in the scale.data 
#' slot for the SCT assay: FCN1, ADAMTS4, ECEL1
  }
}
