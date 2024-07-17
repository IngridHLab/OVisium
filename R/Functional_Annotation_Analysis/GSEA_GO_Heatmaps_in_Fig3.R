# Generate Figure 3 in the manuscript

library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Get_Gene_Annotation.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/DoMultiBarHeatmap.R", sep = "/"))

file.name <- "OVisium_SCT_merged" 
cluster.ident <- "harmony_SCT_res_0.6"
setwd(rds.dir)
load("Variable_features_filt_SCT_log2counts+1_harmony.RData")

#' prepare the data slot for heatmap
#' for scale.data analysis, no filtering will be applied on the log2FC
data.sub.filt@assays[["SCT"]]@scale.data <- as.matrix(data.sub.filt@assays[["SCT"]]@data)
data.sub.filt <- ScaleData(data.sub.filt)

feature.list <- list(vfeatures.filt) # vfeatures 
names(feature.list) <- c("Variable_features_filt") # "Variable_features"
test <- c("MAST") # "wilcox", "MAST"

data <- data.sub.filt

#' rename slide to batch
Idents(data) <- "slide"
data <- RenameIdents(data, "V10S21-048" = "1")
data <- RenameIdents(data, "V10S21-050" = "2")
data <- RenameIdents(data, "V10T06-031" = "3")
data <- RenameIdents(data, "V10T06-110" = "4")
data <- RenameIdents(data, "V10U22-108" = "5")
Idents(data) <- factor(Idents(data),
                       levels = c(1:5))
data[["batch"]] <- Idents(data)

#' reorder patient
Idents(data) <- "patient"
Idents(data) <- factor(Idents(data),
                       levels = c(1:12))
data[["patient"]] <- Idents(data)

#' heatmap settings
clust.col <- 
  c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#F0027F","#BF5B17",
             "#666666","#1B9E77","#D95F02","#386CB0","black")
             
names(clust.col) <-c("0_FTE","1_FTE", "3_Mix", "2_Stroma",
                     "5_Stroma", "6_Stroma", "7_Stroma",
                     "8_Stroma", "9_Stroma", "10_Stroma",
                     "11_Stroma")
patient.col <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(12) 
names(patient.col) <-  c(1:12)
batch.col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(5) 
names(batch.col) <- c(1:5)

anno.col = list(
  Cluster = clust.col,
  Patient = patient.col,
  Batch = batch.col
)

puryel = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00"))

#' Seurat assay data and annotation for individual spots
Idents(data) <- "Visium_clusters"
cellInfo <- data.frame(Visium_clusters=Idents(data))
cellInfo$Visium_clusters <- factor(cellInfo$Visium_clusters, 
                                   levels = c("0_FTE", "1_FTE", "3_Mix", 
                                              "2_Stroma", "5_Stroma", 
                                              "6_Stroma", "7_Stroma",
                                              "8_Stroma", "9_Stroma", 
                                              "10_Stroma", "11_Stroma"))
cellInfo = cellInfo %>% dplyr::arrange(Visium_clusters)
ha = HeatmapAnnotation(Cluster = data@meta.data$Visium_clusters, 
                       Patient = data@meta.data$patient,
                       Batch = data@meta.data$batch,
                       col = anno.col,
                       simple_anno_size = unit(1, "cm"),
                       annotation_name_gp= gpar(fontsize = 18),
                       annotation_legend_param = list(
                         Cluster = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14)),
                        
                         Patient = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14)),
                         
                         Batch = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14))
                       )
)

combined.markers.filt.merged.top20 <- 
  readr::read_csv("~/OVisium/DE_functional_analysis/OVisium_SCT_merged/harmony_SCT_res_0.6/Variable_features_filt/MAST/FindMarkers/Combined_filtered_markers_v1/2024-04-20_harmony_sample_log2SCTcounts1/csv/Cluster.all.once.top20.csv")

combined.markers.filt.merged.top20$id1 <- factor(combined.markers.filt.merged.top20$id1, 
                               levels = c("0_FTE", "1_FTE", "3_Mix", 
                                          "2_Stroma", "5_Stroma", 
                                          "6_Stroma", "7_Stroma",
                                          "8_Stroma", "9_Stroma", 
                                          "10_Stroma", "11_Stroma"))
combined.markers.filt.merged.top20$id2 <- factor(combined.markers.filt.merged.top20$id2, 
                               levels = c("0_FTE", "1_FTE", "3_Mix", 
                                          "2_Stroma", "5_Stroma", 
                                          "6_Stroma", "7_Stroma",
                                          "8_Stroma", "9_Stroma", 
                                          "10_Stroma", "11_Stroma"))

#' plot complex heatmap for panel A
mat.merged <- data[["SCT"]]@data[combined.markers.filt.merged.top20$gene, ] %>% as.matrix()
## scale the rows
mat.merged <- t(scale(t(mat.merged)))

set.seed(1220)
png(file = paste0("Figure_3A_heatmap_20240528.png"), width = 8000 , height = 8000, res = 500)
p <-Heatmap(mat.merged,
            name = "Log2 Expression",
            cluster_columns = F,
            column_order = row.names(cellInfo),
            show_column_dend = F,
            show_column_names = F,
            column_title_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 20),
            column_title = paste("FindMarkers All Clusters Top 20"),
            col = puryel,
            cluster_rows = F,
            row_names_side = "left",
            show_row_dend = F,
            row_names_gp = gpar(fontsize = 7),
            row_title_gp = gpar(fontsize = 18),
            row_title_rot = 0,
            row_title_side = "right",
            row_split = combined.markers.filt.merged.top20$id1,
            top_annotation = ha,
            heatmap_legend_param = list(
              legend_direction = "horizontal", 
              legend_width = unit(10, "cm"),
              lebal_gp = gpar(fontsize = 30)),
            use_raster = TRUE,
            raster_resize_mat = TRUE,
            raster_quality = 4)
  p <- draw(p, heatmap_legend_side="bottom")
  dev.off()



#' Functional annotation heatmap for Panel B & C
#' 
my_merge <- function(df1, df2){                                
  full_join(df1, df2, by = c("ID"))
}

GSEA_list <- readRDS("~/OVisium/DE_functional_analysis/OVisium_SCT_merged/harmony_SCT_res_0.6/Variable_features_filt/MAST/FindMarkers/Combined_filtered_markers_v1/2024-04-20_harmony_sample_log2SCTcounts1/GSEA/c5.go.v2023.2/rds/c5.go.v2023.2_Cluster.all.list.rds_result.rds")

result.up.list <- list()
result.down.list <- list()

for (i in 1:length(GSEA_list)) {
  cluster.name <- names(GSEA_list[i])
  result <- GSEA_list[[i]]@result[c("ID", "NES")]
  result$Cluster <- cluster.name
  result.up <- result %>% dplyr::filter(NES>0) 
  result.down <- result %>% dplyr::filter(NES<0) #%>% slice_max(order_by = abs(NES), n=5)
  colnames(result.up)[2] <- cluster.name
  rownames(result.up) <- NULL
  colnames(result.down)[2] <- cluster.name
  rownames(result.down) <- NULL
  result.up.list[[cluster.name]] <- result.up
  result.down.list[[cluster.name]] <- result.down
}

all.up <- Reduce(my_merge, result.up.list)
all.up[is.na(all.up)] <- 0
write.table(all.up, "Figure_3B.csv", quote = F, row.names = F, col.names = T, sep = ",")
#' manually combine Cluster columns
all.up <- readr::read_csv("Figure_3B_all_annotations.csv")
all.up$Cluster <- factor(all.up$Cluster, levels = c("0_FTE",  
                                                    "2_Stroma", "5_Stroma", 
                                                    "6_Stroma", "7_Stroma",
                                                    "8_Stroma", "9_Stroma", 
                                                    "10_Stroma", "11_Stroma"))
all.up.top5 <- all.up %>% dplyr::group_by(Cluster) %>%
  dplyr::slice_head(n = 5)
all.up.top5 <- column_to_rownames(all.up.top5, var = "ID")
all.up.top5[,1] <- NULL


all.down <- Reduce(my_merge, result.down.list)
all.down[is.na(all.down)] <- 0
write.table(all.down, "Figure_3C.csv", quote = F, row.names = F, col.names = T, sep = ",")
all.down <- readr::read_csv("Figure_3C_all_annotations.csv")
all.down$Cluster <- factor(all.down$Cluster, levels = c("0_FTE", "1_FTE", "2_Stroma", "5_Stroma", "7_Stroma",
                                                    "8_Stroma", "10_Stroma", "11_Stroma"))
all.down.top5 <- all.down %>% dplyr::group_by(Cluster) %>%
  dplyr::slice_head(n = 5)
all.down.top5 <- column_to_rownames(all.down.top5, var = "ID")
all.down.top5[,1] <- NULL


gsea.anno.col = list(
  Cluster = clust.col
)
gsea.puryel = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00"))

Cluster = factor(c("0_FTE","1_FTE", "3_Mix", "2_Stroma",
                   "5_Stroma", "6_Stroma", "7_Stroma",
                   "8_Stroma", "9_Stroma", "10_Stroma",
                   "11_Stroma"), levels = c("0_FTE", "1_FTE", "3_Mix", 
                                              "2_Stroma", "5_Stroma", 
                                              "6_Stroma", "7_Stroma",
                                              "8_Stroma", "9_Stroma", 
                                              "10_Stroma", "11_Stroma"))

gsea.ha <- HeatmapAnnotation(Cluster = Cluster, 
                             col = gsea.anno.col,
                             simple_anno_size = unit(0.8, "cm"),
                             annotation_name_gp= gpar(fontsize = 18),
                             annotation_legend_param = list(
                               Cluster = list(
                                 title_gp = gpar(fontsize = 14, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 14))
                               )
                             )

set.seed(1220)
png(file = paste0("Figure_3B_heatmap_20240529.png"), width = 5400 , height = 5100, res = 460)
p <-Heatmap(all.up.top5,
            name = "NES",
            cluster_columns = F,
            show_column_dend = F,
            show_column_names = F,
            column_title_gp = gpar(fontsize = 15),
            column_names_gp = gpar(fontsize = 20),
            column_title = paste("Top 5 Positive Enriched Annotations"),
            col = gsea.puryel,
            cluster_rows = F,
            row_names_side = "right",
            show_row_dend = F,
            row_names_gp = gpar(fontsize = 12),
            row_names_max_width = unit(24, "cm"),
            row_title_gp = gpar(fontsize = 18),
            row_title_rot = 0,
            row_title_side = "right",
            top_annotation = gsea.ha,
            heatmap_legend_param = list(
              legend_direction = "horizontal", 
              legend_width = unit(8, "cm")),
            use_raster = TRUE,
            raster_resize_mat = TRUE,
            raster_quality = 4)
p <- draw(p, heatmap_legend_side="bottom", annotation_legend_side="right")
dev.off()

set.seed(1220)
png(file = paste0("Figure_3C_heatmap_20240529.png"), width = 4750 , height = 3200, res = 415)
p <-Heatmap(all.down.top5,
            name = "NES",
            cluster_columns = F,
            show_column_dend = F,
            show_column_names = F,
            column_title_gp = gpar(fontsize = 15),
            column_names_gp = gpar(fontsize = 20),
            column_title = paste("Top 5 Negative Enriched Annotations"),
            col = gsea.puryel,
            cluster_rows = F,
            row_names_side = "right",
            show_row_dend = F,
            row_names_gp = gpar(fontsize = 12),
            row_names_max_width = unit(24, "cm"),
            row_title_gp = gpar(fontsize = 18),
            row_title_rot = 0,
            row_title_side = "right",
            top_annotation = gsea.ha,
            heatmap_legend_param = list(
              legend_direction = "horizontal", 
              legend_width = unit(8, "cm")),
            use_raster = TRUE,
            raster_resize_mat = TRUE,
            raster_quality = 4)
p <- draw(p, heatmap_legend_side="bottom", annotation_legend_side="right")
dev.off()
