library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))


file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"

#' Normalized by median UMI counts or seuqncing depth
data.deg <- readRDS(paste(rds.dir, paste0(file.name, "_deg.rds"), sep = "/"))
#' Batch corrected by HarmonyMatrix
load("Variable_features_filt_SCT_log2counts+1_harmony.RData")

#' Choose the data file 
data <- dat.deg  #' data <- data.sub.filt
source(paste(home, "OVisium/manuscript/gitHub/ComplexHeatmap_settings.R", sep = "/"))

#' copy csv files to one place
#' top 20 gene list
setwd(paste(deg.dir, file.name, cluster.ident, "Variable_features_filt", "MAST", "FindMarkers", sep = "/"))
combined_filtered_markers_files <- list.files(path = "./Combined_filtered_markers_v1/csv", pattern = "*top20.csv")
markers_list <- list()
for (f in 1:length(combined_filtered_markers_files)) {
  n <- gsub(".top20.csv", "", combined_filtered_markers_files[[f]])
  n <- gsub("Cluster.", "", n)
  combined_filtered_markers <- 
    readr::read_csv(paste("./Combined_filtered_markers_v1/csv", combined_filtered_markers_files[[f]], sep = "/"), col_names = T)
  markers_list[[n]] <- combined_filtered_markers
}

#' plot complex heatmap
for (i in 1:length(markers_list)) {
  cluster_markers <- markers_list[[i]]
  cluster_markers$id1 <- factor(cluster_markers$id1, 
                                    levels = c("0_FTE", "1_FTE", "3_Mix", 
                                               "2_Stroma", "5_Stroma", 
                                               "6_Stroma", "7_Stroma",
                                               "8_Stroma", "9_Stroma", 
                                               "10_Stroma", "11_Stroma"))
  cluster_markers$id2 <- factor(cluster_markers$id2, 
                                levels = c("0_FTE", "1_FTE", "3_Mix", 
                                           "2_Stroma", "5_Stroma", 
                                           "6_Stroma", "7_Stroma",
                                           "8_Stroma", "9_Stroma", 
                                           "10_Stroma", "11_Stroma"))
  
#' bulk analysis clusters by tissue
  all.pat.data.bulk <- 
    AverageExpression(data, features = cluster_markers$gene,
                      group.by = c("Visium_clusters","Tissue_origin"),
                      assays = "SCT", return.seurat = T)
  all.pat.data.bulk.m <- as.matrix(GetAssayData(all.pat.data.bulk, 
                                                slot = "scale.data", 
                                                assay = "SCT"))
  
  all.pat.ha <- 
    HeatmapAnnotation(Cluster = factor(str_sub(colnames(all.pat.data.bulk.m), end = -10), 
                                       levels = c("0_FTE", "1_FTE", "3_Mix", "2_Stroma",
                                                  "5_Stroma", "6_Stroma", "7_Stroma",
                                                  "8_Stroma", "9_Stroma", "10_Stroma",
                                                  "11_Stroma")),
                      Tissue = gsub(".*_", "", colnames(all.pat.data.bulk.m)),
                      col = anno.col)
  
  set.seed(1220)
  png(file = paste("./Combined_filtered_markers_v1/csv/AverageExpression", names(markers_list[i]), 
                   "top20_clusterbytissue.png", sep = "_"), 
      width = 6000 , height = 10000, res = 400) 
  p <-Heatmap(all.pat.data.bulk.m, 
              name = "Expression", 
              cluster_columns = T,
              show_column_dend = T,
              show_column_names = T,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 12),
              column_title = "All Cases",
              col = puryel, # or col = rev(redblu) OR puryel
              cluster_rows = T,
              cluster_row_slices = F,
              row_names_side = "left",
              show_row_dend = T,
              row_names_gp = gpar(fontsize = 8),
              row_title_gp = gpar(fontsize = 12),
              row_title_rot = 0,
              row_title_side = "right", 
              row_split = cluster_markers$id1, #id2 for individual
              top_annotation = all.pat.ha,
              use_raster = TRUE,
              raster_quality = 4)
  p <- draw(p)
  dev.off()
  
  #' bulk analysis patients by tissue
  all.pat.data.bulk.2 <- 
    AverageExpression(data, 
                      features = cluster_markers$gene,
                      group.by = c("patient","origin"),
                      assays = "SCT", return.seurat = T)
  all.pat.data.bulk.2.m <- as.matrix(GetAssayData(all.pat.data.bulk.2, 
                                                  slot = "scale.data", 
                                                  assay = "SCT"))
  all.pat.ha.2 <- 
    HeatmapAnnotation(Patient = gsub("_.*", "", colnames(all.pat.data.bulk.2.m)),
                      Tissue = gsub(".*_", "", colnames(all.pat.data.bulk.2.m)),
                      col = anno.col)
  
  set.seed(1220)
  png(file = paste("./Combined_filtered_markers_v1/csv/AverageExpression", names(markers_list[i]), 
                   "top20_patientbytissue.png", sep = "_"), 
      width = 6000 , height = 10000, res = 400) 
  p <-Heatmap(all.pat.data.bulk.2.m, 
              name = "Expression", 
              cluster_columns = T,
              show_column_dend = T,
              show_column_names = T,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 12),
              column_title = "All Cases",
              col = puryel,
              cluster_rows = T,
              cluster_row_slices = F,
              row_names_side = "left",
              show_row_dend = T,
              row_names_gp = gpar(fontsize = 8),
              row_title_gp = gpar(fontsize = 12),
              row_title_rot = 0,
              row_title_side = "right", 
              row_split = cluster_markers$id1,
              top_annotation = all.pat.ha.2,
              use_raster = TRUE,
              raster_quality = 4)
  p <- draw(p)
  dev.off()
  
  #' bulk analysis individual samples by slides
  all.pat.data.bulk.3 <- 
    AverageExpression(data, 
                      features = cluster_markers$gene,
                      group.by = c("library", "slide"),
                      assays = "SCT", return.seurat = T)
  all.pat.data.bulk.3.m <- as.matrix(GetAssayData(all.pat.data.bulk.3, 
                                                  slot = "scale.data", 
                                                  assay = "SCT"))
  all.pat.ha.3 <- 
    HeatmapAnnotation(Sample = factor(str_sub(colnames(all.pat.data.bulk.3.m), end = -12), 
                                      levels = c("PO2", "PO2_2re", "PO6", "PO7_A1", "PO7_2", "PO9_B1", 
                                                 "PO13", "PO13a", "PO20", "PO20a", "PO20_2", "PO28", "PO35", 
                                                 "PO37", "PO40", "PO40_2", "PO41", "PO45")),
                      Batch = gsub(".*_", "", colnames(all.pat.data.bulk.3.m)),  
                      col = anno.col)
  
  set.seed(1220)
  png(file = paste("./Combined_filtered_markers_v1/csv/AverageExpression", names(markers_list[i]), 
                   "top20_samplebyslide.png", sep = "_"), 
      width = 6000 , height = 10000, res = 400) 
  p <-Heatmap(all.pat.data.bulk.3.m, 
              name = "Expression", 
              cluster_columns = T,
              show_column_dend = T,
              show_column_names = T,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 12),
              column_title = "All Cases",
              col = puryel,
              cluster_rows = T,
              cluster_row_slices = F,
              row_names_side = "left",
              show_row_dend = T,
              row_names_gp = gpar(fontsize = 8),
              row_title_gp = gpar(fontsize = 12),
              row_title_rot = 0,
              row_title_side = "right", 
              row_split = cluster_markers$id1,
              top_annotation = all.pat.ha.3,
              use_raster = TRUE,
              raster_quality = 4)
  p <- draw(p)
  dev.off()
  
  #' bulk analysis individual cluster by run
  all.pat.data.bulk.4 <- 
    AverageExpression(data, 
                      features = cluster_markers$gene,
                      group.by = c("library", "run"),
                      assays = "SCT", return.seurat = T)
  all.pat.data.bulk.4.m <- as.matrix(GetAssayData(all.pat.data.bulk.4, 
                                                  slot = "scale.data", 
                                                  assay = "SCT"))
  all.pat.ha.4 <- 
    HeatmapAnnotation(Sample = factor(str_sub(colnames(all.pat.data.bulk.4.m), end = -14), 
                                      levels = c("PO2", "PO2_2re", "PO6", "PO7_A1", "PO7_2", "PO9_B1", 
                                                 "PO13", "PO13a", "PO20", "PO20a", "PO20_2", "PO28", "PO35", 
                                                 "PO37", "PO40", "PO40_2", "PO41", "PO45")),
                      SeqRun = factor(str_sub(colnames(all.pat.data.bulk.4.m), start = -12), 
                                      levels = c("CTG_2021_031","CTG_2021_075", "CTG_2021_099")),
                      col = anno.col)
  
  set.seed(1220)
  png(file = paste("./Combined_filtered_markers_v1/csv/AverageExpression", names(markers_list[i]), 
                   "top20_samplebyseqrun.png", sep = "_"), 
      width = 6000 , height = 10000, res = 400) 
  p <-Heatmap(all.pat.data.bulk.4.m, 
              name = "Expression", 
              cluster_columns = T,
              show_column_dend = T,
              show_column_names = T,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 12),
              column_title = "All Cases",
              col = puryel,
              cluster_rows = T,
              cluster_row_slices = F,
              row_names_side = "left",
              show_row_dend = T,
              row_names_gp = gpar(fontsize = 8),
              row_title_gp = gpar(fontsize = 12),
              row_title_rot = 0,
              row_title_side = "right", 
              row_split = cluster_markers$id1,
              top_annotation = all.pat.ha.4,
              use_raster = TRUE,
              raster_quality = 4)
  p <- draw(p)
  dev.off()
  
  #' bulk analysis individual slides
  all.pat.data.bulk.5 <- 
    AverageExpression(data, 
                      features = cluster_markers$gene,
                      group.by = c("slide"),
                      assays = "SCT", return.seurat = T)
  all.pat.data.bulk.5.m <- as.matrix(GetAssayData(all.pat.data.bulk.5, 
                                                  slot = "scale.data", 
                                                  assay = "SCT"))
  all.pat.ha.5 <- 
    HeatmapAnnotation(Batch = colnames(all.pat.data.bulk.5.m),
                      col = anno.col)
  
  set.seed(1220)
  png(file = paste("./Combined_filtered_markers_v1/csv/AverageExpression", names(markers_list[i]), 
                   "top20_byslide.png", sep = "_"), 
      width = 6000 , height = 10000, res = 400) 
  p <-Heatmap(all.pat.data.bulk.5.m, 
              name = "Expression", 
              cluster_columns = T,
              show_column_dend = T,
              show_column_names = T,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 12),
              column_title = "All Cases",
              col = puryel,
              cluster_rows = T,
              cluster_row_slices = F,
              row_names_side = "left",
              show_row_dend = T,
              row_names_gp = gpar(fontsize = 8),
              row_title_gp = gpar(fontsize = 12),
              row_title_rot = 0,
              row_title_side = "right", 
              row_split = cluster_markers$id1,
              top_annotation = all.pat.ha.5,
              use_raster = TRUE,
              raster_quality = 4)
  p <- draw(p)
  dev.off()
  }
