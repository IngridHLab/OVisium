#' https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
#' https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html

clust.col <- 
  c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#F0027F","#BF5B17",
             "#666666","#1B9E77","#D95F02","#386CB0","black")
             
names(clust.col) <-c("0_FTE","1_FTE", "3_Mix", "2_Stroma",
                     "5_Stroma", "6_Stroma", "7_Stroma",
                     "8_Stroma", "9_Stroma", "10_Stroma",
                     "11_Stroma")
patient.col <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(12) 
names(patient.col) <- factor(1:12)
sample.col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(18) 
names(sample.col) <- c("PO2", "PO2_2re", "PO6", "PO7_A1", "PO7_2", "PO9_B1", 
                       "PO13", "PO13a", "PO20", "PO20a", "PO20_2", "PO28", "PO35", 
                       "PO37", "PO40", "PO40_2", "PO41", "PO45")
batch.col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(5) 
names(batch.col) <- c("V10S21-048" ,"V10S21-050","V10T06-031","V10T06-110","V10U22-108") 

anno.col = list(
  Cluster = clust.col,
  Sample = sample.col,
  Patient = patient.col,
  Tissue = c(Fimbrial = "#FC2947", Otherside = "#FE6244", Proximal = "#FFDEB9"),
  Age = c("35+" = "#EBE76C", "40+" = "#F0B86E", "45+" = "#ED7B7B", "50+" = "#836096"),
  Batch = batch.col
)

puryel = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00"))

#' Seurat assay data and annotation for individual spots
cellInfo <- data.frame(Visium_clusters=Idents(data))
cellInfo$Visium_clusters <- factor(cellInfo$Visium_clusters, 
                                   levels = c("0_FTE", "1_FTE", "3_Mix", 
                                              "2_Stroma", "5_Stroma", 
                                              "6_Stroma", "7_Stroma",
                                              "8_Stroma", "9_Stroma", 
                                              "10_Stroma", "11_Stroma"))
cellInfo = cellInfo %>% dplyr::arrange(Visium_clusters)
ha = HeatmapAnnotation(Cluster = data@meta.data$Visium_clusters, 
                       Sample = data@meta.data$library,
                       Patient = data@meta.data$patient,
                       Tissue = data@meta.data$origin,
                       Age = data@meta.data$age,
                       Batch = data@meta.data$slide,
                       col = anno.col,
                       annotation_name_gp= gpar(fontsize = 18),
                       annotation_legend_param = list(
                         Cluster = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14)),
                         Sample = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14)),
                         Patient = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14)),
                         Tissue = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14)),
                         Age = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14)),
                         Batch = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14))
                       )
)