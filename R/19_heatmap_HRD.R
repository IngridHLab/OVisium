#' Utilize Scillus a Seurat wrapper for enhanced processing and visualization of
#' scRNA-seq data https://scillus.netlify.app/
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/DoMultiBarHeatmap.R", sep = "/"))


file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"
data <- readRDS(paste(rds.dir, paste0(file.name, "_deg.rds"), sep = "/"))

#' get cell cycle annotation
setwd(deg.dir)
exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes)

#' select variable features
setwd(deg.dir)
HRD_gene_signature <- readr::read_csv("HRD_gene_signature.csv", col_names = FALSE)
HRD_gene_signature<-t(HRD_gene_signature) %>% base::as.character() %>% na.omit() 


#' Original Seurat Heatmap on all clusters
clust.col <- 
  c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#F0027F","#BF5B17",
             "#666666","#1B9E77","#D95F02","#386CB0","black")
             
names(clust.col) <-c("0_FTE","1_FTE", "3_Mix", "2_Stroma",
                     "5_Stroma", "6_Stroma", "7_Stroma",
                     "8_Stroma", "9_Stroma", "10_Stroma",
                     "11_Stroma")
anno.col = list(
  Cluster = clust.col,
  Tissue = c(Fimbrial = "#FC2947", Otherside = "#FE6244", Proximal = "#FFDEB9"),
  Phase = c("S" = "#071952", "G1" = "#0B666A", "G2M" = "#97FEED"),
  Variant = c(c.3048_3052dup = "#900C3F", c.68_69del = "#FF4949", c.1082_1092del = "#FF8D29"),
  Age = c("35+" = "#EBE76C", "40+" = "#F0B86E", "45+" = "#ED7B7B", "50+" = "#836096"))

col_fun = circlize::colorRamp2(c(-2.5, 0, 2.5), c("#FF00FF", "black", "#FFFF00"))

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)

cols.use = list(
  Visium_clusters = clust.col,
  origin = c(Fimbrial = "#FC2947", Otherside = "#FE6244", Proximal = "#FFDEB9"),
  Phase = c("S" = "#071952", "G1" = "#0B666A", "G2M" = "#97FEED"),
  variant = c(c.3048_3052dup = "#900C3F", c.68_69del = "#FF4949", c.1082_1092del = "#FF8D29"),
  age = c("35+" = "#EBE76C", "40+" = "#F0B86E", "45+" = "#ED7B7B", "50+" = "#836096"))

#' ComplexHeatmap
#' only 144 HRD genes have SCT residuals
#' for individual patients
for (p in c(1,3,6,10)) {
  Idents(data) <- "patient"
  pat.data <- SubsetSTData(data, idents = p)
  patient <- paste("Case", p, sep = "_")

pat.data <- GetResidual(pat.data, features = HRD_gene_signature, assay = "SCT")
base::setdiff(HRD_gene_signature, row.names(pat.data@assays[["SCT"]]@scale.data))
pat.data.m <- as.matrix(GetAssayData(pat.data, slot = "scale.data", assay = "SCT"))
pat.data.m.HRD <- subset(pat.data.m, rownames(pat.data.m) %in% HRD_gene_signature)
HRD_gene_signature_sub <- row.names(pat.data.m.HRD)
quantile(pat.data.m.HRD, c(0.1, 0.95))

ha = HeatmapAnnotation(Cluster = pat.data@meta.data$Visium_clusters, 
                       Tissue = pat.data@meta.data$origin,
                       Phase = pat.data@meta.data$Phase,
                       Variant = pat.data@meta.data$variant,
                       Age = pat.data@meta.data$age, 
                       col = anno.col)
set.seed(1220)
png(file = paste(patient, "cHeatmap", length(HRD_gene_signature_sub), 
                 "HRD_signature.png", sep = "_"), 
                 width = 5000 , height = 5000, res = 400) 
p <-Heatmap(pat.data.m.HRD, 
        name = "Expression", 
        column_km = 3,
        column_km_repeats = 1000,
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 20),
        col = col_fun,
        cluster_rows = TRUE,
        row_km = 3,
        row_km_repeats = 1000,
        row_names_side = "left",
        show_row_dend = TRUE,
        row_names_gp = gpar(fontsize = 6),
        top_annotation = ha,
        use_raster = TRUE,
        raster_quality = 4)
p <- draw(p, row_title = "HRD signatures", column_title = patient,
          row_dend_side = "right")
dev.off()

#' extract heatmap order to plot heatmap without clustering.
hm_row_ord <-row_order(p) %>% list_c()
HRD_gene_signature_sub_reord <- HRD_gene_signature_sub[c(hm_row_ord)]
p <- DoMultiBarHeatmap(pat.data,
                       features = HRD_gene_signature_sub_reord, 
                       group.by = "origin",
                       cols.use = cols.use,
                       additional.group.by = c("Visium_clusters"),
                       additional.group.sort.by = c("Visium_clusters"),
                       disp.max = 3) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 6)) + 
  guides(color="none") +
  scale_fill_gradientn(colors = rev(mapal))

png(file = paste(patient, "mbHeatmap", length(HRD_gene_signature_sub), 
                 "HRD_signature.png", sep = "_"), 
    width = 5000 , height = 5000, res = 400) 
print( 
  p )
dev.off()
}



#' subset those patients who have pair samples
Idents(data) <- "patient"
pat.data <- SubsetSTData(data, idents = c(1,3,6,10))

#' 86 genes are not in the scale.data and 39 genes are not in the counts data 
pat.data <- GetResidual(pat.data, features = HRD_gene_signature, assay = "SCT")
base::setdiff(HRD_gene_signature, row.names(pat.data@assays[["SCT"]]@scale.data))
pat.data.m <- as.matrix(GetAssayData(pat.data, slot = "scale.data", assay = "SCT"))
pat.data.m.HRD <- subset(pat.data.m, rownames(pat.data.m) %in% HRD_gene_signature)
HRD_gene_signature_sub <- row.names(pat.data.m.HRD)
quantile(pat.data.m.HRD, c(0.1, 0.95))
  
  ha = HeatmapAnnotation(Cluster = pat.data@meta.data$Visium_clusters, 
                         Tissue = pat.data@meta.data$origin,
                         Phase = pat.data@meta.data$Phase,
                         Variant = pat.data@meta.data$variant,
                         Age = pat.data@meta.data$age, 
                         col = anno.col)
  set.seed(1220)
  png(file = paste("All_Cases_cHeatmap", length(HRD_gene_signature_sub), 
                   "HRD_signature.png", sep = "_"), 
      width = 5000 , height = 5000, res = 400) 
  p <-Heatmap(pat.data.m.HRD, 
              name = "Expression", 
              column_km = 3,
              column_km_repeats = 1000,
              cluster_columns = TRUE,
              show_column_dend = FALSE,
              show_column_names = FALSE,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 20),
              col = col_fun,
              cluster_rows = TRUE,
              row_km = 3,
              row_km_repeats = 1000,
              row_names_side = "left",
              show_row_dend = TRUE,
              row_names_gp = gpar(fontsize = 6),
              top_annotation = ha,
              use_raster = TRUE,
              raster_quality = 4)
  p <- draw(p, row_title = "HRD signatures", column_title = "Combined 4 Cases",
            row_dend_side = "right")
  dev.off()
  
  #' Extract row order and plot the most de genes for multibarheatmap
  hm_row_ord <- row_order(p)[c("2","3")] %>% list_c()
  hm_row_ord <- row_order(p) %>% list_c()
  HRD_gene_signature_sub_reord <- HRD_gene_signature_sub[c(hm_row_ord)]
  P <- DoMultiBarHeatmap(pat.data,
                         features = HRD_gene_signature_sub_reord, 
                         group.by = "Tissue_origin",
                         cols.use = cols.use,
                         additional.group.by = c("Visium_clusters"),
                         additional.group.sort.by = c("Visium_clusters"),
                         disp.max = 3) +
    theme(text = element_text(size = 10)) +
    theme(axis.text = element_text(size = 10)) + 
    guides(color="none") +
    scale_fill_gradientn(colors = rev(mapal))
  
  png(file = paste("All_Cases_mbHeatmap", length(HRD_gene_signature_sub_reord), 
                   "HRD_signature.png", sep = "_"), 
      width = 3000 , height = 3000, res = 300) 
  print( 
    P )
  dev.off()
  
  #' vlnplot of all features
  png(file = paste("All_Cases_vlnplot", length(HRD_gene_signature_sub_reord),
                   "HRD_signature_by_cluster.png", sep = "_"), 
      width = 8000, height = 8000, res = 100)
  VlnPlot(pat.data, assay = "SCT", 
          group.by = "Tissue_origin",
          features = HRD_gene_signature_sub_reord,
          ncol = 12,
          pt.size = 0) & 
    xlab("") & 
    ylab("")
  dev.off()
  
  #' plot only those are different between fimbrial and proximal
  dfeatures <- c("CPE", "DPYSL3", "LAMB2","CRYAB","FBLN1","STAT2", "FOXO3", "CTSC")
  png(file = paste("All_Cases_vlnplot", length(dfeatures),
                   "HRD_signature_by_cluster.png", sep = "_"), 
      width = 800, height = 1200, res = 100)
  VlnPlot(pat.data, assay = "SCT", 
          group.by = "Visium_clusters",
          split.by = "Tissue_origin",
          features = dfeatures,
          ncol = 2,
          pt.size = 0,
          log = T) +
    theme(legend.position = "bottom") &
    xlab("") &
    ylab("")
  dev.off()