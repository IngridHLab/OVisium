library(fs)
library(gplots)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

cell_col <-c("#A6CEE3","#1F78B4","#B8E986", "#7ED321", "#417505", "#FFFF99", "#E31A1C", 
             "#FB9A99", "#CAB2D6", "#6A3D9A", "#FDBF6F", "#B15928")

names(cell_col) <- c("ciliated epithelial cell",
                     "secretory cell",
                     "fibroblast",
                     "myofibroblast cell",
                     "smooth muscle cell", 
                     "pericyte",
                     "blood vessel endothelial cell",
                     "endothelial cell of lymphatic vessel",
                     "B cell",
                     "mature NK T cell",
                     "mast cell",
                     "macrophage")
# cell_col <-list(CellType = cell_col)
#' load and process scRNAseq reference data
#' prepare counts from scRNA data
#' prepare nUMI object
#' prepare celltype matrix
#' last create single cell reference object
#' Note: the integrated data for some reason doesn't work for all cell types
#' which might due to SCT counts but the merged data works
set.seed(1220)
setwd(rds.dir)
healthy.tubes = readRDS("healthyTubes_SCT_vIntegrated.rds")
Idents(healthy.tubes) <- "cell_type"
gene_names <- rownames(healthy.tubes)
Idents(healthy.tubes) <- factor(Idents(healthy.tubes),
                           levels = c("ciliated epithelial cell",
                                      "secretory cell",
                                      "fibroblast",
                                      "myofibroblast cell",
                                      "smooth muscle cell", 
                                      "pericyte",
                                      "blood vessel endothelial cell",
                                      "endothelial cell of lymphatic vessel",
                                      "B cell",
                                      "mature NK T cell",
                                      "mast cell",
                                      "macrophage"))

#' load out DEG data
setwd(deg.dir)
features.id <- readRDS("features_EntrezID.rds")
setwd("~/OVisium/R/Deconvolution/Correlation/All_DEGs")
deg.file <- list.files(pattern = "*all.list.rds")
name <- gsub(".list.*", "", deg.file)
markers.list <- readRDS(deg.file) 
markers<-bind_rows(markers.list) 

#' Convert gene name to entrezid
common <- inner_join(markers, features.id, by = c("gene" = "SYMBOL"), 
                     relationship = "many-to-many")

#' Check genes that don't have EntrezID
warning(paste(setdiff(markers$gene, common$gene), 
              "skipped; ", sep = " "))

#' order by foldchange and remove duplicate 
common_unique <-common[!duplicated(common$gene),]
common_unique_gene<- common_unique$ENSEMBL
common_unique_gene <- common_unique_gene[common_unique_gene %in% gene_names]

#' Perform pairwise correlations and hierarchical clustering for these gene sets
#' heatmap2 setting
random.matrix <- matrix(runif(500, min = -1, max = 1), nrow = 50)
quantile.range <- quantile(random.matrix, probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["35%"], quantile.range["75%"], 0.06)
color.palette <- colorRampPalette(c("#0571b0","#f7f7f7","#ca0020"))(length(palette.breaks)-1)
clustFunction<-function(x)
  hclust(as.dist(1-cor(t(as.matrix(x)), method = "pearson")), method = "average")
heatmapPearson<- function(correlations)
  heatmap.2(x = correlations,
            col = color.palette,
            breaks = palette.breaks,
            trace = "none", symm = T,
            cexRow = 1, cexCol = 1,
            hclustfun = clustFunction)

#' Compute the Pearson correlation coefficient on the logarithm of normalized 
#' expression of genes from “sc.features” in all cells from scRNA-seq subset
features <- common_unique_gene
sc_data <- ScaleData(healthy.tubes, features = features)
sc_data <- RunPCA(sc_data, features = features, npcs = 25)
elbow_plot<- ElbowPlot(sc_data, ndims = 50)
ggsave(filename = "visium_DEGs_all_elbow_plot.png", plot = elbow_plot, width = 3.75, height = 3.75, dpi = 300)
sc.features <- rownames(sc_data@reductions[["pca"]]@feature.loadings)
Idents(sc_data) <- "cell_type"
pca_plot <- DimPlot(sc_data, reduction = "pca", pt.size = 0.1, label = F, cols = cell_col)
pca_plot<- pca_plot & NoLegend()
ggsave(filename = "visium_DEGs_all_pca_plot.png", plot = pca_plot, width = 3.75, height = 3.75, dpi = 300)

#' calculate pearson correlation
mat<-as.data.frame(sc_data@assays[["integrated"]]@data[sc.features,]) %>% rownames_to_column("ENSEMBL")
mat <- inner_join(mat, features.id, relationship = "many-to-many")
row.names(mat) <- mat$SYMBOL
mat <- mat[-c(1,58815,58816)]

correlations_DEGs_log <- cor(method = "pearson",
                             log2(t(expm1(as.matrix(mat)))+1))

png(file = "visium_DEGs_all_corr_plot.png", width = 9000, height = 9000, res = 300)
heatmapPearson(correlations_DEGs_log)
dev.off()

#' subset DEGs with correlation >0.5
idx <- which(correlations_DEGs_log > 0.5 & lower.tri(correlations_DEGs_log, diag = F), arr.ind = TRUE)
high_corr <-cbind(rownames(idx), colnames(correlations_DEGs_log)[idx[, 2]]) %>% as.data.frame()
high_corr <-high_corr[!duplicated(high_corr$V1),]

#' plot those highly correlated genes on heatmap
Idents(healthy.tubes) <- "cell_type"
high_corr <- inner_join(high_corr, features.id, by = c("V1" = "SYMBOL"), 
                     relationship = "many-to-many")
saveRDS(high_corr, file = "High_corr_all_DEGs.rds")

#' perform PCA 
sc_data <- ScaleData(healthy.tubes, features = high_corr$ENSEMBL)
sc_data <- RunPCA(sc_data, features = high_corr$ENSEMBL, npcs = 25)
elbow_plot<- ElbowPlot(sc_data, ndims = 25)
ggsave(filename = "visium_DEGs_all_hi_corr_elbow_plot.png", plot = elbow_plot, width = 3.75, height = 3.75, dpi = 300)
pca_plot<- DimPlot(sc_data, reduction = "pca", pt.size = 0.1, label = F, cols = cell_col)
legend <- get_legend(pca_plot)
pca_plot<- pca_plot & NoLegend()
ggsave(filename = "visium_DEGs_all_hi_corr_pca_plot.png", plot = pca_plot, width = 3.75, height = 3.75, dpi = 300)
ggsave(filename = "cell_type_legend.png", plot = legend, width = 2, height = 3.75, dpi = 300)

data.bulk <- 
  AverageExpression(sc_data, 
                    features = high_corr$ENSEMBL,
                    group.by = c("cell_type"),
                    assays = "integrated", return.seurat = T)
data.bulk.m<- as.matrix(GetAssayData(data.bulk, slot = "scale.data", 
                                                assay = "integrated"))
rownames(data.bulk.m) <- high_corr$V1
data.bulk.m <- data.bulk.m[, c("ciliated epithelial cell",
                              "secretory cell",
                              "fibroblast",
                              "myofibroblast cell",
                              "smooth muscle cell", 
                              "pericyte",
                              "blood vessel endothelial cell",
                              "endothelial cell of lymphatic vessel",
                              "B cell",
                              "mature NK T cell",
                              "mast cell",
                              "macrophage")]

ha <- 
  HeatmapAnnotation(Cell_Type = factor(colnames(data.bulk.m),
                                       levels = c("ciliated epithelial cell",
                                                  "secretory cell",
                                                  "fibroblast",
                                                  "myofibroblast cell",
                                                  "smooth muscle cell", 
                                                  "pericyte",
                                                  "blood vessel endothelial cell",
                                                  "endothelial cell of lymphatic vessel",
                                                  "B cell",
                                                  "mature NK T cell",
                                                  "mast cell",
                                                  "macrophage")),
                    col = list(Cell_Type = cell_col),
                    simple_anno_size = unit(1, "cm"))
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)

set.seed(1220)
png(file = "visium_DEGs_all_hi_corr_sc_heatmap.png", width = 5000, height = 9000, res = 300)
p <-Heatmap(data.bulk.m, 
            name = "Expression", 
            cluster_columns = F,
            show_column_dend = F,
            show_column_names = F,
            cluster_column_slices = TRUE,
            column_title_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 20),
            column_title = "scRNAseq 12 Major Cell Types",
            col = rev(mapal),
            cluster_rows = clustFunction,
            cluster_row_slices = F,
            row_names_side = "right",
            show_row_dend = T,
            row_names_gp = gpar(fontsize = 6),
            row_title_gp = gpar(fontsize = 12),
            row_title_rot = 0,
            row_title_side = "right",
            top_annotation = ha,
            use_raster = TRUE,
            raster_quality = 4)
p <- draw(p)
dev.off()

#' Extract row order and plot cell type signatures in Visium data
hm_row_ord <-row_order(p)
high_corr_reord <- high_corr$V1[c(hm_row_ord)]%>%as.data.frame()
colnames(high_corr_reord) <- "gene"
high_corr_reord$cell_type <- ""
high_corr_reord$cell_type[c(1:96)]<-"endothelial cell" 
high_corr_reord$cell_type[c(97:245)]<-"Stromal cell" 
high_corr_reord$cell_type[c(246:405)]<-"ciliated epithelial cell" 
high_corr_reord$cell_type[c(406:475)]<-"secretory cell" 
high_corr_reord$cell_type[c(476:551)]<-"immune cell" 
high_corr_reord$cell_type <- factor(high_corr_reord$cell_type, 
                                    levels = c("ciliated epithelial cell",
                                               "secretory cell" ,
                                               "Stromal cell",
                                               "endothelial cell",
                                               "immune cell" ))
write.table(high_corr_reord, 
            file = "Celltype_signatures.csv", 
            sep =",", quote = F, row.names = F, col.names = T)


#' Average expression in our data
load(paste0(rds.dir, "/Variable_features_filt_SCT_log2counts+1_harmony.RData"))

#' Choose the data file 
# data <- data.deg  
data <- data.sub.filt
source(paste(home, "OVisium/manuscript/gitHub/ComplexHeatmap_settings.R", sep = "/")) 
#' Between clusters
all.pat.data.bulk <- 
  AverageExpression(data, features = high_corr_reord$gene,
                    group.by = c("Visium_clusters"),
                    assays = "SCT", return.seurat = T)
all.pat.data.bulk.m <- as.matrix(GetAssayData(all.pat.data.bulk, 
                                              slot = "scale.data", 
                                              assay = "SCT"))

all.pat.ha <- 
  HeatmapAnnotation(Cluster = factor(colnames(all.pat.data.bulk.m), 
                                     levels = c("0_FTE", "1_FTE", "3_Mix", "2_Stroma",
                                                "5_Stroma", "6_Stroma", "7_Stroma",
                                                "8_Stroma", "9_Stroma", "10_Stroma",
                                                "11_Stroma")),
                    col = anno.col,
                    simple_anno_size = unit(1, "cm"))

#' Original order
set.seed(1220)
png(file = "visium_DEGs_all_hi_corr_visium_heatmap.png", width =5000, height = 9000, res = 300) 
p <-Heatmap(all.pat.data.bulk.m, 
            name = "Expression", 
            cluster_columns = F,
            show_column_dend = F,
            show_column_names = F,
            cluster_column_slices = TRUE,
            column_title_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 20),
            column_title = "OVisium 11 Spatial Clusters",
            col = rev(mapal),
            cluster_rows = F,
            row_names_side = "right",
            show_row_dend = T,
            row_names_gp = gpar(fontsize = 6),
            row_title_gp = gpar(fontsize = 20),
            row_title_rot = 0,
            row_title_side = "right",
            row_split = high_corr_reord$cell_type,
            top_annotation = all.pat.ha,
            use_raster = TRUE,
            raster_quality = 4)
p <- draw(p)
dev.off()

#' Pearson clustering
set.seed(1220)
png(file = "visium_DEGs_all_hi_corr_visium_heatmap_cluster_pearson.png", width = 5000, height = 9000, res = 300) 
p <-Heatmap(all.pat.data.bulk.m, 
            name = "Expression", 
            cluster_columns = F,
            show_column_dend = F,
            show_column_names = F,
            cluster_column_slices = TRUE,
            column_title_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 20),
            column_title = "OVisium 11 Spatial Clusters",
            col = rev(mapal),
            cluster_rows = clustFunction,
            cluster_row_slices = F,
            row_names_side = "right",
            show_row_dend = T,
            row_names_gp = gpar(fontsize = 6),
            row_title_gp = gpar(fontsize = 12),
            row_title_rot = 0,
            row_title_side = "right",
            top_annotation = all.pat.ha,
            use_raster = TRUE,
            raster_quality = 4)
p <- draw(p)
dev.off()

