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
file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"
feature.ident <- c("Variable_features_filt")
#' Other features : 
#' "Variable_a500UMI_features","Scaled_features","Scaled_a500UMI_features","All_features"
#' Different list based on the filtering
#' Test methods to identify DE genes
test.ident <- c("MAST")
markers.list.dir <- "FindMarkers/Combined_filtered_markers_v1/rds/Correlation"

input.dir <- paste(deg.dir, 
                   file.name, 
                   cluster.ident, 
                   feature.ident, 
                   test.ident,
                   markers.list.dir,
                   sep = "/")

setwd(input.dir)
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

#' split into DEGs by fold change
common_up <- common[common$avg_log2FC > 0,]
common_de <- common[common$avg_log2FC < 0,]

#' order by foldchange and remove duplicate 
common_up <- common_up[order(-common_up$avg_log2FC),]
common_up <-common_up[!duplicated(common_up$gene),]
visium_DEGs_up <- head(common_up, 500)
visium_DEGs_up <- visium_DEGs_up[order(visium_DEGs_up$id1),]
visium_DEGs_up_gene <- visium_DEGs_up$ENSEMBL
names(visium_DEGs_up_gene) <- visium_DEGs_up$id1
visium_DEGs_up_gene <- visium_DEGs_up_gene[visium_DEGs_up_gene %in% gene_names]


common_de <- common_de[order(-abs(common_de$avg_log2FC)),]
common_de <-common_de[!duplicated(common_de$gene),]
visium_DEGs_de <- head(common_de, 500)
visium_DEGs_de <- visium_DEGs_de[order(visium_DEGs_de$id1),]
visium_DEGs_de_gene <- visium_DEGs_de$ENSEMBL
names(visium_DEGs_de_gene) <- visium_DEGs_de$id1
visium_DEGs_de_gene <- visium_DEGs_de_gene[visium_DEGs_de_gene %in% gene_names]


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
features <- visium_DEGs_up_gene
sc_data <- ScaleData(healthy.tubes, features = features)
sc_data <- RunPCA(sc_data, features = features, npcs = 25)
# elbow_plot<- ElbowPlot(sc_data, ndims = 50)
#ggsave(filename = "visium_DEGs_up_elbow_plot.png", plot = elbow_plot, width = 3.75, height = 3.75, dpi = 300)
sc.features<- rownames(sc_data@reductions[["pca"]]@feature.loadings)
Idents(sc_data) <- "cell_type"
pca_plot <- DimPlot(sc_data, reduction = "pca", pt.size = 0.1, label = F, cols = cell_col)
#legend <- get_legend(pca_plot)
pca_plot<- pca_plot & NoLegend()
ggsave(filename = "visium_DEGs_up_pca_plot.png", plot = pca_plot, width = 3.75, height = 3.75, dpi = 300)
#ggsave(filename = "visium_DEGs_up_legend.png", plot = legend, width = 2, height = 3.75, dpi = 300)

mat<-as.data.frame(sc_data@assays[["integrated"]]@data[sc.features,]) %>% rownames_to_column("ENSEMBL")
mat <- inner_join(mat, features.id, relationship = "many-to-many")
row.names(mat) <- mat$SYMBOL
mat <- mat[-c(1,58815,58816)]

correlations_DEGs_log <- cor(method = "pearson",
                             log2(t(expm1(as.matrix(mat)))+1))

png(file = "visium_DEGs_up_corr_plot.png", width = 9000, height = 9000, res = 300)
heatmapPearson(correlations_DEGs_log)
dev.off()

#' subset DEGs with correlation >0.6
idx <- which(correlations_DEGs_log > 0.5 & lower.tri(correlations_DEGs_log, diag = F), arr.ind = TRUE)
DEGs.cor <-cbind(rownames(idx), colnames(correlations_DEGs_log)[idx[, 2]]) %>% as.data.frame()
DEGs.cor.up <-DEGs.cor[!duplicated(DEGs.cor$V1),]
saveRDS(DEGs.cor.up, file = "High_corr_up_genes.rds")


#' downregulated genes
features <- visium_DEGs_de_gene
sc_data <- ScaleData(healthy.tubes, features = features)
sc_data <- RunPCA(sc_data, features = features, npcs = 25)
#elbow_plot<- ElbowPlot(sc_data, ndims = 50)
#ggsave(filename = "visium_DEGs_de_elbow_plot.png", plot = elbow_plot, width = 3.75, height = 3.75, dpi = 300)
sc.features<- rownames(sc_data@reductions[["pca"]]@feature.loadings)
Idents(sc_data) <- "cell_type"
pca_plot<- DimPlot(sc_data, reduction = "pca", pt.size = 0.1, label = F, cols = cell_col)
#legend <- get_legend(pca_plot)
pca_plot<- pca_plot & NoLegend()
ggsave(filename = "visium_DEGs_de_pca_plot.png", plot = pca_plot, width = 3.75, height = 3.75, dpi = 300)
#ggsave(filename = "visium_DEGs_de_legend.png", plot = legend, width = 2, height = 3.75, dpi = 300)

mat<-as.data.frame(sc_data@assays[["integrated"]]@data[sc.features,]) %>% rownames_to_column("ENSEMBL")
mat <- inner_join(mat, features.id, relationship = "many-to-many")
row.names(mat) <- mat$SYMBOL
mat <- mat[-c(1,58815,58816)]

correlations_DEGs_log <- cor(method = "pearson",
                             log2(t(expm1(as.matrix(mat)))+1))

png(file = "visium_DEGs_de_corr_plot.png", width = 9000, height = 9000, res = 300)
heatmapPearson(correlations_DEGs_log)
dev.off()

#' subset DEGs with correlation >0.6
idx <- which(correlations_DEGs_log > 0.5 & lower.tri(correlations_DEGs_log, diag = F), arr.ind = TRUE)
DEGs.cor <-cbind(rownames(idx), colnames(correlations_DEGs_log)[idx[, 2]]) %>% as.data.frame()
DEGs.cor.down <-DEGs.cor[!duplicated(DEGs.cor$V1),]
saveRDS(DEGs.cor.up, file = "High_corr_down_genes.rds")

#' plot those highly correlated genes on heatmap
Idents(healthy.tubes) <- "cell_type"
high_corr <- rbind(DEGs.cor.up, DEGs.cor.down)
high_corr <-high_corr[!duplicated(high_corr$V1),]
high_corr <- inner_join(high_corr, features.id, by = c("V1" = "SYMBOL"), 
                     relationship = "many-to-many")
high_corr_genes <- high_corr$ENSEMBL

#' perform PCA 
features <- high_corr_genes
sc_data <- ScaleData(healthy.tubes, features = features)
sc_data <- RunPCA(sc_data, features = features, npcs = 25)
elbow_plot<- ElbowPlot(sc_data, ndims = 25)
ggsave(filename = "visium_DEGs_hi_corr_elbow_plot.png", plot = elbow_plot, width = 3.75, height = 3.75, dpi = 300)
pca_plot<- DimPlot(sc_data, reduction = "pca", pt.size = 0.1, label = F, cols = cell_col)
legend <- get_legend(pca_plot)
pca_plot<- pca_plot & NoLegend()
ggsave(filename = "visium_DEGs_hi_corr_pca_plot.png", plot = pca_plot, width = 3.75, height = 3.75, dpi = 300)
ggsave(filename = "cell_type_legend.png", plot = legend, width = 2, height = 3.75, dpi = 300)

data.bulk <- 
  AverageExpression(sc_data, 
                    features = features,
                    group.by = c("cell_type"),
                    assays = "integrated", return.seurat = T)
data.bulk.m<- as.matrix(GetAssayData(data.bulk, slot = "scale.data", 
                                                assay = "integrated"))
data.bulk.m <- data.bulk.m %>% as.data.frame() %>% rownames_to_column("gene")
data.bulk.m <- inner_join(data.bulk.m, features.id, by = c("gene" = "ENSEMBL"), 
                        relationship = "many-to-many")
data.bulk.m <- column_to_rownames(data.bulk.m, var = "SYMBOL")
data.bulk.m <- data.bulk.m[-c(1,14,15)] %>% as.matrix()
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
                    col = list(Cell_Type = cell_col))
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256)

png(file = "visium_DEGs_hi_corr_heatmap.png", width = 4000, height = 8000, res = 300)
p <-Heatmap(data.bulk.m, 
            name = "Expression", 
            cluster_columns = F,
            show_column_dend = F,
            show_column_names = F,
            cluster_column_slices = TRUE,
            column_title_gp = gpar(fontsize = 12),
            column_names_gp = gpar(fontsize = 20),
            column_title = "scRNAseq 12 Major Cell Types",
            col = rev(mapal),
            cluster_rows = clustFunction,
            cluster_row_slices = F,
            row_names_side = "right",
            show_row_dend = T,
            row_names_gp = gpar(fontsize = 8),
            row_title_gp = gpar(fontsize = 12),
            row_title_rot = 0,
            row_title_side = "right",
            top_annotation = ha,
            use_raster = TRUE,
            raster_quality = 4)
p <- draw(p)
dev.off()
