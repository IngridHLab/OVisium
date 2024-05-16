#' Perform different gene expression analysis using Seurat functions:
#' FindAllMarkers(): Gene expression markers for all identity classes;
#' FindMarker(): Gene expression markers of identity classes 
#' features 

library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/17_geneAnnotations.R", sep = "/"))

file.name <- "OVisium_SCT_merged" 
cluster.ident <- "harmony_SCT_res_0.6"
setwd(rds.dir)
#' load("Variable_features_filt_SCT_harmony_data.RData")
load("Variable_features_filt_SCT_log2counts+1_harmony.RData")
#' load("Variable_features_filt_SCT_log1pcounts_harmony.RData")


#' prepare the data slot for DEA using FindMarkers
#' data.sub@assays[["SCT"]]@data <- Matrix(log1p(2^(as.matrix(data.sub@assays[["SCT"]]@data))-1), sparse = T)
#' data.sub <- ScaleData(data.sub, assay = "SCT")
#' for scale.data analysis, no filtering will be applied on the log2FC
#' subset only fimbrial samples
Idents(data.sub.filt) <- "Tissue_origin"
data.sub.filt <- SubsetSTData(data.sub.filt, idents = c("Fimbrial"))
data.sub.filt@assays[["SCT"]]@scale.data <- as.matrix(data.sub.filt@assays[["SCT"]]@data)
data.sub.filt <- ScaleData(data.sub.filt)

#' In order to test all gene, we set min.diff.pct = 0.1, logfc.threshold=-Inf
#' choose different mutation to perform comparison
id1 = "BRCA1" 
id2 = "BRCA2"
data <- data.sub.filt
features <- vfeatures.filt

#' patient no. 4 developed HGSC 6 years after the salpingectomy

  out.dir <- paste(deg.dir, file.name, cluster.ident, 
                   "Variable_features_filt", 
                   "MAST", 
                   "FindMarkers", 
                   "Mutation",
                   sep = "/") 
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
      
  #' Find all markers between fimbrial vs proximal in PCA-Harmony
  #' log2FC threshold default 0.1
      Idents(data.sub.filt) <- "mutation"
      markers <- FindMarkers(data.sub.filt, 
                               assay = "SCT",
                               slot = "scale.data",
                               fc.name = "avg_log2FC",
                               features = features,
                               ident.1 = id1, 
                               ident.2 = id2,
                               min.pct = -Inf,
                               logfc.threshold = -Inf,
                               min.diff.pct = 0.1,
                               test.use = "MAST",
                               recorrect_umi = F) %>% 
        rownames_to_column() %>% 
        dplyr::rename("gene" = 1) %>% 
        mutate(rank=-log10(p_val_adj+2.225074e-308)*avg_log2FC, 
               id1 = id1, 
               id2 = id2) %>% 
        dplyr::filter(abs(rank) > -log10(10E-5)*-Inf) %>%
        dplyr::relocate(id1, id2, .after = gene) 
        
       markers <- markers[order(-markers$avg_log2FC),] %>%
         left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                   by = c("gene" = "gene_name"))
       
       markers.top <- markers %>% 
         dplyr::filter(abs(avg_log2FC) >= 0.25) %>%
         dplyr::slice_max(order_by = abs(rank), n = 20)
       markers.top <- markers.top[order(-markers.top$avg_log2FC),] %>%
         left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                   by = c("gene" = "gene_name"))
      
       write.table(markers, 
                    file = "BRCA1_vs_BRCA2.csv", 
                    sep =",", quote = F, row.names = F, col.names = T)
       write.table(markers.top, 
                   file = "BRCA1_vs_BRCA2_top20.csv", 
                   sep =",", quote = F, row.names = F, col.names = T) 
  
        #' plot de gene using volcano plot
          #' enhancevolcano package 
          #' https://github.com/kevinblighe/EnhancedVolcano
            png(file = "BRCA_1vs2_volcano.png", 
                  width = 3500, height = 3500, res = 400)
            print(
              EnhancedVolcano(markers,
                              lab = markers$gene,
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              title = "BRCA1 vs BRCA2",
                              subtitle = "Fimbrial",
                              pCutoff = 0.01,
                              FCcutoff = log2(1.5),                              
                              pointSize = 3.0,
                              labSize = 6.0,
                              legendLabels=c('Not sig.',
                                             'avg_Log2(FC)',
                                             'p_adj_val',
                                             'p_adj_val & avg_Log2(FC)')) +
                labs(x = "avg_Log2(Fold Change)", y = "-log10 (Adjust P_value)") +
                facet_wrap(. ~ id1, ncol = 3)
              )
            dev.off()

out.dir <- paste(deg.dir, file.name, cluster.ident, 
                 "Variable_features_filt", sep = "/") 
saveRDS(markers, file="BRCA.DEA.list.rds")
saveRDS(markers.top, file="BRCA.DEA.top20.list.rds")

#' complexheatmap with all genes
#' Original Seurat Heatmap on all clusters
#' Compare clusters in proximal vs fimbrial across all 4 patients

clust.col <- 
  c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#F0027F","#BF5B17",
             "#666666","#1B9E77","#D95F02","#386CB0","black")
             
names(clust.col) <-c("0_FTE","1_FTE", "3_Mix", "2_Stroma",
                     "5_Stroma", "6_Stroma", "7_Stroma",
                     "8_Stroma", "9_Stroma", "10_Stroma",
                     "11_Stroma")
patient.col <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(12) 
names(patient.col) <- factor(1:12)
sample.col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(14) 
names(sample.col) <- c("PO2", "PO6", "PO7_A1", "PO9_B1", 
                       "PO13", "PO13a", "PO20", "PO20a", "PO28", "PO35", 
                       "PO37", "PO40", "PO41", "PO45")
batch.col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(5) 
names(batch.col) <- c("V10S21-048" ,"V10S21-050","V10T06-031","V10T06-110","V10U22-108") 
variant.col <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(9) 
names(variant.col) <- c("c.3048_3052dup" ,"c.1082_1092del","c.68_69del",
                         "c.1796_1800del","c.1687C>T" ,"c.2830A>T","c.7007G>A",
                        "c.81-1588_134+1725del3367","c.3626del") 

anno.col = list(
  Mutation = c(BRCA1 = "red", BRCA2 = "blue"),
  Variant = variant.col,
  Cluster = clust.col,
  Sample = sample.col,
  Patient = patient.col,
  Tissue = c(Fimbrial = "#FC2947", Otherside = "#FE6244", Proximal = "#FFDEB9"),
  Age = c("35+" = "#EBE76C", "40+" = "#F0B86E", "45+" = "#ED7B7B", "50+" = "#836096"),
  Batch = batch.col
)

puryel = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00"))
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

#' Seurat assay data and annotation for individual spots
Idents(data.sub.filt) <- "Visium_clusters"
cellInfo <- data.frame(Visium_clusters=Idents(data.sub.filt))
cellInfo$Visium_clusters <- factor(cellInfo$Visium_clusters, 
                            levels = c("0_FTE", "1_FTE", "3_Mix", 
                                       "2_Stroma", "5_Stroma", 
                                       "6_Stroma", "7_Stroma",
                                       "8_Stroma", "9_Stroma", 
                                       "10_Stroma", "11_Stroma"))
cellInfo = cellInfo %>% dplyr::arrange(Visium_clusters)

ha = HeatmapAnnotation(Mutation = data.sub.filt@meta.data$mutation,
                       Variant = data.sub.filt@meta.data$variant,
                       Cluster = data.sub.filt@meta.data$Visium_clusters, 
                       Sample = data.sub.filt@meta.data$library,
                       Patient = data.sub.filt@meta.data$patient,
                       Tissue = data.sub.filt@meta.data$origin,
                       Age = data.sub.filt@meta.data$age,
                       Batch = data.sub.filt@meta.data$slide,
                       col = anno.col,
                       annotation_name_gp= gpar(fontsize = 18),
                       annotation_legend_param = list(
                         Mutation = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14)),
                         Variant = list(
                           title_gp = gpar(fontsize = 14, 
                                           fontface = "bold"), 
                           labels_gp = gpar(fontsize = 14)),
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

#' plot regular seurat heatmap
png(file = "BRCA1_vs_BRCA2_top20_heatmap.png", 
    width = 3000, height = 3000, res = 250) 
print(
  DoHeatmap(data.sub.filt, 
            assay = "SCT",
            features = markers.top$gene,
            size = 6, 
            angle = 45) +
    theme(text = element_text(size = 10)) +
    theme(axis.text = element_text(size = 12)) + 
    guides(color="none") +
    labs(fill = "Expression")
)
dev.off() 

#' plot complex heatmap
mat<- data.sub.filt[["SCT"]]@data[markers.top$gene, ] %>% as.matrix()
## scale the rows
mat<- t(scale(t(mat)))
dend1 <- cluster_between_groups(mat, data.sub.filt@meta.data$mutation)

set.seed(1220)
png(file = "BRCA1.vs.BRCA2.top20.cheatmap.png", 
    width = 5000 , height = 5000, res = 300) 
p <-Heatmap(mat, 
            name = "Expression", 
            column_split = 2,
            cluster_columns = dend1,
            column_order = row.names(cellInfo),
            show_column_dend = F,
            show_column_names = F,
            column_title_gp = gpar(fontsize = 12),
            column_title = "BRCA1 vs BRCA2",
            col = puryel,
            cluster_rows = F,
            row_names_side = "left",
            show_row_dend = F,
            row_names_gp = gpar(fontsize = 18),
            row_title_gp = gpar(fontsize = 12),
            row_title_rot = 0,
            row_title_side = "right",
            cluster_row_slices = F,
            top_annotation = ha,
            use_raster = TRUE,
            raster_resize_mat = TRUE,
            raster_quality = 4)
p <- draw(p)
dev.off()


#' bulk analysis of BRCA status and variant
mutation.data.bulk <- 
  AverageExpression(data.sub.filt, 
                    features = markers$gene,
                    group.by = c("mutation","variant"),
                    assays = "SCT", return.seurat = T)
mutation.data.bulk.m <- as.matrix(GetAssayData(mutation.data.bulk, 
                                                   slot = "scale.data", 
                                                   assay = "SCT"))
mutation.ha <- 
  HeatmapAnnotation(Mutation = gsub("_.*", "", colnames(mutation.data.bulk.m)),
                    Variant = str_sub(colnames(mutation.data.bulk.m), start = 7),
                    col = anno.col)

set.seed(1220)
png(file = "BRCA1_vs_BRCA2_141_mutation_variant_bulk.png", width = 4000, height = 7000, res = 400) 
p <-Heatmap(mutation.data.bulk.m, 
            name = "Expression", 
            cluster_columns = T,
            clustering_distance_columns = "euclidean",
            clustering_method_columns = "average",
            show_column_dend = T,
            show_column_names = T,
            cluster_column_slices = TRUE,
            column_title_gp = gpar(fontsize = 20),
            column_title = "BRCA1 vs BRCA2",
            col = rev(mapal),
            cluster_rows = F,
            row_names_side = "left",
            show_row_dend = F,
            row_names_gp = gpar(fontsize = 8),
            row_title_gp = gpar(fontsize = 12),
            row_title_rot = 90,
            row_title_side = "right", 
            top_annotation = mutation.ha,
            use_raster = TRUE,
            raster_quality = 4)
p <- draw(p)
dev.off()

#' bulk with BRCA status and patient
mutation.data.bulk <- 
  AverageExpression(data.sub.filt, 
                    features = markers$gene,
                    group.by = c("mutation","patient"),
                    assays = "SCT", return.seurat = T)
mutation.data.bulk.m <- as.matrix(GetAssayData(mutation.data.bulk, 
                                               slot = "scale.data", 
                                               assay = "SCT"))
mutation.ha <- 
  HeatmapAnnotation(Mutation = gsub("_.*", "", colnames(mutation.data.bulk.m)),
                    Patient =  gsub(".*_", "", colnames(mutation.data.bulk.m)),
                    col = anno.col)

set.seed(1220)
png(file = "BRCA1_vs_BRCA2_141_mutation_patient_bulk.png", width = 4000, height = 7000, res = 400) 
p <-Heatmap(mutation.data.bulk.m, 
            name = "Expression", 
            cluster_columns = T,
            clustering_distance_columns = "euclidean",
            clustering_method_columns = "average",
            show_column_dend = T,
            show_column_names = T,
            cluster_column_slices = TRUE,
            column_title_gp = gpar(fontsize = 20),
            column_title = "BRCA1 vs BRCA2",
            col = rev(mapal),
            cluster_rows = F,
            row_names_side = "left",
            show_row_dend = F,
            row_names_gp = gpar(fontsize = 8),
            row_title_gp = gpar(fontsize = 12),
            row_title_rot = 90,
            row_title_side = "right", 
            top_annotation = mutation.ha,
            use_raster = TRUE,
            raster_quality = 4)
p <- draw(p)
dev.off()

