#' Perform different gene expression analysis using Seurat functions:
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
data.sub.filt@assays[["SCT"]]@scale.data <- as.matrix(data.sub.filt@assays[["SCT"]]@data)
data.sub.filt <- ScaleData(data.sub.filt)

#' In order to test all gene, we set min.diff.pct = 0, logfc.threshold=-Inf
#' only compare the fimbrial and proximal from the same side (excluded the
#' proximal tissue from the other side) of patient no. 1, 3, 6, 10
de <- list()
de.top <- list()

#' choose different tissue origin to perform comparison
id1 = "Fimbrial" 
id2 = "Proximal"
data <- data.sub.filt
features <- vfeatures.filt

#' patient no. 4 developed HGSC 6 years after the salpingectomy
for (p in c(1,3,6,10)) {
  Idents(data) <- "patient"
  pat.data <- SubsetSTData(data, idents = p)
  patient <- paste("Patient", p, sep = "_")
  out.dir <- paste(deg.dir, file.name, cluster.ident, 
                   "Variable_features_filt", 
                   "MAST", 
                   "FindMarkers", 
                   "Patient_Pair_core",
                   sep = "/") 
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  
      
  #' Find all markers between fimbrial vs proximal in PCA-Harmony
  #' log2FC threshold default 0.1
      Idents(pat.data) <- "Tissue_origin"
      f.markers <- FindMarkers(pat.data, 
                               assay = "SCT",
                               slot = "scale.data",
                               fc.name = "avg_log2FC",
                               features = features,
                               ident.1 = id1, 
                               ident.2 = id2,
                               min.pct = -Inf,
                               logfc.threshold = -Inf,
                               min.diff.pct = -Inf,
                               test.use = "MAST",
                               recorrect_umi = F) %>% 
        rownames_to_column() %>% 
        dplyr::rename("gene" = 1) %>% 
        mutate(rank=-log10(p_val_adj+2.225074e-308)*avg_log2FC, 
               id1 = id1, 
               id2 = id2,
               patient = patient) %>%
        dplyr::relocate(patient, id1, id2, .after = gene) 
        
        f.markers.top <- f.markers %>% 
          dplyr::filter(avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75)) %>% 
          dplyr::group_by(patient) %>%
          dplyr::slice_max(order_by = abs(rank), n = 20)
        
      
        write.table(f.markers, 
                    file = paste(patient, "fimbrial_vs_proximal.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        write.table(f.markers.top, 
                    file = paste(patient, p, "fimbrial_vs_proximal.top.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
  
        #' plot de gene using volcano plot
         de[[paste("Pat", p, sep = " ")]] <- f.markers

        
        for (n in seq_along(names(de))) {
          #' enhancevolcano package 
          #' https://github.com/kevinblighe/EnhancedVolcano
            if ( n != 2 ) {
              png(file = paste0(patient, "_volcano_", n, ".png"), 
                  width = 3500, height = 3500, res = 400)
            } else {
              png(file = paste0(patient, "_volcano_", n, ".png"), 
                  width = 8000, height = 8000, res = 500)
            }
            
            print(
              EnhancedVolcano(de[[n]],
                              lab = de[[n]]$gene,
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              title = patient,
                              subtitle = names(de[n]),
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
        }
        
        #' save top de for violin plot later across all patients
        de.top[[paste("Pat", p, sep = " ")]] <- f.markers.top
        
        
        
    }


saveRDS(de, file="Patients.DEA.list.rds")
saveRDS(de.top, file="Patients.DEA.top.list.rds")

#' complexheatmap with top 20 genes from each patients
#' Original Seurat Heatmap on all clusters
#' Compare clusters in proximal vs fimbrial across all 4 patients
Idents(data) <- "patient"
paired.pat.data <- SubsetSTData(data, idents = c(1,3,6,10))

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
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

#' Seurat assay data and annotation for individual spots
cellInfo <- data.frame(Visium_clusters=Idents(paired.pat.data))
cellInfo$Visium_clusters <- factor(cellInfo$Visium_clusters, 
                                   levels = c("0_FTE", "1_FTE", "3_Mix", 
                                              "2_Stroma", "5_Stroma", 
                                              "6_Stroma", "7_Stroma",
                                              "8_Stroma", "9_Stroma", 
                                              "10_Stroma", "11_Stroma"))
cellInfo = cellInfo %>% dplyr::arrange(Visium_clusters)
ha = HeatmapAnnotation(Cluster = paired.pat.data@meta.data$Visium_clusters, 
                       Sample = paired.pat.data@meta.data$library,
                       Patient = paired.pat.data@meta.data$patient,
                       Tissue = paired.pat.data@meta.data$origin,
                       Age = paired.pat.data@meta.data$age,
                       Batch = paired.pat.data@meta.data$slide,
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

#' bulk analysis of proximal and fimbrial between all 4 patients
DEA.top.list <- readRDS("Patients.DEA.top.list.rds")
markers<-bind_rows(DEA.top.list) 

#' remove duplicates
markers <- 
  markers[order(-markers$avg_log2FC),]
markers <-
  markers[!duplicated(markers$gene),]
markers <- 
  markers[order(markers$patient,-markers$avg_log2FC),]
markers$patient <- factor(markers$patient, levels = c("Patient_1", "Patient_3", 
                                                      "Patient_6", "Patient_10"))

paried.pat.data.bulk <- 
  AverageExpression(paired.pat.data, 
                    features = markers$gene,
                    group.by = c("patient","Tissue_origin"),
                    assays = "SCT", return.seurat = T)
paried.pat.data.bulk.m <- as.matrix(GetAssayData(paried.pat.data.bulk, 
                                                   slot = "scale.data", 
                                                   assay = "SCT"))
paired.pat.ha <- 
  HeatmapAnnotation(Tissue = gsub(".*_", "", colnames(paried.pat.data.bulk.m)),
                    Patient = gsub("_.*", "", colnames(paried.pat.data.bulk.m)),
                    col = anno.col)

set.seed(1220)
png(file = paste("Paired_Cases_cHeatmap", nrow(paried.pat.data.bulk.m), 
                 "biomarkers_tissue_bulk.png", sep = "_"), 
    width = 3000 , height = 6000, res = 400) 
p <-Heatmap(paried.pat.data.bulk.m, 
            name = "Expression", 
            cluster_columns = T,
            show_column_dend = T,
            show_column_names = T,
            cluster_column_slices = TRUE,
            column_title_gp = gpar(fontsize = 20),
            column_title = "Combined 4 Cases",
            col = rev(mapal),
            cluster_rows = F,
            row_names_side = "left",
            show_row_dend = F,
            row_names_gp = gpar(fontsize = 8),
            row_title_gp = gpar(fontsize = 12),
            row_title_rot = 90,
            row_title_side = "right", 
            row_split = markers$patient,
            top_annotation = paired.pat.ha,
            use_raster = TRUE,
            raster_quality = 4)
p <- draw(p)
dev.off()
