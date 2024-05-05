#' Perform different gene expression analysis using Seurat functions:
#' FindAllMarkers(): Gene expression markers for all identity classes;
#' FindMarker(): Gene expression markers of identity classes 
#' features 
#' https://divingintogeneticsandgenomics.com/post/do-you-really-understand-log2fold-change-in-single-cell-rnaseq-data/
#' 

library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/17_geneAnnotations.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/DoMultiBarHeatmap.R", sep = "/"))

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
Idents(data.sub.filt) <- "Visium_clusters"
data.sub.filt <- SubsetSTData(data.sub.filt, idents = c("2_Stroma",
                                                 "5_Stroma", "6_Stroma", "7_Stroma",
                                                 "8_Stroma", "9_Stroma", "10_Stroma",
                                                 "11_Stroma"))
data.sub.filt@assays[["SCT"]]@scale.data <- as.matrix(data.sub.filt@assays[["SCT"]]@data)
data.sub.filt <- ScaleData(data.sub.filt)

feature.list <- list(vfeatures.filt) # vfeatures 
names(feature.list) <- c("Variable_features_filt") # "Variable_features"
test <- c("MAST") # "wilcox", "MAST"

FindMarkers_filter_combined_Stroma <- function(data, feature.list, slot, test, out.name) {
  
  #' Original Seurat Heatmap on all clusters
  clust.col <- 
    c("#FFFF99","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02","#386CB0","black")
               
  names(clust.col) <-c("2_Stroma",
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
                                     levels = c("2_Stroma", "5_Stroma", 
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
  
  #' Find all markers of all graph-based clusters in PCA-Harmony
  for(f in seq_along(feature.list)) {
    for(t in seq_along(test)) {
      features <- feature.list[[f]]
      out.dir <- paste(deg.dir, file.name, cluster.ident, 
                       names(feature.list[f]), test[t], 
                       "FindMarkers", 
                       "Combined_filtered_markers_v1",
                       paste(Sys.Date(), out.name, sep = "_"),
                       sep = "/") 
      dir.create(out.dir, recursive = T)
      setwd(out.dir)
      Idents(data) <- "Visium_clusters"
      
      #' Find differential expressed markers between 8 clusters
      #' Stroma clusters 
      #' min.diff.pct at 0.25 has no gene found, set 10 instead
      combined.markers.filt.list <- list()
      for (x in seq_along(levels(data))) {
        gene.list <- list()
        for (y in seq_along(1:8)) {
          id1=levels(data)[[x]]
          id2=levels(data)[[y]]
          if (gsub(".*_", "", id1) != gsub(".*_", "", id2)) {
            log2fc <- log2(2.5)
            diff <- 0.25
          } else { 
            log2fc <- log2(1.5)
            diff <- 0.1
          }
          if (id1 != id2) {
            if (y > 8) break 
            try(gene.list[[id2]] <- FindMarkers(data, 
                                                assay = "SCT",
                                                slot = slot,
                                                fc.name = "avg_log2FC",
                                                features = features,
                                                ident.1 = id1, 
                                                ident.2 = id2,
                                                min.pct = -Inf,
                                                logfc.threshold = -Inf,
                                                min.diff.pct = diff,
                                                test.use = test[t],
                                                recorrect_umi = F) %>% 
                  rownames_to_column() %>% 
                  dplyr::rename("gene" = 1) %>% 
                  mutate(rank=-log10(p_val_adj+2.225074e-308)*avg_log2FC, 
                         id1 = id1, 
                         id2 = id2) %>% 
                  dplyr::filter(abs(avg_log2FC) > log2fc,
                                abs(rank) > -log10(10E-5)*log2fc) %>%
                  dplyr::relocate(id1, id2, .after = gene) 
            )
            y <- y + 1
            
          } else {
            next
          }
        }
        
        sub.out.dir <- c("rds","csv","heatmap","cheatmap")
        for (name in sub.out.dir) {
          path <- paste0(out.dir, "/", name)
          dir.create(path)
        }
        
        gene.list[sapply(gene.list, is.null)] <- NULL
        saveRDS(gene.list, file = paste0("rds/Cluster.",levels(data)[[x]],".list.rds"))
        
        all.markers.filt <- bind_rows(gene.list) 
        all.markers.filt$id1 <- factor(all.markers.filt$id1, 
                                       levels = c("2_Stroma", "5_Stroma", 
                                                  "6_Stroma", "7_Stroma",
                                                  "8_Stroma", "9_Stroma", 
                                                  "10_Stroma", "11_Stroma"))
        all.markers.filt$id2 <- factor(all.markers.filt$id2, 
                                       levels = c("2_Stroma", "5_Stroma", 
                                                  "6_Stroma", "7_Stroma",
                                                  "8_Stroma", "9_Stroma", 
                                                  "10_Stroma", "11_Stroma"))
        
        all.markers.filt <- all.markers.filt[order(-all.markers.filt$avg_log2FC),]
        write.table(all.markers.filt, 
                    file = paste0("csv/Cluster.",id1,".dup.csv"),
                    sep =",", quote = F, row.names = F, col.names = T)
        
        combined.markers.filt <-all.markers.filt[!duplicated(all.markers.filt$gene),] %>%
          left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                    by = c("gene" = "gene_name"))
        write.table(combined.markers.filt, 
                    file = paste0("csv/Cluster.",id1,".once.csv"),
                    sep =",", quote = F, row.names = F, col.names = T)
        
        combined.markers.filt <- combined.markers.filt[order(combined.markers.filt$id2),]
        combined.markers.filt.list[[id1]] <- combined.markers.filt
        
        combined.markers.filt.top20 <- combined.markers.filt %>% 
          dplyr::group_by(id2) %>%
          dplyr::slice_max(order_by = rank, n = 20)
        
        write.table(combined.markers.filt.top20, 
                    file = paste0("csv/Cluster.",id1,".top20.csv"),
                    sep =",", quote = F, row.names = F, col.names = T)
        
        #' plot regular seurat heatmap
        png(file = paste0("heatmap/Cluster.",id1,".top20.heatmap.png"), 
            width = 6000, height = 8000, res = 300) 
        print(
          DoHeatmap(data, 
                    assay = "SCT",
                    features = combined.markers.filt.top20$gene,
                    size = 6, 
                    angle = 45) +
            theme(text = element_text(size = 10)) +
            theme(axis.text = element_text(size = 12)) + 
            guides(color="none") +
            labs(fill = "Expression")
        )
        dev.off() 
        
        #' plot complex heatmap
        mat<- data[["SCT"]]@data[combined.markers.filt.top20$gene, ] %>% as.matrix()
        ## scale the rows
        mat<- t(scale(t(mat)))
        
        set.seed(1220)
        png(file = paste0("cheatmap/Cluster.",id1,".top20.cheatmap.png"), 
            width = 6000 , height = 8000, res = 400) 
        p <-Heatmap(mat, 
                    name = "Expression", 
                    cluster_columns = F,
                    column_order = row.names(cellInfo),
                    show_column_dend = F,
                    show_column_names = F,
                    column_title_gp = gpar(fontsize = 12),
                    column_title = paste("FindMarkers", combined.markers.filt$id1[1], "Top 20",sep = " "),
                    col = puryel,
                    cluster_rows = F,
                    row_names_side = "left",
                    show_row_dend = F,
                    row_names_gp = gpar(fontsize = 8),
                    row_title_gp = gpar(fontsize = 12),
                    row_title_rot = 0,
                    row_title_side = "right", 
                    row_split = combined.markers.filt.top20$id2,
                    cluster_row_slices = F,
                    top_annotation = ha,
                    use_raster = TRUE,
                    raster_resize_mat = TRUE,
                    raster_quality = 4)
        p <- draw(p)
        dev.off()
      }
      
      combined.markers.filt.list[sapply(combined.markers.filt.list, is.null)] <- NULL
      saveRDS(combined.markers.filt.list, file = paste0("rds/Cluster.all.list.rds"))
      
      combined.markers.filt.merged <- 
        bind_rows(combined.markers.filt.list[c("2_Stroma", 
                                               "5_Stroma", "6_Stroma", "7_Stroma",
                                               "8_Stroma", "9_Stroma", "10_Stroma", 
                                               "11_Stroma")])
      write.table(combined.markers.filt.merged, 
                  file = paste0("csv/Cluster.all.dup.csv"),
                  sep =",", quote = F, row.names = F, col.names = T)
      
      combined.markers.filt.merged <- 
        combined.markers.filt.merged[order(-combined.markers.filt.merged$avg_log2FC),]
      combined.markers.filt.merged <-
        combined.markers.filt.merged[!duplicated(combined.markers.filt.merged$gene),]
      combined.markers.filt.merged <- 
        combined.markers.filt.merged[order(combined.markers.filt.merged$id1,
                                           -combined.markers.filt.merged$avg_log2FC),]
      
      write.table(combined.markers.filt.merged, 
                  file = paste0("csv/Cluster.all.once.csv"),
                  sep =",", quote = F, row.names = F, col.names = T)
      
      combined.markers.filt.merged.top20 <- combined.markers.filt.merged %>% 
        dplyr::group_by(id1) %>%
        dplyr::filter(avg_log2FC > log2(2)) %>%
        dplyr::slice_max(order_by = rank, n = 20) 
      
      write.table(combined.markers.filt.merged.top20, 
                  file = paste0("csv/Cluster.all.once.top20.csv"),
                  sep =",", quote = F, row.names = F, col.names = T)
     
       #' plot regular seurat heatmap
      png(file = paste0("heatmap/Cluster.all.top20.heatmap.png"), 
          width = 6000, height = 8000, res = 300) 
      print(
        DoHeatmap(data, 
                  assay = "SCT",
                  features = combined.markers.filt.merged.top20$gene,
                  size = 6, 
                  angle = 45) +
          theme(text = element_text(size = 10)) +
          theme(axis.text = element_text(size = 12)) + 
          guides(color="none") +
          labs(fill = "Expression")
      )
      dev.off() 
      
      #' plot complex heatmap
      mat.merged <- data[["SCT"]]@data[combined.markers.filt.merged.top20$gene, ] %>% as.matrix()
      ## scale the rows
      mat.merged <- t(scale(t(mat.merged)))

      set.seed(1220)
      png(file = paste0("cheatmap/Cluster.all.top20.cheatmap.png"), 
          width = 8000 , height = 8000, res = 450) 
      p <-Heatmap(mat.merged, 
                  name = "Expression", 
                  cluster_columns = F,
                  column_order = row.names(cellInfo),
                  show_column_dend = F,
                  show_column_names = F,
                  column_title_gp = gpar(fontsize = 12),
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
                  use_raster = TRUE,
                  raster_resize_mat = TRUE,
                  raster_quality = 4)
      p <- draw(p)
      dev.off()
    }
  }

}

FindAllMarkers_filter_combined_Stroma <- function(data, slot, feature.list, test, out.name) {
  
  #' Original Seurat Heatmap on all clusters
  clust.col <- 
    c("#FFFF99","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02","#386CB0","black")
               
  names(clust.col) <-c("2_Stroma",
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
  
  blured <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(256))
  puryel = circlize::colorRamp2(c(-2, 0, 2), c("#FF00FF", "black", "#FFFF00"))

  #' Seurat assay data and annotation for individual spots
  cellInfo <- data.frame(Visium_clusters=Idents(data))
  cellInfo$Visium_clusters <- factor(cellInfo$Visium_clusters, 
                                     levels = c("2_Stroma", "5_Stroma", 
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
                         col = anno.col)
  
  #' Find all markers of all graph-based clusters in PCA-Harmony
  for(f in seq_along(feature.list)) {
    for(t in seq_along(test)) {
      features <- feature.list[[f]]
      out.dir <- paste(deg.dir, file.name, cluster.ident, 
                       names(feature.list[f]), test[t], 
                       "FindAllMarkers", 
                       "Combined_filtered_markers_v1",
                       paste(Sys.Date(), out.name, sep = "_"),
                       sep = "/") 
      dir.create(out.dir, recursive = T)
      setwd(out.dir)
      Idents(data) <- "Visium_clusters"
      
      #' identify cluster markers
      log2fc <- log2(1.5) #'log2(1.5) or 0.25
      diff <- 0.25
      all.markers <- FindAllMarkers(data, 
                                    assay = "SCT", 
                                    slot = slot,
                                    fc.name = "avg_log2FC",
                                    features = features,
                                    min.diff.pct = diff,
                                    logfc.threshold = -Inf,
                                    test.use = test[t],
                                    recorrect_umi = F) %>% 
        mutate(rank=-log10(p_val_adj+2.225074e-308)*avg_log2FC) %>%
        left_join(y = unique(annotations[, c("gene_name", "description")]), 
                  by = c("gene" = "gene_name")) %>% 
        dplyr::filter(avg_log2FC > log2fc, 
                      abs(rank) > -log10(10E-5)*log2fc)
        
      all.markers$cluster <- factor(all.markers$cluster, 
                                    levels = c("2_Stroma", "5_Stroma", 
                                               "6_Stroma", "7_Stroma",
                                               "8_Stroma", "9_Stroma", 
                                               "10_Stroma", "11_Stroma"))
     
      # top 20 positive up-regulated genes
      all.markers.top20 <- all.markers %>% 
        dplyr::group_by(cluster) %>%
        dplyr::slice_max(order_by = rank, n = 20) 
      
      sub.out.dir <- c("csv","heatmap","cheatmap")
      for (name in sub.out.dir) {
        path <- paste0(out.dir, "/", name)
        dir.create(path)
      }
        
      write.table(all.markers, 
                  file = "csv/All.clusters.markers.csv",
                  sep =",", quote = F, row.names = F, col.names = T)
      
      write.table(all.markers.top20, 
                  file = "csv/All.clusters.markers.top20.csv",
                  sep =",", quote = F, row.names = F, col.names = T)

          
        #' plot regular seurat heatmap
        png(file = paste0("heatmap/All.clusters.markers.top20.heatmap.png"), 
            width = 6000, height = 8000, res = 300) 
        print(
          DoHeatmap(data, 
                    assay = "SCT",
                    features = all.markers.top20$gene,
                    size = 6, 
                    angle = 45) +
            theme(text = element_text(size = 10)) +
            theme(axis.text = element_text(size = 12)) + 
            guides(color="none") +
            labs(fill = "Expression")
        )
        dev.off()
        
      
        #' plot complex heatmap
        mat <- data[["SCT"]]@data[all.markers.top20$gene, ] %>% as.matrix()
        ## scale the rows
        mat <- t(scale(t(mat)))
        #' make the color scale
        #scale <- quantile(mat, c(0.1, 0.95)) %>% unname() %>% abs() %>% round(digits = 0)
        #puryel = circlize::colorRamp2(c(-scale[1], 0, scale[2]), c("#FF00FF", "black", "#FFFF00"))
        
        
        set.seed(1220)
        png(file = paste0("cheatmap/All.clusters.markers.top20.cheatmap.png"), 
            width = 6000 , height = 8000, res = 300) 
        p <-Heatmap(mat, 
                    name = "Expression", 
                    cluster_columns = F,
                    column_order = row.names(cellInfo),
                    show_column_dend = F,
                    show_column_names = F,
                    column_title_gp = gpar(fontsize = 12),
                    column_title = "FindAllMarkers All Clusters Top 20",
                    col = puryel,
                    cluster_rows = F,
                    row_names_side = "left",
                    show_row_dend = F,
                    row_names_gp = gpar(fontsize = 8),
                    row_title_gp = gpar(fontsize = 12),
                    row_title_rot = 0,
                    row_title_side = "right", 
                    row_split = all.markers.top20$cluster,
                    cluster_row_slices = F,
                    top_annotation = ha,
                    use_raster = TRUE,
                    raster_resize_mat = TRUE,
                    raster_quality = 4)
        p <- draw(p)
        dev.off()
        
      }
    }
  }


#' batch corrected on data with combat on sample as covariate
FindMarkers_filter_combined_Stroma(data = data.sub.filt, 
                               slot = "scale.data", # data or scale.data
                               feature.list = feature.list, 
                               test = test,
                               out.name = "harmony_sample_log2SCTcounts1_Stroma")

FindAllMarkers_filter_combined_Stroma(data = data.sub.filt, 
                                  slot = "scale.data",
                                  feature.list = feature.list, 
                                  test = test,
                                  out.name = "harmony_sample_log2SCTcounts1_Stroma")
