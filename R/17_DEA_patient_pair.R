#' Perform different gene expression analysis using Seurat functions:
#' FindAllMarkers(): Gene expression markers for all identity classes;
#' FindMarker(): Gene expression markers of identity classes 
#' features 

library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/17_geneAnnotations.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/get_marker.R", sep = "/"))

file.name <- "OVisium_SCT_merged" 
cluster.ident <- "harmony_SCT_res_0.6"
setwd(rds.dir)
data <- readRDS(data, file = paste0(file.name, "_deg.rds"))


#' In order to test all gene, we set min.diff.pct = 0, logfc.threshold=-Inf
#' only compare the fimbrial and proximal from the same side (excluded the
#' proximal tissue from the other side) of patient no. 1, 3, 6, 10
DE <- list()
DE.TOP <- list()

#' choose different tissue origin to perform comparison
ident.1 = "Otherside" 
ident.2 = "Proximal"

#' patient no. 4 developed HGSC 6 years after the salpingectomy
for (p in c(6)) {
  Idents(data) <- "patient"
  pat.data <- SubsetSTData(data, idents = p)
  patient <- paste("Patient", p, sep = "_")
  
  #' QC of features vs UMI 
  gene_attr <- 
    data.frame(nUMI = Matrix::rowSums(pat.data@assays[["SCT"]]@counts), 
              nSpots = Matrix::rowSums(pat.data@assays[["SCT"]]@counts > 0))
  bks <- c(min(gene_attr$nUMI),1, 10, 50, 100, 200, 500, 1000, 10000,
           max(gene_attr$nUMI))
  setwd(deg.dir)
  png(file = paste(patient, "SCT_nUMI_per_feature.png", sep = "_"), width = 2500, 
      height = 1300, res = 300) 
  print(
  ggplot() +
    geom_histogram(data = gene_attr, aes(nUMI), 
                   fill = "red", alpha = 0.7, bins = 50) +
    scale_x_log10(breaks=bks,labels=bks) +
    xlab("SCT Normalized UMI Counts (log10 scale)") +
    ggtitle("Total counts per feature")
  )
  dev.off()
  
  #' choose features for DEA
  vfeatures <- VariableFeatures(pat.data)
  test <- c("MAST", "roc")
 
    for(t in seq_along(test)) {
      features <- vfeatures
      out.dir <- paste(deg.dir, file.name, cluster.ident, 
                       "Variable_features", test[t], 
                       paste(patient, ident.1, ident.2, sep = "_"), sep = "/") 
      dir.create(out.dir, recursive = T)
      setwd(out.dir)
      
  #' Find all markers between fimbrial vs proximal in PCA-Harmony
  #' log2FC threshold default 0.1
      Idents(pat.data) <- "origin"
      f.markers <- FindMarkers(pat.data, assay = "SCT", 
                               features = features,
                               ident.1 = ident.1,
                               ident.2 = ident.2,
                               min.diff.pct = 0,
                               logfc.threshold=-Inf,
                               recorrect_umi = FALSE,
                               test.use = test[t]) %>% 
        rownames_to_column() %>% 
        dplyr::rename("gene" = 1) %>%
        mutate(cluster = ident.1) 
      
      if (!is.null(f.markers$p_val)) {
        f.markers.filt <- f.markers %>% 
          mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>% 
          dplyr::filter(!grepl("^MT-", gene), 
                        !grepl("^RP[SL]", gene),
                        !grepl("^MTRNR", gene),
                        !grepl("^LINC", gene)) %>% 
          left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                    by = c("gene" = "gene_name"))
        
        f.markers.top <- f.markers.filt %>% 
          dplyr::filter(p_val_adj < 0.01, 
                        avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75)) 
      } else {
        print("This is ROC analysis result")
        f.markers.filt <- f.markers %>% 
          dplyr::filter(!grepl("^MT-", gene), 
                        !grepl("^RP[SL]", gene),
                        !grepl("^MTRNR", gene),
                        !grepl("^LINC", gene)) %>% 
          left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                    by = c("gene" = "gene_name"))
        
        f.markers.top <- f.markers.filt %>% 
          dplyr::filter(power < 0.45 | power > 0.55, 
                        avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75))
      }
      
        write.table(f.markers.filt, 
                    file = paste(patient, ident.1, "filt.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        write.table(f.markers.top, 
                    file = paste(patient, ident.1, "top.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
  
        #' FTE clusters fimbria 
        #' min.diff.pct=0.25 gives only 1 gene
        #' use min.diff.pct=0.1
        FTE.f.markers <- FindMarkers(pat.data, assay = "SCT",
                                     features = features,
                                     ident.1 = "0_FTE", 
                                     ident.2 = "1_FTE",
                                     group.by = "Visium_clusters",
                                     subset.ident = ident.1,
                                     recorrect_umi = FALSE,
                                     min.diff.pct = 0,
                                     logfc.threshold=-Inf,
                                     test.use = test[t]) %>% 
          rownames_to_column() %>% 
          dplyr::rename("gene" = 1) %>%
          mutate(cluster = "0_FTE") 
        
        if (!is.null(FTE.f.markers$p_val)) {
          FTE.f.markers.filt <- FTE.f.markers %>% 
            mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          FTE.f.markers.top <- FTE.f.markers.filt %>% 
            dplyr::filter(p_val_adj < 0.01, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75)) 
        } else {
          print("This is ROC analysis result")
          FTE.f.markers.filt <- FTE.f.markers %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          FTE.f.markers.top <- FTE.f.markers.filt %>% 
            dplyr::filter(power < 0.45 | power > 0.55, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75))
        }
        
        
        write.table(FTE.f.markers.filt, 
                    file = paste(patient, ident.1, "Cluster_0vs1_filt.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        write.table(FTE.f.markers.top, 
                    file = paste(patient, ident.1, "Cluster_0vs1_top.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        
        
        #' Immune clusters fimbria 
        Immune.f.markers <- FindMarkers(pat.data, assay = "SCT",
                                        features = features,
                                        ident.1 = "10_Stroma",
                                        ident.2 = "11_Stroma",
                                        group.by = "Visium_clusters",
                                        subset.ident = ident.1,
                                        recorrect_umi = FALSE,
                                        min.diff.pct = 0,
                                        logfc.threshold=-Inf,
                                        test.use = test[t]) %>% 
          rownames_to_column() %>% 
          dplyr::rename("gene" = 1) %>%
          mutate(cluster="10_Stroma") 
        
        if (!is.null(Immune.f.markers$p_val)) {
          Immune.f.markers.filt <- Immune.f.markers %>% 
            mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          Immune.f.markers.top <- Immune.f.markers.filt %>% 
            dplyr::filter(p_val_adj < 0.01, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75)) 
        } else {
          print("This is ROC analysis result")
          Immune.f.markers.filt <- Immune.f.markers %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          Immune.f.markers.top <- Immune.f.markers.filt %>% 
            dplyr::filter(power < 0.45 | power > 0.55, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75))
        }
        
        
        write.table(Immune.f.markers.filt, 
                    file = paste(patient, ident.1, "Cluster_10vs11_filt.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        write.table(Immune.f.markers.top, 
                    file = paste(patient, ident.1, "Cluster_10vs11_top.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        
        
    #' Find fimbrial vs proximal markers in each cluster 
        Idents(pat.data) <- "Visium_clusters"
        all.f.markers <- get_marker(dataset = pat.data, assay = "SCT",
                                    features = features,
                                    clusters = c("0_FTE", "1_FTE", "3_Mix", "2_Stroma", 
                                                 "5_Stroma", "6_Stroma", "7_Stroma",
                                                 "8_Stroma", "9_Stroma", "10_Stroma",
                                                 "11_Stroma"),
                                    comparison = "origin",
                                    ident.1 = ident.1,
                                    ident.2 = ident.2,
                                    recorrect_umi = FALSE, 
                                    min.diff.pct = 0,
                                    logfc.threshold=-Inf,
                                    test.use = test[t]) 
        
        all.f.markers$cluster<-
          factor(all.f.markers$cluster, 
                 levels = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                            "5_Stroma", "6_Stroma", "7_Stroma",
                            "8_Stroma", "9_Stroma", "10_Stroma",
                            "11_Stroma"))
        
        if (!is.null(all.f.markers$p_val)) {
          all.f.markers.filt <- all.f.markers %>% 
            mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          all.f.markers.top <- all.f.markers.filt %>% 
            dplyr::filter(p_val_adj < 0.01, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75)) 
        } else {
          print("This is ROC analysis result")
          all.f.markers.filt <- all.f.markers %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          all.f.markers.top <- all.f.markers.filt %>% 
            dplyr::filter(power < 0.45 | power > 0.55, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75))
        }
        
        
        write.table(all.f.markers.filt, 
                    file = paste(patient, ident.1, "Clusters_filt.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        write.table(all.f.markers.top, 
                    file = paste(patient, ident.1, "Clusters_top.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        
        #' Find differential expressed markers between 2 FTE clusters 
        #' min.diff.pct at 0.25 gives only 16 genes, set 0.1 instead
        FTE.markers <- FindMarkers(pat.data, assay = "SCT",
                                   features = features,
                                   ident.1 = "0_FTE", 
                                   ident.2 = "1_FTE",
                                   recorrect_umi = FALSE,
                                   min.diff.pct = 0,
                                   logfc.threshold=-Inf,
                                   test.use = test[t]) %>% 
          rownames_to_column() %>% 
          dplyr::rename("gene" = 1) %>%
          mutate(cluster="0_FTE") 
        
        if (!is.null(FTE.markers$p_val)) {
          FTE.markers.filt <- FTE.markers %>% 
            mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          FTE.markers.top <- FTE.markers.filt %>% 
            dplyr::filter(p_val_adj < 0.01, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75)) 
        } else {
          print("This is ROC analysis result")
          FTE.markers.filt <- FTE.markers %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          FTE.markers.top <- FTE.markers.filt %>% 
            dplyr::filter(power < 0.45 | power > 0.55, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75))
        }
        
        
        write.table(FTE.markers.filt, 
                    file = paste(patient, "Cluster_0vs1_filt.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        write.table(FTE.markers.top, 
                    file = paste(patient, "Cluster_0vs1_top.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        
        #' between 2 Immune clusters 
        Immune.markers <- FindMarkers(pat.data, assay = "SCT",
                                      features = features,
                                      ident.1 = "10_Stroma",
                                      ident.2 = "11_Stroma",
                                      recorrect_umi = FALSE,
                                      min.diff.pct = 0,
                                      logfc.threshold=-Inf,
                                      test.use = test[t]) %>% 
          rownames_to_column() %>% 
          dplyr::rename("gene" = 1) %>%
          mutate(cluster = "10_Stroma")
        
        if (!is.null(Immune.markers$p_val)) {
          Immune.markers.filt <- Immune.markers %>% 
            mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          Immune.markers.top <- Immune.markers.filt %>% 
            dplyr::filter(p_val_adj < 0.01, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75)) 
        } else {
          print("This is ROC analysis result")
          Immune.markers.filt <- Immune.markers %>% 
            dplyr::filter(!grepl("^MT-", gene), 
                          !grepl("^RP[SL]", gene),
                          !grepl("^MTRNR", gene),
                          !grepl("^LINC", gene)) %>% 
            left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
                      by = c("gene" = "gene_name"))
          
          Immune.markers.top <- Immune.markers.filt %>% 
            dplyr::filter(power < 0.45 | power > 0.55, 
                          avg_log2FC >= log2(1.5) | avg_log2FC <= log2(0.75))
        }
        
        
        write.table(Immune.markers.filt, 
                    file = paste(patient, "Cluster_10vs11_filt.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        write.table(Immune.markers.top, 
                    file = paste(patient, "Cluster_10vs11_top.csv", 
                                 sep = "_"),
                    sep =",", quote = F, row.names = F, col.names = T)
        
        
        #' plot de gene using volcano plot
         de <- list()
         de[[paste(ident.1, "vs", ident.2, sep = " ")]] <- f.markers.filt
         de[[paste(ident.1, "vs", ident.2, "in individual clusters", sep = " ")]] <- all.f.markers.filt
         de[[paste("Cluster 0_FTE vs 1_FTE in", ident.1, sep = " ")]] <- FTE.f.markers.filt
         de[[paste("Cluster10_Stroma vs 11_Stroma in", ident.1, sep = " ")]] <- Immune.f.markers.filt
         de[["Cluster 0_FTE vs 1_FTE"]] <- FTE.markers.filt
         de[["Cluster 10_Stroma vs 11_Stroma"]] <- Immune.markers.filt
         DE[[paste("Pat", p, test[t], ident.1, ident.2, sep = "_")]] <- de
            
        
        for (n in seq_along(names(de))) {
          #' enhancevolcano package 
          #' https://github.com/kevinblighe/EnhancedVolcano
          if (!is.null(de[[n]]$p_val)) {
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
                facet_wrap(. ~ cluster, ncol = 3)
              )
            dev.off()
          } else {
            #' result from ROC result
            #' myAUC : area under curve, at least 0.5
            #' classification power : (abs(AUC-0.5)) * 2), 
            #' 0.5 means no predictive power;
            #' 1 or 0 is positive and negative power
            if ( n != 2 ) {
              png(file = paste0(patient, "_volcano_", n, ".png"), 
                  width = 3500, height = 3500, res = 400)
            } else {
              png(file = paste0(patient, "_volcano_", n, ".png"), 
                  width = 8000, height = 8000, res = 500)
            }
            de[[n]]$power.log <- 1/(10^de[[n]]$power)
            keyvals <- ifelse(de[[n]]$power > 0.45 & de[[n]]$power < 0.55, "grey30",
              ifelse(de[[n]]$avg_log2FC < -log2(1.5), 'royalblue',
              ifelse(de[[n]]$avg_log2FC > log2(1.5), 'red2',"forestgreen")))
                     keyvals[is.na(keyvals)] <- 'forestgreen'
                       names(keyvals)[keyvals == 'red2'] <- 'Up'
                       names(keyvals)[keyvals == 'forestgreen'] <- 'No change'
                       names(keyvals)[keyvals == 'grey30'] <- 'No power'
                       names(keyvals)[keyvals == 'royalblue'] <- 'Down'
            print(
              EnhancedVolcano(de[[n]],
                              lab = de[[n]]$gene,
                              x = 'avg_log2FC',
                              y = 'power.log',
                              selectLab = de[[n]]$gene[which(names(keyvals) %in% c('Up', 'Down'))],
                              ylim = c(0, 1),
                              title = patient,
                              subtitle = names(de[n]),
                              pCutoff = c(1/(10^0.45), 1/(10^0.55)),
                              FCcutoff = log2(1.5),
                              pointSize = 3.0,
                              labSize = 6.0,
                              colCustom = keyvals) +
                labs(x = "avg_Log2(Fold Change)", y = "Seurat ROC power") +
                annotate("text", x= max(de[[n]]$avg_log2FC + 0.8), 
                         y = 0.5, label="No prediction power") +
                facet_wrap(. ~ cluster, ncol = 3)
            )
            dev.off()
          }
        }
        
        #' save top de for violin plot later across all patients
        de.top <- list()
        de.top[[paste(ident.1, "vs", ident.2, sep = " ")]] <- f.markers.top
        de.top[[paste(ident.1, "vs", ident.2, "in individual clusters", sep = " ")]] <- all.f.markers.top
        de.top[[paste("Cluster 0_FTE vs 1_FTE in", ident.1, sep = " ")]] <- FTE.f.markers.top
        de.top[[paste("Cluster10_Stroma vs 11_Stroma in", ident.1, sep = " ")]] <- Immune.f.markers.top
        de.top[["Cluster 0_FTE vs 1_FTE"]] <- FTE.markers.top
        de.top[["Cluster 10_Stroma vs 11_Stroma"]] <- Immune.markers.top
        DE.TOP[[paste("Pat", p, test[t], ident.1, ident.2, sep = "_")]] <- de.top
    }
}

out.dir <- paste(deg.dir, file.name, cluster.ident, 
                 "Variable_features", sep = "/") 
saveRDS(DE, file="Patients.DEA.filt.list.rds")
saveRDS(DE.TOP, file="Patients.DEA.top.list.rds")


#' violin plot for interesting genes
png(file = paste0(patient, "_violin_origins_", n, ".png"), 
    width = 5000, height = 5000, res = 400)
VlnPlot(pat.data, 
        group.by = 'Tissue_origin',
        features = features, 
        ncol=3, 
        pt.size = 0)
dev.off()

png(file = paste0(patient, "_violin_clusters_origins_", n, ".png"), 
    width = 5000, height = 5000, res = 400)
print(VlnPlot(pat.data, 
              idents = c("0_FTE", "1_FTE", "10_Stroma", "11_Stroma"),
              split.by = 'Tissue_origin',
              features = features, 
              ncol=3, 
              pt.size = 0)
)
dev.off()

   
  
