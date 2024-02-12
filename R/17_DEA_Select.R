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
Idents(data) <- "Tissue_origin"
f.data <- SubsetSTData(data, idents = "Fimbrial")


#' select features
afeatures <- row.names(data@assays[["SCT"]]@counts)
vfeatures <- VariableFeatures(data)
sfeatures <- row.names(data@assays[["SCT"]]@scale.data)
gene_attr <- 
  data.frame(nUMI = Matrix::rowSums(data@assays[["SCT"]]@counts), 
             nSpots = Matrix::rowSums(data@assays[["SCT"]]@counts > 0))
bks <- c(min(gene_attr$nUMI),1, 10, 50, 100, 200, 500, 1000, 10000, 
         max(gene_attr$nUMI))
setwd(deg.dir)
png(file = "SCT_nUMI_per_feature.png", width = 2500, height = 1300, res = 300) 
print(
ggplot() +
  geom_histogram(data = gene_attr, aes(nUMI), 
                 fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10(breaks=bks,labels=bks) +
  xlab("SCT Normalized UMI Counts (log10 scale)") +
  ggtitle("Total counts per feature")
)
dev.off()
hfeatures <- row.names(gene_attr[which(gene_attr$nUMI>500),]) 
shfeatures <- unique(sfeatures[sfeatures%in%hfeatures])
vhfeatures <- unique(vfeatures[vfeatures%in%hfeatures])
feature.list <- list(vfeatures, vhfeatures, sfeatures, shfeatures, afeatures)
names(feature.list) <- 
  c("Variable_features","Variable_a500UMI_features",
    "Scaled_features","Scaled_a500UMI_features","All_features")
test <- c("wilcox", "MAST")



#' Find all markers of all graph-based clusters in PCA-Harmony
for(i in seq_along(feature.list)) {
  for(t in seq_along(test)) {
features <- feature.list[[i]]
out.dir <- paste(deg.dir, file.name, cluster.ident, 
                 names(feature.list[i]), test[t], sep = "/") 
dir.create(out.dir, recursive = T)
setwd(out.dir)
Idents(data) <- "Visium_clusters"
all.markers <- FindAllMarkers(data, assay = "SCT", 
                              features = features,
                              min.diff.pct = 0.25,
                              logfc.threshold = log2(1.5),
                              test.use = test[t]) %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 
  
all.markers.filt <- all.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

all.markers.filt.top100 <- all.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = rank, n = 100) 

write.table(all.markers.filt, 
            file = "Clusters.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.markers.filt.top100, 
            file = "Clusters.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)



#' Find differential expressed markers between 2 clusters
#' FTE clusters 
#' min.diff.pct at 0.25 has no gene found, set 10 instead
FTE.markers <- FindMarkers(data, assay = "SCT",
                           features = features,
                           ident.1 = "0_FTE", 
                           ident.2 = "1_FTE",
                           min.diff.pct = 0.10,
                           logfc.threshold = log2(1.5),
                           test.use = test[t]) %>% 
  rownames_to_column() %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC), cluster="0_FTE") %>%
  dplyr::rename("gene" = 1) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 

FTE.markers.filt <- FTE.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene), 
                p_val_adj  < 0.01)

FTE.markers.filt.top100 <- FTE.markers.filt %>% 
  dplyr::slice_max(order_by = rank, n = 100) 

write.table(FTE.markers.filt, 
            file = "Clusters.0vs1.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)


write.table(FTE.markers.filt.top100, 
            file = "Clusters.0vs1.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)


#' FTE clusters fimbria 
Idents(f.data) <- "Visium_clusters"
FTE.f.markers <- FindMarkers(f.data, assay = "SCT",
                             features = features,
                             ident.1 = "0_FTE", 
                             ident.2 = "1_FTE",
                             group.by = "Visium_clusters",
                             recorrect_umi = FALSE,
                             min.diff.pct = 0.10,
                             logfc.threshold = log2(1.5),
                             test.use = test[t]) %>% 
  rownames_to_column() %>% 
  mutate(rank = -log10(p_val)*sign(avg_log2FC), cluster = "0_FTE") %>%
  dplyr::rename("gene" = 1) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 

FTE.f.markers.filt <- FTE.f.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

FTE.f.markers.filt.top100 <- FTE.f.markers.filt %>% 
  dplyr::slice_max(order_by = rank, n = 100) 

write.table(FTE.f.markers.filt, 
            file = "Clusters.0vs1.fimbrial.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(FTE.f.markers.filt.top100, 
            file = "Clusters.0vs1.fimbrial.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)

#' Immune clusters 
immune.markers <- FindMarkers(data, assay = "SCT",
                              features = features,
                              ident.1 = "11_Stroma",
                              ident.2 = "10_Stroma",
                              min.diff.pct = 0.25,
                              logfc.threshold = log2(1.5),
                              test.use = test[t]) %>% 
  rownames_to_column() %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC), cluster = "10_Stroma") %>%
  dplyr::rename("gene" = 1) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name"))

immune.markers.filt <- immune.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene), 
                p_val_adj  < 0.01)

immune.markers.filt.top100 <- immune.markers.filt %>% 
  dplyr::slice_max(order_by = rank, n = 100) 

write.table(immune.markers.filt, 
            file = "Cluster.10vs11.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(immune.markers.filt.top100, 
            file = "Cluster.10vs11.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)

#' Immune clusters fimbria 
Immune.f.markers <- FindMarkers(f.data, assay = "SCT",
                                features = features,
                                ident.1 = "10_Stroma",
                                ident.2 = "11_Stroma",
                                group.by = "Visium_clusters",
                                recorrect_umi = FALSE,
                                min.diff.pct = 0.25,
                                logfc.threshold = log2(1.5),
                                test.use = test[t]) %>% 
  rownames_to_column() %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC), cluster="10_Stroma") %>%
  dplyr::rename("gene" = 1) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 

Immune.f.markers.filt <- Immune.f.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

Immune.f.markers.filt.top100 <- Immune.f.markers.filt %>% 
  dplyr::slice_max(order_by = rank, n = 100) 

write.table(Immune.f.markers.filt, 
            file = "Cluster.10vs11.fimbrial.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(Immune.f.markers.filt.top100, 
            file = "Cluster.10vs11.fimbrial.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)


#' Find fimbrial vs proximal markers in each cluster 
all.f.markers <- get_marker(dataset = data, assay = "SCT",
                            features = features,
                            clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                         "5_Stroma", "6_Stroma", "7_Stroma",
                                         "8_Stroma", "9_Stroma", "10_Stroma",
                                         "11_Stroma"),
                            comparison = "Tissue_origin",
                            ident.1 = "Fimbrial",
                            recorrect_umi = FALSE, 
                            min.diff.pct = 0.25,
                            logfc.threshold = log2(1.5),
                            test.use = test[t]) %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name"))

all.f.markers.filt <- all.f.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

all.f.markers.filt.top100 <- all.f.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = rank, n = 100)

write.table(all.f.markers.filt, 
            file = "Clusters.FimbrialvsProximal.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.f.markers.filt.top100, 
            file = "Clusters.FimbrialvsProximal.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)



#' Find BRCA1 markers in each cluster 
all.b1.markers <- get_marker(dataset = data, assay = "SCT",
                             features = features,
                            clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                         "5_Stroma", "6_Stroma", "7_Stroma",
                                         "8_Stroma", "9_Stroma", "10_Stroma",
                                         "11_Stroma"),
                            comparison = "mutation",
                            ident.1 = "BRCA1",
                            recorrect_umi = FALSE, 
                            min.diff.pct = 0.25,
                            logfc.threshold = log2(1.5),
                            test.use = test[t]) %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 

all.b1.markers.filt <- all.b1.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

all.b1.markers.filt.top100 <- all.b1.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = rank, n = 100)

write.table(all.b1.markers.filt, 
            file = "Clusters.BRCA1vs2.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.b1.markers.filt.top100, 
            file = "Clusters.BRCA1vs2.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T) 

#' Find BRCA1 fimbrial markers in each cluster 
all.b1.f.markers <- get_marker(dataset = f.data, assay = "SCT",
                             features = features,
                             clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                          "5_Stroma", "6_Stroma", "7_Stroma",
                                          "8_Stroma", "9_Stroma", "10_Stroma",
                                          "11_Stroma"),
                             comparison = "mutation",
                             ident.1 = "BRCA1",
                             recorrect_umi = FALSE, 
                             min.diff.pct = 0.25,
                             logfc.threshold = log2(1.5),
                             test.use = test[t]) %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 

all.b1.f.markers.filt <- all.b1.f.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

all.b1.f.markers.filt.top100 <- all.b1.f.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = rank, n = 100)

write.table(all.b1.f.markers.filt, 
            file = "Clusters.BRCA1vs2.fimbrial.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.b1.f.markers.filt.top100, 
            file = "Clusters.BRCA1vs2.fimbrial.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T) 

#' Find variant c.3048_3052dup (patients PO2 and PO40) markers in each cluster 
all.v1.markers <- get_marker(dataset = data, assay = "SCT",
                             features = features,
                             clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                          "5_Stroma", "6_Stroma", "7_Stroma",
                                          "8_Stroma", "9_Stroma", "10_Stroma",
                                          "11_Stroma"),
                             comparison = "variant",
                             ident.1 = "c.3048_3052dup",
                             recorrect_umi = FALSE, 
                             min.diff.pct = 0.25,
                             logfc.threshold = log2(1.5),
                             test.use = test[t]) %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 

all.v1.markers.filt <- all.v1.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

all.v1.markers.filt.top100 <- all.v1.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = rank, n = 100) 

write.table(all.v1.markers.filt, 
            file = "Clusters.c3048_3052dupvsRest.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.v1.markers.filt.top100, 
            file = "Clusters.c3048_3052dupvsRest.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T) 

#' Find variant c.3048_3052dup (patients PO2 and PO40) 
#' fimbrial markers in each cluster 
all.v1.f.markers <- get_marker(dataset = f.data, assay = "SCT",
                             features = features,
                             clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                          "5_Stroma", "6_Stroma", "7_Stroma",
                                          "8_Stroma", "9_Stroma", "10_Stroma",
                                          "11_Stroma"),
                             comparison = "variant",
                             ident.1 = "c.3048_3052dup",
                             recorrect_umi = FALSE, 
                             min.diff.pct = 0.25,
                             logfc.threshold = log2(1.5),
                             test.use = test[t]) %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 

all.v1.f.markers.filt <- all.v1.f.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

all.v1.f.markers.filt.top100 <- all.v1.f.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = rank, n = 100) 

write.table(all.v1.f.markers.filt, 
            file = "Clusters.c3048_3052dupvsRest.fimbrial.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.v1.f.markers.filt.top100, 
            file = "Clusters.c3048_3052dupvsRest.fimbrial.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T) 

#' Find variant c.1082_1092del (patients PO6, PO13 and PO20) 
#' markers in each cluster 
all.v2.markers <- get_marker(dataset = data, assay = "SCT",
                             features = features,
                             clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                          "5_Stroma", "6_Stroma", "7_Stroma",
                                          "8_Stroma", "9_Stroma", "10_Stroma",
                                          "11_Stroma"),
                             comparison = "variant",
                             ident.1 = "c.1082_1092del",
                             recorrect_umi = FALSE, 
                             min.diff.pct = 0.25,
                             logfc.threshold = log2(1.5),
                             test.use = test[t]) %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 

all.v2.markers.filt <- all.v2.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

all.v2.markers.filt.top100 <- all.v2.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = rank, n = 100)

write.table(all.v2.markers.filt, 
            file = "Clusters.c1082_1092delvsRest.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.v2.markers.filt.top100, 
            file = "Clusters.c1082_1092delvsRest.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)

#' Find variant c.1082_1092del (patients PO6, PO13 and PO20) 
#' fimbrial markers in each cluster 
all.v2.f.markers <- get_marker(dataset = f.data, assay = "SCT",
                               features = features,
                               clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                            "5_Stroma", "6_Stroma", "7_Stroma",
                                            "8_Stroma", "9_Stroma", "10_Stroma",
                                            "11_Stroma"),
                               comparison = "variant",
                               ident.1 = "c.1082_1092del",
                               recorrect_umi = FALSE, 
                               min.diff.pct = 0.25,
                               logfc.threshold = log2(1.5),
                               test.use = test[t]) %>% 
  mutate(rank=-log10(p_val)*sign(avg_log2FC)) %>%
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) 

all.v2.f.markers.filt <- all.v2.f.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                p_val_adj < 0.01)

all.v2.f.markers.filt.top100 <- all.v2.f.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = rank, n = 100)

write.table(all.v2.f.markers.filt, 
            file = "Clusters.c1082_1092delvsRest.fimbrial.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.v2.f.markers.filt.top100, 
            file = "Clusters.c1082_1092delvsRest.fimbrial.markers.filt.rank.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)
}
}

#' Generate volcano plot for DEA result
png(file = paste0(patient, "_volcano_", n, ".png"), 
    width = 8000, height = 8000, res = 500)
EnhancedVolcano(Clusters_markers_filt,
                lab = Clusters_markers_filt$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "patient",
                subtitle = names(de[n]),
                pCutoff = 0.01,
                FCcutoff = 1,                              
                pointSize = 3.0,
                labSize = 6.0,
                legendLabels=c('Not sig.',
                               'avg_Log2(FC)',
                               'p_adj_val',
                               'p_adj_val & avg_Log2(FC)')) +
  labs(x = "avg_Log2(Fold Change)", y = "-log10 (Adjust P_value)") +
  facet_wrap(. ~ cluster, ncol = 3)
dev.off()