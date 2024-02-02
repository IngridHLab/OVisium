
#' FindConservedMarkers(): Finds markers that are conserved between the groups;
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/17_geneAnnotations.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/get_marker.R", sep = "/"))

file.name <- "OVisium_SCT_merged" 
cluster.ident <- "harmony_SCT_res_0.6"

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
feature.list <- list(afeatures, vfeatures, vhfeatures, sfeatures, shfeatures)
names(feature.list) <- 
  c("Variable_features","Variable_a500UMI_features",
    "Scaled_features","Scaled_a500UMI_features","All_features")
test <- c("wilcox", "MAST", "DESeq2")

#' Find conserved markers between sample/library in each cluster
#' Note: there is no conserved marker in cluster 10_Stroma (skipped) 
all.cons.markers <- 
  get_conserved_marker(dataset = data, assay = "SCT",
                       features = features,
                       clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                    "5_Stroma", "6_Stroma", "7_Stroma",
                                    "8_Stroma", "9_Stroma", "11_Stroma"),
                       group = "sample",
                       min.pct = 0.25,
                       logfc.threshold = log2(1.5),
                       test.use = test[t]) %>% 
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) %>%
  mutate(rank=-log10(minimump_p_val)*(avg_log2FC))

all.cons.markers.filt <- all.cons.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                minimump_p_val < 0.01)

all.cons.markers.filt.top100 <- all.cons.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = abs(avg_log2FC), n = 100)

write.table(all.cons.markers.filt, 
            file = "Clusters.sample.conserved.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.cons.markers.filt.top100, 
            file = "Clusters.sample.conserved.markers.filt.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)

#' Find conserved markers between fimbrial and proximal in each cluster 
all.f.cons.markers <- 
  get_conserved_marker(dataset = data, assay = "SCT",
                       features = features,
                       clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                    "5_Stroma", "6_Stroma", "7_Stroma",
                                    "8_Stroma", "9_Stroma", "10_Stroma",
                                    "11_Stroma"),
                       group = "Tissue_origin",
                       min.pct = 0.25,
                       logfc.threshold = log2(1.5),
                       test.use = test[t]) %>% 
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) %>%
  mutate(rank=-log10(minimump_p_val)*(avg_log2FC))

all.f.cons.markers.filt <- all.f.cons.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                minimump_p_val < 0.01)

all.f.cons.markers.filt.top100 <- all.f.cons.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = abs(avg_log2FC), n = 100)

write.table(all.f.cons.markers.filt, 
            file = "Clusters.fimbrial.conserved.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.f.cons.markers.filt.top100, 
            file = "Clusters.fimbrial.conserved.markers.filt.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)

#' Find conserved markers between BRCA1/2 in each cluster 
all.m.cons.markers <- 
  get_conserved_marker(dataset = data, assay = "SCT",
                       features = features,
                       clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                    "5_Stroma", "6_Stroma", "7_Stroma",
                                    "8_Stroma", "9_Stroma", "10_Stroma",
                                    "11_Stroma"),
                       group = "mutation",
                       min.pct = 0.25,
                       logfc.threshold = log2(1.5),
                       test.use = test[t]) %>% 
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) %>%
  mutate(rank=-log10(minimump_p_val)*(avg_log2FC))

all.m.cons.markers.filt <- all.m.cons.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                minimump_p_val < 0.01)

all.m.cons.markers.filt.top100 <- all.m.cons.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = abs(avg_log2FC), n = 100) 

write.table(all.m.cons.markers.filt, 
            file = "Clusters.BRCA12.conserved.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.m.cons.markers.filt.top100, 
            file = "Clusters.BRCA12.conserved.markers.filt.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)

#' Find conserved markers between different variants in each cluster 
all.v.cons.markers <- 
  get_conserved_marker(dataset = data, assay = "SCT",
                       features = features,
                       clusters = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                    "5_Stroma", "6_Stroma", "7_Stroma",
                                    "8_Stroma", "9_Stroma", "10_Stroma",
                                    "11_Stroma"),
                       group = "variant",
                       min.pct = 0.25,
                       logfc.threshold = log2(1.5),
                       test.use = test[t]) %>% 
  left_join(. , y = unique(annotations[, c("gene_name", "description")]), 
            by = c("gene" = "gene_name")) %>%
  mutate(rank=-log10(minimump_p_val)*(avg_log2FC))

all.v.cons.markers.filt <- all.v.cons.markers %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                minimump_p_val < 0.01)

all.v.cons.markers.filt.top100 <- all.v.cons.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = abs(avg_log2FC), n = 100)

write.table(all.v.cons.markers.filt, 
            file = "Clusters.variants.conserved.markers.filt.csv",
            sep =",", quote = F, row.names = F, col.names = T)

write.table(all.v.cons.markers.filt.top100, 
            file = "Clusters.variants.conserved.markers.filt.top100.csv",
            sep =",", quote = F, row.names = F, col.names = T)