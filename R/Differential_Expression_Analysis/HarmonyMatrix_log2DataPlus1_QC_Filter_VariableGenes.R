#' Perform filtering on the variable features.
#' QC on SCT raw counts
#' filtering based SCT QC
library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

file.name <- "OVisium_SCT_merged" 
cluster.ident <- "harmony_SCT_res_0.6"
setwd(rds.dir)
#' use the data which generate PCA and for Harmony-clustering
data <- readRDS(file = paste0(file.name, "_clust.rds"))

#' select variable features and try regressing out batch via Harmony
#' based on "sample"
vfeatures <- VariableFeatures(data, assay = "SCT")
m = log(as.matrix(data@assays[["SCT"]]@counts[vfeatures, ])+1, 2)
meta_data <- as.data.frame(as.matrix(data@meta.data)) 
m.harmony <- HarmonyMatrix(
  data_mat  = m,       # Matrix with coordinates for each cell (row) along many PCs (columns)
  meta_data = meta_data, # Dataframe with information for each cell (row)
  vars_use  = "sample", # Column in meta_data that defines dataset for each cell
  do_pca    = FALSE      # Since we are providing PCs, do not run PCA
)

#' m.harmony <- read_csv("~/OVisium/SRIQ/All_clusters_vfeatureas_SCTcounts_log2+1_harmony_samples.csv") %>% as.matrix()
data@assays[["SCT"]]@data <- Matrix(as.matrix(m.harmony), sparse = T)
data@assays[["SCT"]]@counts <- Matrix(pmax(round(2^as.matrix(data@assays[["SCT"]]@data)-1, digits = 0),0), sparse = T)

#' subset only 11 clusters
source(paste(home, "OVisium/R/Help_functions/Subset_Rename_Clusters.R", sep = "/"))

#' filter the variable features based on distribution
gene_attr <- 
  data.frame(nUMI = Matrix::rowSums(data.sub@assays[["SCT"]]@counts), 
             nSpots = Matrix::rowSums(data.sub@assays[["SCT"]]@counts > 0))

#' check low expressed genes and bleed-over genes
#' https://github.com/satijalab/seurat/issues/4958
#' https://github.com/satijalab/seurat/issues/2496

counts <- data.frame(GetAssayData(data.sub, slot="counts", assay="SCT"))
log2counts1 <- data.frame(GetAssayData(data.sub, slot="data", assay="SCT"))

counts$max <- apply(counts, 1, max, na.rm = T) 
counts$pct <- rowMeans(counts>0 )*100
counts <- counts %>% relocate(c(max, pct))

log2counts1$pct98 <- apply(log2counts1, 1, quantile, probs = 0.98, na.rm = T)
log2counts1 <- log2counts1 %>% relocate(pct98)

#' plot distribution of our genes
bks <- c(min(gene_attr$nUMI), 300, 500, 1000, 10000, 100000, 1000000,
         max(gene_attr$nUMI))
p1 <-
  ggplot(data = gene_attr, aes(nUMI)) +
  geom_histogram(aes(y = ..density..), fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10(breaks=bks,labels=bks) +
  geom_vline(xintercept = 500, color = "yellow", lwd = 2) +
  geom_density(color = "blue", lwd = 1.2) +
  xlab("Total Harmony SCT counts (log10 scale)") +
  ggtitle("5894 genes with total counts >500")

bks <- c(min(counts$max),3,4,10,50,100,200,500,1000,max(counts$max))
p2 <-
  ggplot(data = counts, aes(max)) +
  geom_histogram(aes(y = ..density..), fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10(breaks=bks,labels=bks) +
  geom_vline(xintercept = 4, color = "yellow", lwd = 2) +
  geom_density(color = "blue", lwd = 1.2) +
  xlab("Max Harmony SCT counts (log10 scale)") +
  ggtitle("3843 genes with maximum counts >4")

p3 <-
  ggplot(data = counts, aes(pct)) +
  geom_histogram(aes(y = ..density..), fill = "red", alpha = 0.7, bins = 50) +
  geom_vline(xintercept = 1, color = "yellow", lwd = 2) +
  geom_density(color = "blue", lwd = 1.2) +
  xlab("Percentage of Cells %") +
  ggtitle("6117 genes with fraction of spots >1%")

p4 <-
  ggplot(data = log2counts1, aes(pct98)) +
  geom_histogram(aes(y = ..density..), fill = "red", alpha = 0.7, bins = 50) +
  geom_vline(xintercept = log2(2.5), color = "yellow", lwd = 2) +
  geom_vline(xintercept = log2(64), color = "green", lwd = 2) +
  geom_density(color = "blue", lwd = 1.2) +
  xlab("98 Percentile of Harmony log2(SCTcounts+1)") +
  ggtitle("2760 genes with 98 percentile log2(SCTcounts+1)>log2(2.5)")

setwd(paste(deg.dir, file.name, cluster.ident, sep = "/"))
png(file = "QC_per_vfeature_SCT_log2counts1_harmony.png", width = 5000, height = 2600, res = 300) 
ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, labels = c("A","B","C","D")) 
dev.off()

#' subset the vfeatures
vfeatures.500 <- row.names(gene_attr[gene_attr$nUMI>500,]) 
vfeatures.max <- row.names(counts[counts$max>4,])
vfeatures.pct <- row.names(counts[counts$pct>1,])
vfeatures.pct98 <- row.names(log2counts1[log2counts1$pct98>log2(2.5),])
vfeatures.pct98.bleed <- row.names(log2counts1[log2counts1$pct98>=6,]) 
#' skip filter out bleedover genes
[1] "RPS27"   "EEF1A1"  "RPS12"   "TAGLN"   "RPL41"   "RPLP1"   "RPL13"   "MYL9"    "RPL10"   "MT-CO1"  "MT-CO2"  "MT-ATP6" "MT-CO3" 
[14] "MT-ND3"  "MT-ND4" 

#' Last remove all mit
gene_attr_core <- gene_attr %>% 
  dplyr::filter(!grepl("^MT-", row.names(.)), 
                !grepl("^RP[SL]", row.names(.)),
                !grepl("^MTRNR", row.names(.)),
                !grepl("^LINC", row.names(.)))
vfeatures.core <- row.names(gene_attr_core)
vfeatures.filt <- Reduce(intersect,list(vfeatures.max, 
                                        vfeatures.pct, 
                                        vfeatures.pct98, 
                                        vfeatures.core))

#' plot those after filter away low expressed genes and high bleed-over expressed genes 
gene_attr_filt <- 
  dplyr::filter(gene_attr, row.names(gene_attr) %in% vfeatures.filt)

bks <- c(min(gene_attr_filt$nUMI), 5000, 10000, 100000, max(gene_attr_filt$nUMI))
png(file = "nUMI_per_vfeature_filt_SCT_log2counts1_harmony.png", width = 2500, height = 1300, res = 300) 
print(
  ggplot(data = gene_attr_filt, aes(nUMI)) +
    geom_histogram(aes(y = ..density..), fill = "red", alpha = 0.7, bins = 50) +
    scale_x_log10(breaks=bks,labels=bks) +
    geom_density(color = "blue", lwd = 1.2) +
    scale_x_log10(breaks=bks,labels=bks) +
    xlab("Total Harmony SCT counts (log10 scale)") +
    ggtitle("2644 variable Genes after filtering")
)
dev.off()

#' subset only the filtered variable features for combat
setwd(rds.dir) 
data.sub.filt <- SubsetSTData(data.sub, features = vfeatures.filt)
save(data.sub, data.sub.filt, vfeatures.filt, 
     file = "Variable_features_filt_SCT_log2counts+1_harmony.RData")

#' export the filtered harmony matrix for SRIQ
setwd(anno.dir)
m.harmony.filt = as.matrix(data.sub.filt@assays[["SCT"]]@data)
write.csv(m.harmony.filt, "11clusters_filt_vfeatureas_SCTcounts_log2+1_harmony_samples.csv", quote = F)
Idents(data.sub.filt) <- "sample"
cellInfo.1 <- data.frame(cluster=Idents(data.sub.filt))
Idents(data.sub.filt) <- cluster.ident
cellInfo.2 <- data.frame(sample=Idents(data.sub.filt))
cellInfo <- cbind(cellInfo.1, cellInfo.2) %>% rownames_to_column("barcode")
write.csv(cellInfo, "Filtered_spots_sample_cluster_annotation.csv", quote = F, row.names = F)


#' in this step the following genes have log2.pct98 > = 6 
[1] "RPS27"   "EEF1A1"  "RPS12"   "TAGLN"   "RPL41"   "RPLP1"   "RPL13"   "MYL9"    "RPL10"   "MT-CO1"  "MT-CO2"  "MT-ATP6" "MT-CO3" 
[14] "MT-ND3"  "MT-ND4"
