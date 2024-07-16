#' Box plot/Violin plot to compare QC matrix between samples
#' https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html
#' https://ludvigla.github.io/STUtility_web_site/Quality_Control.html
#' https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Plot_filtered_QC

library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))

#' Collect all genes coded on the mitochondrial genome and ribosomal proteins
#' plot the mito and ribo genes expression level spatially
merged <- readRDS(paste(rds.dir, "OVisium_merged.rds", sep = "/"))
Idents(merged) <- "sample"
num <- length(levels(merged))
merged[["Mito.percent"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged[["Ribo.percent"]] <- PercentageFeatureSet(merged, pattern = "^RP[SL]")
#' Percentage hemoglobin genes - includes all genes starting with HB except HBP.
merged[["Hb.percent"]] <- PercentageFeatureSet(merged, pattern = "^HB[^(P)]")


#' spatial plots
P <- list()
measures <- c("nFeature_RNA", "nCount_RNA", "Mito.percent", "Ribo.percent", 
              "Hb.percent")
for (i in seq_along(1:length(measures))) {
  
  P[[i]] <- ST.FeaturePlot(merged, sampleids = 1:18, features = measures[i], 
                           ncol = 5, show.sb = F, pt.size = 0.8) &
    theme(plot.title = element_text(hjust = 0.5, vjust = 2, 
                                    size=30)) &
    theme(text=element_text(size=30)) &
    theme(legend.title = element_blank(), legend.key.size = unit(2, 'cm'), 
          legend.key.width = unit(1.5, 'cm'))
}

setwd(qc.dir)
png(file = "Spatial_nFeature_RNA.png", 
    width = 8500, height = 7000, res = 500)
print(
  P[[1]] + ggtitle("nFeature_RNA")
)
dev.off()
png(file = "Spatial_nCount_RNA.png", 
    width = 8500, height = 7000, res = 500)
print(
  P[[2]] + ggtitle("nCount_RNA")
)
dev.off()
png(file = "Spatial_mito_percent.png", 
    width = 8500, height = 7000, res = 500)
print(
  P[[3]] + ggtitle("Mitochrondrial Expression %")
)
dev.off()
png(file = "Spatial_ribo_percent.png", 
    width = 8500, height = 7000, res = 500)
print(
  P[[4]] + ggtitle("Ribosomal Expression %")
)
dev.off()
png(file = "Spatial_hb_percent.png", 
    width = 8500, height = 7000, res = 500)
print(
  P[[5]] + ggtitle("Homoglobin Expression %")
)
dev.off()


#' violin plots
group_by <- "sample"
df <- FetchData(merged, vars = c(group_by, measures))
df$sample <- as.factor(df$sample)
p <- list()
for (i in seq_along(1:length(measures))) {
  
  p[[i]] <- ggplot(df, aes(x = .data[[group_by]],
                           y = .data[[measures[i]]],
                           fill = .data[[group_by]])) +
    scale_fill_manual(values = col_vector) +
    geom_violin() +
    geom_boxplot(alpha = 0.2) +
    labs(x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 2, 
                                    size=30),
          legend.position = "none",
          text=element_text(size=30))
}

setwd(qc.dir)
png(file = paste(qc.dir, "Vln_nFeature_RNA.png", sep = "/"),  
    width = 8500, height = 4500, res = 500) 
print(
  p[[1]] + ggtitle("nFeatures_RNA")
)
dev.off()
png(file = paste(qc.dir, "Vln_nCount_RNA.png", sep = "/"),  
    width = 8500, height = 4500, res = 500) 
print(
  p[[2]] + ggtitle("nCount_RNA")
)
dev.off()
png(file = paste(qc.dir, "Vln_mito_percent.png", sep = "/"),  
    width = 8500, height = 4500, res = 500) 
print(
  p[[3]] + ggtitle("Mitochrondrial expression %") 
)
dev.off()
png(file = paste(qc.dir, "Vln_ribo_percent.png", sep = "/"),  
    width = 8500, height = 4500, res = 500) 
print(
  p[[4]] + ggtitle("Ribosomal expression %") 
)
dev.off()
png(file = paste(qc.dir, "Vln_hb_percent.png", sep = "/"),  
    width = 8500, height = 4500, res = 500) 
print(
  p[[5]] + ggtitle("Homoglobin expression %") 
)
dev.off()

#' Bar plot to show distribution of data
gene_attr <- data.frame(nUMI = Matrix::rowSums(merged@assays$RNA@counts), 
                        nSpots = Matrix::rowSums(merged@assays$RNA@counts > 0))
p1 <- ggplot() +
  geom_histogram(data = merged[[]], aes(nFeature_RNA), 
                 fill = "red", alpha = 0.7, bins = 50) +
  xlab("nFeature_RNA") +
  ggtitle("Unique features per spot")
p2 <- ggplot() +
  geom_histogram(data = merged[[]], aes(nCount_RNA), 
                 fill = "red", alpha = 0.7, bins = 50) +
  xlab("nCount_RNA") +
  ggtitle("Total counts per spot")
p3 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nUMI), 
                 fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10() +
  xlab("nUMI (log10 scale)") +
  ggtitle("Total counts per feature")
p4 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nSpots),
                 fill = "red", alpha = 0.7,  bins = 50) +
  xlab("nSpot (nCount_RNA > 0)") +
  ggtitle("Total spots per feature")
setwd(qc.dir)
png(file = "QC_RNA_histo.png", width = 5000, height = 2600, res = 300) 
ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, labels = c("A","B","C","D"))
dev.off()

#' Data distribution density plots
#' Add number of genes per UMI for each spot to metadata
merged$log10FeaturesPerUMI <- log10(merged$nFeature_RNA) / log10(merged$nCount_RNA)
metadata <- merged@meta.data

#' Visualize the number of spot counts per sample
p1 <-metadata %>% 
    ggplot(aes(x=origin, fill=origin)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) +
    ylab("nSpots") + 
    xlab("") + 
  theme(plot.margin = margin(1,1,1,1, "cm"))

#' Visualize the number UMIs/transcripts per spot
p2 <- metadata %>% 
    ggplot(aes(color=origin, x=nCount_RNA, fill= origin)) + 
    geom_density(alpha = 0.2) +
    scale_x_log10() + 
    theme_classic() +
    xlab("nCount_RNA") +
    ylab("Spot density") +
    geom_vline(xintercept = 500) + 
  theme(plot.margin = margin(1,1,2,1, "cm"))

#' Visualize the distribution of genes detected per spot via histogram
p3 <- metadata %>% 
    ggplot(aes(color=origin, x=nFeature_RNA, fill= origin)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    xlab("nFeature_RNA") +
    ylab("Spot density") +
    geom_vline(xintercept = 500) + 
  theme(plot.margin = margin(1,1,2,1, "cm"))

#' Visualize the overall complexity of the gene expression by visualizing the 
#' genes detected per UMI (novelty score)
p4 <- metadata %>%
    ggplot(aes(x=log10FeaturesPerUMI, color = origin, fill=origin)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    xlab("log10FeaturesPerUMI") +
    ylab("Spot density") +
    geom_vline(xintercept = 0.83) + 
  theme(plot.margin = margin(1,1,2,1, "cm"))

#' Visualize the distribution of mitochondrial gene expression detected per spot
p5 <- metadata %>% 
    ggplot(aes(color=origin, x=Mito.percent, fill=origin)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    xlab("Mitochondrial %") +
    ylab("Spot density") +
    geom_vline(xintercept = 15) + 
  theme(plot.margin = margin(1,1,2,1, "cm"))

#' Visualize the correlation between genes detected and number of UMIs and 
#' determine whether strong presence of spots with low numbers of genes/UMIs
p6 <- metadata[which(metadata$log10FeaturesPerUMI>0.83),] %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=Mito.percent)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "red") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    xlab("nCount_RNA") +
    ylab("nFeature_RNA") +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 500) +
    facet_wrap(~origin) + 
  theme(plot.margin = margin(1,1,1,1, "cm"))

png(file = "QC_RNA_dens.png", width = 5000, height = 5000, res = 300) 
ggarrange(p2,p3,p4,p5,p1,p6, ncol = 2, nrow = 3, 
          labels = c("A","B","C","D","E","F"),
          font.label = list(size = 25), vjust = 1.2)
dev.off()

#' Additionally, we can also see which genes contribute the most to such reads. 
#' We can for instance plot the percentage of counts per gene.

# Compute the relative expression of each gene per spot Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way
# will take a very long time.
par(mar = c(4, 8, 2, 1))
C <- merged@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
png(file = "Box_mostExpressedFeatures.png", width = 5000, height = 3000, res = 400) 
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, 
        xlab = "% total count per spot",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()

saveRDS(merged, file = paste(rds.dir, "OVisium_merged_qc.rds", sep = "/"))

#' prepare supplementary figure 2
A <- plot_grid((P[[1]] + ggtitle("Number of Features or Genes")), p[[1]], 
               ncol = 1, rel_heights = c(2,1))
B <- plot_grid((P[[2]] + ggtitle("Number of UMI Counts")), p[[2]],
               ncol = 1, rel_heights = c(2,1))
C <- plot_grid((P[[3]] + ggtitle("Mitochrondrial Gene Expression %")), p[[3]], 
               ncol = 1, rel_heights = c(2,1))
D <- plot_grid((P[[4]] + ggtitle("Ribosomal Gene Expression %")), p[[4]], 
               ncol = 1, rel_heights = c(2,1))

setwd(qc.dir)
png(file = "QC_spatial_vln.png", 
    width = 8500, height = 11300, res = 300)
ggarrange(A,NULL,B,NULL, NULL, NULL,C,NULL,D, ncol = 3, nrow = 3, 
          labels = c("A","","B","","","","C","","D"), 
          font.label = list(size = 40), 
          widths = c(1, 0.05, 1), 
          heights = c(1, 0.05, 1))
dev.off()
