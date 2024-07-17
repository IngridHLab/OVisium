#' Prepare reference data and Visium data for Seurat label transfer: 
#' Mapping approach that can be used to “anchor” diverse datasets together, 
#' including different types of single cell data (transcriptomic, epigenomic, 
#' and proteomic) and single cell and spatial data.
#' https://satijalab.org/seurat/articles/spatial_vignette.html#integration-with-single-cell-data
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/10_choose_pca_dims.R", sep = "/"), 
       print.eval = TRUE, local = TRUE)

#' load and process scRNAseq reference data and QC in the same way 
#' as the Visium data
#' some issue with the scRNAseq data
setwd(paste(dec.dir, "scRNAseq_Healthy_FT_Hammoud_2022", sep = "/"))
healthy.tubes = readRDS("healthy.tubes.rds")

#' Convert ensembl_ID to gene symbol based on the aggregated feature matrix
#' Collect all genes coded on the mitochondrial genome and ribosomal proteins
#' plot the mito and ribo genes expression level spatially
metafeatures <- 
  healthy.tubes@assays[["RNA"]]@meta.features %>% 
  rownames_to_column() %>% 
  dplyr::rename("ensemble_id"=1)
features <- read.table(paste(sr.dir, "Outs", "AGG_BRCA_18_20231018", 
                             "outs/filtered_feature_bc_matrix/features.tsv.gz", 
                             sep = "/"), sep = "\t") %>% 
  dplyr::rename("ENSEMBL" = 1, "SYMBOL" = 2, "TYPE" = 3)
features$SYMBOL <- make.unique(features$SYMBOL)
mito.features <- metafeatures %>%
dplyr::filter((grepl("^MT-", feature_name)))
ribo.features <- metafeatures %>%
dplyr::filter((grepl("^RP[SL]", feature_name)))


healthy.tubes[["Mito.percent"]] <-
PercentageFeatureSet(healthy.tubes, features = mito.features$ensemble_id)
healthy.tubes[["Ribo.percent"]] <-
PercentageFeatureSet(healthy.tubes, features = ribo.features$ensemble_id)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))

Idents(healthy.tubes) <- "Dataset"
group_by <- "Dataset"
measures <- c("nFeature_RNA", "nCount_RNA", "Mito.percent", "Ribo.percent")
metadata <- Fetch_Meta(healthy.tubes) 
p <- list()
for (i in seq_along(1:length(measures))) {
  p[[i]] <- ggplot(metadata, aes(x = .data[[group_by]],
                           y = .data[[measures[i]]],
                           fill = .data[[group_by]])) +
    scale_fill_manual(values = col_vector) +
    geom_violin() +
    geom_boxplot(alpha = 0.2) +
    theme(legend.position = "none") +
    labs(x = group_by, y = NULL) +
    ggtitle(measures[i]) +
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
}
png(file = "healthy_tubes_QC_RNA_vln.png", 
    width = 5000, height = 2600, res = 300) 
print(
  patchwork::wrap_plots(p)
)
dev.off()

#' Bar plot to show distribution of data
p1 <- ggplot() +
  geom_histogram(data = healthy.tubes[[]], aes(nFeature_RNA), 
                 fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Unique genes per cell")
p2 <- ggplot() +
  geom_histogram(data = healthy.tubes[[]], aes(nCount_RNA), 
                 fill = "red", alpha = 0.7, bins = 50) +
  ggtitle("Total counts per cell")

gene_attr <- 
  data.frame(nUMI = Matrix::rowSums(healthy.tubes@assays$RNA@counts), 
             nSpots = Matrix::rowSums(healthy.tubes@assays$RNA@counts > 0))

p3 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nUMI), 
                 fill = "red", alpha = 0.7, bins = 50) +
  scale_x_log10() +
  ggtitle("Total counts per gene (log10 scale)")
p4 <- ggplot() +
  geom_histogram(data = gene_attr, aes(nSpots), 
                 fill = "red", alpha = 0.7,  bins = 50) +
  ggtitle("Total cells per gene")

png(file = "healthy_tubes_QC_RNA_histo.png", 
    width = 5000, height = 2600, res = 300) 
ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, labels = c("A","B","C","D"))
dev.off()

#' Add number of genes per UMI for each cell to metadata
healthy.tubes$log10GenesPerUMI <- 
  log10(healthy.tubes$nFeature_RNA) / log10(healthy.tubes$nCount_RNA)
metadata <- Fetch_Meta(healthy.tubes) 

#' Visualize the number of cell counts per sample
p1 <-metadata %>% 
  ggplot(aes(x=tissue, fill=tissue)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ylab("nCells") + 
  xlab("")

#' Visualize the number UMIs/transcripts per cell
p2 <- metadata %>% 
  ggplot(aes(color=tissue, x=nCount_RNA, fill= tissue)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

#' Visualize the distribution of genes detected per cell via histogram
p3 <- metadata %>% 
  ggplot(aes(color=tissue, x=nFeature_RNA, fill= tissue)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  ylab("Cell density") +
  geom_vline(xintercept = 500)

#' Visualize the overall complexity of the gene expression by visualizing the 
#' genes detected per UMI (novelty score)
p4 <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = tissue, fill=tissue)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  xlab("FeaturesPerUMICount (log10)") +
  ylab("Cell density") +
  geom_vline(xintercept = 0.8)

#' Visualize the distribution of mitochondrial gene expression detected per cell
p5 <- metadata %>% 
  ggplot(aes(color=tissue, x=Mito.percent, fill=tissue)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  xlab("Mitochondrial %") +
  ylab("Cell density") +
  geom_vline(xintercept = 15)

#' Visualize the correlation between genes detected and number of UMIs and 
#' determine whether strong presence of cells with low numbers of genes/UMIs
p6 <- metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=Mito.percent)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 500) +
  facet_wrap(~tissue)

png(file = "healthy_tubes_QC_RNA_dens.png", width = 5000, height = 4000, 
    res = 300) 
  ggarrange(p2,p3,p4,p5,p1,p6, ncol = 2, nrow = 3, 
            labels = c("A","B","C","D","E","F"))
dev.off()

#' Perform sct transformation of each donor
ref.data.list <- list()
Idents(healthy.tubes) <- "donor_id"
for (i in levels(Idents(healthy.tubes))) { 
  set.seed(1220)
  setwd(rds.dir)
  #' subset each library 
  data <- subset(healthy.tubes, idents = i)
  dataset <- levels(Idents(data))
  
  #' Filter data base on distribution and violin plots to remove out-liner
  data <- subset(data, subset = Ribo.percent > 5 &
                   log10GenesPerUMI > 0.8)
  
  #' SCTransform and regress out difference caused by technical batch.
  data <- SCTransform(data, verbose = TRUE, return.only.var.genes = T)
  
  #' Save under list for later-on merging    
  ref.data.list[[i]] <- data
}
saveRDS(ref.data.list, file = "healthyTubes_SCT_list.rds")

#' perform canonical
set.seed(1220)
setwd(rds.dir)
reference.m <- merge(x = ref.data.list[[1]],
                      y = ref.data.list[2:length(ref.data.list)],
                      project = "Healthy_FT_Hammoud", merge.data = TRUE)
vfeatures <- rownames(reference.m[["SCT"]]@scale.data)
saveRDS(vfeatures, file = "healthyTubes_vfeatures.rds")

#' perform integration of different reference data sets
ref.data.list <- 
  PrepSCTIntegration(ref.data.list, anchor.features = vfeatures)
anchors.ref<- FindIntegrationAnchors(ref.data.list, 
                                     normalization.method = "SCT",
                                     anchor.features = vfeatures)
saveRDS(anchors.ref, file = "healthyTubes_vfeatures_anchors.rds")
reference.m <- IntegrateData(anchors.ref, normalization.method = "SCT")

#' npcs is decided by 'choose.pcs' function 
npcs <- choose.pcs(data=reference.m, assay = "integrated")
file.rename("pca_elbow.png", "pca_elbow_healthyTubes_SCT_vIntegrated.png")
file.rename("pca_quant_elbow.png", 
            "pca_quant_elbow_healthyTubes_SCT_vIntegrated.png")
file.rename("pca_heatmap.png", "pca_heatmap_healthyTubes_SCT_vIntegrated.png")
set.seed(1220)
reference.m <- RunPCA(reference.m, assay = "integrated", npcs = npcs, 
                 reduction.name = "pca_integrated") %>% 
  RunUMAP(assay = "integrated", reduction = "pca_integrated", 
          dims = 1:npcs, umap.method = "uwot", metric = "cosine", 
          min.dist = 0.01, n.neighbors = 30, reduction.name = "uwot_integrated")
saveRDS(reference.m, file = paste(rds.dir, "healthyTubes_SCT_vIntegrated.rds", 
                             sep = "/"))
