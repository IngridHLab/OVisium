#' https://www.10xgenomics.com/resources/analysis-guides/integrating-10x-visium-and-chromium-data-with-r
#' https://raw.githack.com/dmcable/spacexr/master/vignettes/visium_full_regions.html
#' Spacexr/Robust cell-type decomposition (RCTD):
#' Deconvolution approach that uses a reference-based probabilistic model to 
#' resolve cell types from a single spot containing a mixture of cell types, 
#' infers the cell type proportions with a maximum-likelihood estimation, 
#' and projects them onto a spatial map of cell types.
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

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
sc_counts <- healthy.tubes@assays[["SCT"]]@counts
sc_umis <- healthy.tubes@meta.data[,c(34,35)]
sc_umis <- setNames(sc_umis[[1]], rownames(sc_umis))
Idents(healthy.tubes) <- "cell_type"
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
healthy.tubes[["cell_type"]] <- Idents(healthy.tubes)
cell_types<-healthy.tubes@active.ident
SCreference <- Reference(sc_counts, cell_types, sc_umis)
gc()

#' prepare color codes for cell types
cell_col <- 
  c("#A6CEE3","#1F78B4","#B8E986", "#7ED321", "#417505", "#FFFF99", "#E31A1C", 
             "#FB9A99", "#CAB2D6", "#6A3D9A", "#FDBF6F", "#B15928")
names(cell_col) <- levels(cell_types)             
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))   

#' load Visium data 
file.name <- "OVisium_SCT_merged" 
cluster.ident <- "harmony_SCT_res_0.6"
setwd(rds.dir)
data <- readRDS(data, file = paste0(file.name, "_deg.rds"))
Idents(data) <- "Visium_clusters"

#' Extract barcodes
barcode <- base::colnames(data)
barcode <- gsub("-1_", "-", barcode)

#' Export new clusters
clusters <-  Idents(data)
visium.clusters <- cbind("barcode" = barcode, data.frame("cluster" = clusters, 
                                                         row.names = NULL))

#' load processed spatial data
Idents(data) <- "library"
ids <- levels(data)

#' prepare data list 
VisiumData.list <- list()
RCTD.list <- list()
agg.stack.plot.list <- list()
agg.fill.plot.list <- list()
count.stack.plot.list <- list()
count.fill.plot.list <- list()

#' Perform RCTD on individual samples
for (i in seq_along(ids)) {
  Sys.setenv("OPENBLAS_NUM_THREADS"=4)

  #' Prepare Visium data
  #' Read in the tissue_positions_list.csv file
  #' Optional: To match the orientation of the tissue on the slide with spacexr 
  #' result plots, we will flip and mirror the coordinates:
  setwd(paste(sr.dir, "Outs", ids[[i]], "outs", sep = "/"))
  vis_coords <- read.csv("spatial/tissue_positions.csv")
  vis_coords[,c(3,4)] <- vis_coords[,c(4,3)]
  vis_coords[,4] <- vis_coords[,4]*-1
  write.table(vis_coords, "spatial/tissue_positions.csv", quote=FALSE, 
              row.names=FALSE, col.names = T, sep=',')
  
  #' Load Visium data directly from Space Ranger output directory
  VisiumData <- read.VisiumSpatialRNA("./")
  barcodes <- colnames(VisiumData@counts)
  barcodes <- gsub("1", i, barcodes)
  clusters.barcodes <- subset(visium.clusters, barcode %in% barcodes)
  clusters.barcodes$barcode <- gsub(i, "1", clusters.barcodes$barcode)

  #' convert gene name to ensemble id
  features <- 
    readr::read_tsv("filtered_feature_bc_matrix/features.tsv.gz", 
                    col_names = FALSE) %>% 
    dplyr::rename("ENSEMBL" = 1, "SYMBOL" = 2, "TYPE" = 3)
  VisiumData@counts@Dimnames[1]<-list(features$ENSEMBL)
  VisiumData.list[[i]] <- VisiumData
  
  #' Create the output directory in your working directory
  out.dir <- paste(dec.dir, "Spacexr", ids[[i]], sep = "/") 
  dir.create(out.dir, recursive = TRUE) 
  setwd(out.dir)
  plot_puck_continuous(VisiumData, barcodes = clusters.barcodes$barcode,
                       plot_val = VisiumData@nUMI, 
                       size = 1, 
                       ylimit = c(0,round(quantile(VisiumData@nUMI,0.9))),
                       title = 'plot of nUMI',
                       ylim = c(-100,10))
  ggsave("cell_type_nUMI.png", height=2500, width=2800, units='px', dpi=400)

  #' Create and run RCTD (robust cell type decomposition)
  myRCTD <- create.RCTD(VisiumData, SCreference, max_cores = 4)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  RCTD.list[[i]] <- myRCTD 

  #' Exploring the full mode results
  #' Create variables from the myRCTD object to plot results
  #' Weights for each cell type per barcode
  #' Normalize per spot weights so cell type probabilities sum to 1 for each spot
  weights <- myRCTD@results$weights 
  norm_weights <- normalize_weights(weights) 
  
  #' List of cell types
  cell_type_names <- colnames(norm_weights) 
  #' matrix of the number of spots and 12 cell types
  dim(norm_weights)
  
  #' For each spot barcode (row), you can see the normalized weight for each 
  #' cell type (column)
  weights_table <- as.data.frame(weights[, cell_type_names]) %>% 
    rownames_to_column() %>% dplyr::rename("barcode"=1)
  weights_table <- inner_join(clusters.barcodes, weights_table)
  weights_table$barcode <- gsub("1", i, weights_table$barcode)
  write.table(weights_table, "weight_clusters_table.csv", 
              quote=FALSE, row.names=F, sep=',')
  norm_weights_table <- as.data.frame(norm_weights[, cell_type_names]) %>% 
    rownames_to_column() %>% dplyr::rename("barcode"=1)
  norm_weights_table <- inner_join(clusters.barcodes, norm_weights_table)
  norm_weights_table$barcode <- gsub("1", i, norm_weights_table$barcode)
  write.table(norm_weights_table, "norm_weight_clusters_table.csv", 
              quote=FALSE, row.names=F, sep=',')
  
  agg_cell_types_table_list <- list()
  count_cell_types_table_list <- list()
  for (c in seq_along(levels(clusters.barcodes$cluster))) {
    cluster <- levels(clusters.barcodes$cluster)[c]
    cluster.barcodes <- 
      clusters.barcodes[which(clusters.barcodes$cluster == cluster), ]$barcode
    if (length(cluster.barcodes) < 2)
    {cluster.barcodes <- rep(cluster.barcodes, times = 2)}
    #' aggregate pixel count of each cell type without filtering based on the 
    #' cell type weight
    
    agg_cell_types <- aggregate_cell_types(myRCTD,
                                           barcodes = cluster.barcodes,
                                           doublet_mode = F) 
    agg_cell_types_table <- as.data.frame(agg_cell_types) 
    colnames(agg_cell_types_table)[1] <- cluster
    agg_cell_types_table_list[[cluster]] <- agg_cell_types_table
    
    #' count cell type based on cell type weight threshold 1
    count_cell_types <- count_cell_types(myRCTD,
                                         barcodes = cluster.barcodes,
                                         cell_types = cell_type_names,
                                         cell_type_threshold = 0,
                                         doublet_mode = F)
    count_cell_types_table <- as.data.frame(count_cell_types) 
    colnames(count_cell_types_table)[1] <- cluster
    count_cell_types_table_list[[cluster]] <- count_cell_types_table
  }
  
  #' Plot aggregate results
  agg <- bind_cols(agg_cell_types_table_list) %>% rownames_to_column() %>% 
    dplyr::rename("cell_type"=1)
  agg$cell_type <- factor(agg$cell_type, levels = cell_type_names)
  magg <- reshape2::melt(agg) %>% dplyr::rename("cluster"=2)
  p1 <- ggplot(magg, aes(cluster, value, fill = cell_type)) + 
    geom_bar(position = "stack", stat = "identity") + 
    scale_fill_manual(values = cell_col[cell_type_names]) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("OVisium clusters") + 
    ylab("Aggregated pixel counts") +
    ggtitle(ids[[i]]) 
  agg.stack.plot.list[[i]] <- p1
  
  p2 <- ggplot(magg, aes(cluster, value, fill = cell_type)) + 
    geom_bar(position = "fill", stat = "identity") + 
    scale_fill_manual(values = cell_col[cell_type_names]) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("OVisium clusters") + 
    ylab("Aggregated pixel proportions") +
    ggtitle(ids[[i]]) 
  agg.fill.plot.list[[i]] <- p2
  
  #' Plot count results
  count <- bind_cols(count_cell_types_table_list) %>% rownames_to_column() %>% 
    dplyr::rename("cell_type"=1)
  count$cell_type <- factor(count$cell_type, levels = cell_type_names)
  mcount <- reshape2::melt(count) %>% dplyr::rename("cluster"=2)
  p3 <- ggplot(mcount, aes(cluster, value, fill = cell_type)) + 
    geom_bar(position = "stack", stat = "identity") + 
    scale_fill_manual(values = cell_col[cell_type_names]) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("OVisium clusters") + 
    ylab("Weight filtered pixel counts") +
    ggtitle(ids[[i]]) 
  count.stack.plot.list[[i]] <- p3
  
  p4 <- ggplot(mcount, aes(cluster, value, fill = cell_type)) + 
    geom_bar(position = "fill", stat = "identity") + 
    scale_fill_manual(values = cell_col[cell_type_names]) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("OVisium clusters") + 
    ylab("Weight filtered pixel proportions") +
    ggtitle(ids[[i]]) 
  count.fill.plot.list[[i]] <- p4
  
  #' Export results
  tagg <-transpose(agg[,-1])
  colnames(tagg) <- agg$cell_type
  rownames(tagg) <- colnames(agg)[2:12]
  tagg <-tagg[, cell_type_names] %>% rownames_to_column() %>% dplyr::rename("cluster"=1)
  write.table(tagg, "agg_cell_type_clusters_table.csv", 
              quote=FALSE, col.names=T, row.names=F, sep=',')
  tcount <-transpose(count[,-1])
  colnames(tcount) <- count$cell_type
  rownames(tcount) <- colnames(count)[2:12]
  tcount <-tcount[, cell_type_names] %>% rownames_to_column() %>% dplyr::rename("cluster"=1)
  write.table(tcount, "count_cell_type_clusters_table.csv", 
              quote=FALSE, col.names=T, row.names=F, sep=',')

#' Barplot of confident counts of each cell type in 'full_mode'. 
#' saved as 'cell_type_occur.pdf' 
  plot_cond_occur(cell_type_names, "./", norm_weights, myRCTD@spatialRNA) + 
    theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1)) +
    scale_x_discrete(limits=cell_type_names) +
    scale_fill_manual(breaks = cell_type_names, values = cell_col[cell_type_names])
  ggsave("cell_type_occur.png", height=3000, width=3500, units = "px", dpi=400, 
         bg = "white")

#' plot normalized weight or binary weight for each cell type spatially
  for (t in 1:length(cell_type_names)) {
    #' Plot normalized weight
    plot_puck_continuous(myRCTD@spatialRNA, 
                         clusters.barcodes$barcode, 
                         norm_weights[,cell_type_names[t]], 
                         title =cell_type_names[t], 
                         size=0.5,
                         ylim = c(-100, 10))
    ggsave(paste(cell_type_names[t],'weights.png', sep='_'), 
           height=2500, width=2800, units='px', dpi=400)
    
    #' Plot binary weight 
    #' create a binary weight for each cell type using median weight as cutoff 
    #' above 75th percentile  as 1
    #' below median called as 0
    cell_type_weight <- myRCTD@results$weights[, cell_type_names[t]]
    med <- median(cell_type_weight)
    q3 <- quantile(cell_type_weight, probs=.75)
    cell_type_weight[cell_type_weight>=q3] <- 1
    cell_type_weight[cell_type_weight<med] <- 0 
    plot_puck_continuous(myRCTD@spatialRNA, 
                         clusters.barcodes$barcode, 
                         cell_type_weight, 
                         title = cell_type_names[t], 
                         size = 0.5,
                         ylim = c(-100, 10))
    ggsave(paste(cell_type_names[t], "binary_weights.png", sep = "_"), 
           height=2500, width=2800, units='px', dpi=400)
  }
}
gc()

#' Note! 
#' Perform RCTD on the aggregated data
Sys.setenv("OPENBLAS_NUM_THREADS"=4)
setwd(paste(sr.dir, "Outs", "AGG_BRCA_18_20231018", "outs", sep = "/"))
VisiumData <- read.VisiumSpatialRNA("./")
barcodes <- colnames(VisiumData@counts)
clusters.barcodes <- subset(visium.clusters, barcode %in% barcodes)
features <- 
  readr::read_tsv("filtered_feature_bc_matrix/features.tsv.gz", 
                  col_names = FALSE) %>% 
  dplyr::rename("ENSEMBL" = 1, "SYMBOL" = 2, "TYPE" = 3)
VisiumData@counts@Dimnames[1]<-list(features$ENSEMBL)
VisiumData.list[["agg"]] <- VisiumData

#' Create the output directory in your working directory
out.dir <- paste(dec.dir, "Spacexr", "AGG_BRCA_18", sep = "/") 
dir.create(out.dir, recursive = TRUE) 
setwd(out.dir)
plot_puck_continuous(VisiumData, barcodes = clusters.barcodes$barcode,
                     plot_val = VisiumData@nUMI, 
                     size = 1, 
                     ylimit = c(0,round(quantile(VisiumData@nUMI,0.9))),
                     title = 'plot of nUMI',
                     xlim = c(0, 100))
ggsave("cell_type_nUMI.png", height=2500, width=2800, units='px', dpi=400)

#' Create and run RCTD (robust cell type decomposition)
myRCTD <- create.RCTD(VisiumData, SCreference, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
RCTD.list[["agg"]] <- myRCTD 

#' For each spot barcode (row), you can see the normalized weight for each 
#' cell type (column)
weights <- myRCTD@results$weights 
weights_table <- as.data.frame(weights[,cell_type_names]) %>% 
  rownames_to_column() %>% dplyr::rename("barcode"=1)
weights_table <- inner_join(clusters.barcodes, weights_table)
write.table(weights_table, "weight_clusters_table.csv", 
            quote=FALSE, row.names=F, sep=',')
norm_weights <- normalize_weights(weights) 
norm_weights_table <- as.data.frame(norm_weights[,cell_type_names]) %>% 
  rownames_to_column() %>% dplyr::rename("barcode"=1)
norm_weights_table <- inner_join(clusters.barcodes, norm_weights_table)
write.table(norm_weights_table, "norm_weight_clusters_table.csv", 
            quote=FALSE, row.names=F, sep=',')
cell_type_names <- colnames(norm_weights) 

agg_cell_types_table_list <- list()
count_cell_types_table_list <- list()
for (c in seq_along(levels(clusters.barcodes$cluster))) {
  cluster <- levels(clusters.barcodes$cluster)[c]
  cluster.barcodes <- 
    clusters.barcodes[which(clusters.barcodes$cluster == cluster), ]$barcode
  #' aggregate pixel count of each cell type without filtering based on the 
  #' cell type weight
  agg_cell_types <- aggregate_cell_types(myRCTD,
                                         barcodes = cluster.barcodes,
                                         doublet_mode = F) 
  agg_cell_types_table <- as.data.frame(agg_cell_types) 
  colnames(agg_cell_types_table)[1] <- cluster
  agg_cell_types_table_list[[cluster]] <- agg_cell_types_table
  
  #' count cell type based on cell type weight threshold 0.95
  count_cell_types <- count_cell_types(myRCTD,
                                       barcodes = cluster.barcodes,
                                       cell_types = cell_type_names,
                                       cell_type_threshold = 0,
                                       doublet_mode = F)
  count_cell_types_table <- as.data.frame(count_cell_types) 
  colnames(count_cell_types_table)[1] <- cluster
  count_cell_types_table_list[[cluster]] <- count_cell_types_table
}

#' Plot aggregate results
agg <- bind_cols(agg_cell_types_table_list) %>% rownames_to_column() %>% 
  dplyr::rename("cell_type"=1)
agg$cell_type <- factor(agg$cell_type, levels = cell_type_names)
magg <- reshape2::melt(agg) %>% dplyr::rename("cluster"=2)
p1 <- ggplot(magg, aes(cluster, value, fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = cell_col[cell_type_names]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab("OVisium clusters") + 
  ylab("Aggregated pixel counts") +
  ggtitle("Aggregated")  
agg.stack.plot.list[["agg"]] <- p1

p2 <- ggplot(magg, aes(cluster, value, fill = cell_type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = cell_col[cell_type_names]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab("OVisium clusters") + 
  ylab("Aggregated pixel proportions") +
  ggtitle("Aggregated")  
agg.fill.plot.list[["agg"]] <- p2

#' Plot count results
count <- bind_cols(count_cell_types_table_list) %>% rownames_to_column() %>% 
  dplyr::rename("cell_type"=1)
count$cell_type <- factor(count$cell_type, levels = cell_type_names)
mcount <- reshape2::melt(count) %>% dplyr::rename("cluster"=2)
p3 <- ggplot(mcount, aes(cluster, value, fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = cell_col[cell_type_names]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab("OVisium clusters") + 
  ylab("Weight filtered pixel counts") +
  ggtitle("Aggregated")  
count.stack.plot.list[["agg"]] <- p3

p4 <- ggplot(mcount, aes(cluster, value, fill = cell_type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = cell_col[cell_type_names]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab("OVisium clusters") + 
  ylab("Weight filtered pixel proportions") +
  ggtitle("Aggregated") 
count.fill.plot.list[["agg"]] <- p4

#' Export results
tagg <-transpose(agg[,-1])
colnames(tagg) <- agg$cell_type
rownames(tagg) <- colnames(agg)[2:12]
tagg <-tagg[, cell_type_names] %>% rownames_to_column() %>% dplyr::rename("cluster"=1)
write.table(tagg, "agg_cell_type_clusters_table.csv", 
            quote=FALSE, col.names=T, row.names=F, sep=',')
tcount <-transpose(count[,-1])
colnames(tcount) <- count$cell_type
rownames(tcount) <- colnames(count)[2:12]
tcount <-tcount[, cell_type_names] %>% rownames_to_column() %>% dplyr::rename("cluster"=1)
write.table(tcount, "count_cell_type_clusters_table.csv", 
            quote=FALSE, col.names=T, row.names=F, sep=',')

#' Barplot of confident counts of each cell type in 'full_mode'. 
#' saved as 'cell_type_occur.pdf' 
plot_cond_occur(cell_type_names, "./", norm_weights, myRCTD@spatialRNA) + 
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1)) +
  scale_x_discrete(limits=cell_type_names) +
  scale_fill_manual(breaks = cell_type_names, values = cell_col[cell_type_names])
ggsave("cell_type_occur.png", height=3000, width=3500, units='px', dpi=400, 
       bg = "white")

#' plot normalized weight or binary weight for each cell type spatially
for (t in 1:length(cell_type_names)) {
  #' Plot normalized weight
  plot_puck_continuous(myRCTD@spatialRNA, 
                       clusters.barcodes$barcode, 
                       norm_weights[,cell_type_names[t]], 
                       title =cell_type_names[t], 
                       size=0.5,
                       xlim = c(0, 100))
  ggsave(paste(cell_type_names[t],'weights.png', sep='_'), 
         height=2500, width=2800, units='px', dpi=400)
  
  #' Plot binary weight 
  #' create a binary weight for each cell type using median weight as cutoff 
  #' above 75th percentile  as 1
  #' below median called as 0
  cell_type_weight <- myRCTD@results$weights[, cell_type_names[t]]
  med <- median(cell_type_weight)
  q3 <- quantile(cell_type_weight, probs=.75)
  cell_type_weight[cell_type_weight>=q3] <- 1
  cell_type_weight[cell_type_weight<med] <- 0 
  plot_puck_continuous(myRCTD@spatialRNA, 
                       clusters.barcodes$barcode, 
                       cell_type_weight, 
                       title = cell_type_names[t], 
                       size = 0.5,
                       xlim = c(0, 100))
  ggsave(paste(cell_type_names[t], "binary_weights.png", sep = "_"), 
         height=2500, width=2800, units='px', dpi=400)
}
gc()

#' save and plot all results in one plot
setwd(paste(dec.dir, "Spacexr", sep = "/"))
saveRDS(VisiumData.list, file = "VisiumData_list.rds")
saveRDS(RCTD.list, file = "RCTD_list.rds")
saveRDS(agg.stack.plot.list, file = "agg.stack.plot.list.rds")
saveRDS(agg.fill.plot.list, file = "agg.fill.plot.list.rds")
saveRDS(count.stack.plot.list, file = "count.stack.plot.list.rds")
saveRDS(count.fill.plot.list, file = "count.fill.plot.list.rds")

#' plot figure 7
for (i in 1:18) {
  agg.stack.plot.list[[i]] <- agg.stack.plot.list[[i]] + 
    theme(legend.position = "none", axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + labs(title = i, x = "", y = "")
  agg.fill.plot.list[[i]] <- agg.fill.plot.list[[i]] + 
    theme(legend.position = "none", axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + labs(title = i, x = "", y = "")
  count.stack.plot.list[[i]] <- count.stack.plot.list[[i]] + 
    theme(legend.position = "none", axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + labs(title = i, x = "", y = "")
  count.fill.plot.list[[i]] <- count.fill.plot.list[[i]] + 
    theme(legend.position = "none", axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + labs(title = i, x = "", y = "")
}

#' plot only those have pair sample from the proximal end
png(file = "aggStack_cluster_pair_aggregated.png", 
    width = 8000, height = 7000, res = 500) 
print(gridExtra::grid.arrange(grobs=agg.stack.plot.list[c(1:2,4:5,9:11,15:16,19)], 
                              widths = c(1, 1, 1, 1),
                              height = c(1,1,1),
                              layout_matrix = cbind(c(1:2,19,19),
                                                    c(4:5,19,19),
                                                    c(9:11, NA),
                                                    c(15:16,NA,NA))) 
)
dev.off()

png(file = "aggFill_cluster_pair_aggregated.png", 
    width = 8000, height = 7000, res = 500) 
print(gridExtra::grid.arrange(grobs=agg.fill.plot.list[c(1:2,4:5,9:11,15:16,19)], 
                              widths = c(1, 1, 1, 1),
                              height = c(1,1,1),
                              layout_matrix = cbind(c(1:2,19,19),
                                                    c(4:5,19,19),
                                                    c(9:11, NA),
                                                    c(15:16,NA,NA))) 
)
dev.off()



png(file = "countStack_cluster_pair_aggregated.png", 
    width = 8000, height = 7000, res = 500) 
print(gridExtra::grid.arrange(grobs=count.stack.plot.list[c(1:2,4:5,9:11,15:16,19)], 
                              widths = c(1, 1, 1, 1),
                              height = c(1,1,1),
                              layout_matrix = cbind(c(1:2,19,19),
                                                    c(4:5,19,19),
                                                    c(9:11, NA),
                                                    c(15:16,NA,NA))) 
)
dev.off()

png(file = "countFill_cluster_pair_aggregated.png", 
    width = 8000, height = 7000, res = 500) 
print(gridExtra::grid.arrange(grobs=count.fill.plot.list[c(1:2,4:5,9:11,15:16,19)], 
                              widths = c(1, 1, 1, 1),
                              height = c(1,1,1),
                              layout_matrix = cbind(c(1:2,19,19),
                                                    c(4:5,19,19),
                                                    c(9:11, NA),
                                                    c(15:16,NA,NA))) 
)
dev.off()

png(file = "countFill_cluster_pair_case1.png", 
    width = 2000, height = 3000, res = 500) 
print(gridExtra::grid.arrange(grobs=count.fill.plot.list[c(1:2)], nrow = 2) 
)
dev.off()

#' plot all samples and aggregated data
png(file = "aggStack_cluster_all_aggregate.png", 
    width = 8000, height = 9000, res = 500) 
print(gridExtra::grid.arrange(grobs=agg.stack.plot.list, 
                              widths = c(1, 1, 1, 2.5),
                              layout_matrix = rbind(c(1:3, NA),
                                                    c(4:6, NA),
                                                    c(7:9, NA),
                                                    c(10:12, 19),
                                                    c(13:15, 19),
                                                    c(16:18, NA)))
)
dev.off()
      
png(file = "aggFill_cluster_all_aggregate.png", 
    width = 8000, height = 9000, res = 500) 
print(gridExtra::grid.arrange(grobs=agg.fill.plot.list, 
                              widths = c(1, 1, 1, 2.5),
                              layout_matrix = rbind(c(1:3, NA),
                                                    c(4:6, NA),
                                                    c(7:9, NA),
                                                    c(10:12, 19),
                                                    c(13:15, 19),
                                                    c(16:18, NA)))
)
dev.off()

png(file = "countStack_cluster_all_aggregate.png", 
    width = 8000, height = 9000, res = 500) 
print(gridExtra::grid.arrange(grobs=count.stack.plot.list, 
                              widths = c(1, 1, 1, 2.5),
                              layout_matrix = rbind(c(1:3, NA),
                                                    c(4:6, NA),
                                                    c(7:9, NA),
                                                    c(10:12, 19),
                                                    c(13:15, 19),
                                                    c(16:18, NA)))
)
dev.off()

png(file = "countFill_cluster_all_aggregate.png", 
    width = 8000, height = 9000, res = 500) 
print(gridExtra::grid.arrange(grobs=count.fill.plot.list, 
                              widths = c(1, 1, 1, 2.5),
                              layout_matrix = rbind(c(1:3, NA),
                                                    c(4:6, NA),
                                                    c(7:9, NA),
                                                    c(10:12, 19),
                                                    c(13:15, 19),
                                                    c(16:18, NA)))
)
dev.off()

#' plot piechart for sample 4 and 5
#' coordinates
#' data for sample 4

for (i in 1:18) {
  out.dir <- paste(dec.dir, "Spacexr", ids[[i]], sep = "/") 
  dir.create(out.dir, recursive = TRUE) 
  setwd(out.dir)
  #' coordinates in Visium data
  coords <- data.frame(VisiumData.list[[i]]@coords) %>% rownames_to_column() %>% 
  dplyr::rename("barcode" = 1)
  #' weights in the RCTD result
  weights <- RCTD.list[[i]]@results[["weights"]]
  norm_weights <- normalize_weights(weights)
  cell_type_names <- RCTD.list[[i]]@cell_type_info[["info"]][[2]]
  norm <- as.data.frame(norm_weights[, cell_type_names]) %>% 
    rownames_to_column() %>% dplyr::rename("barcode"=1)
  #' merge and filter away spots without weight (no signal)
  norm <- left_join(norm, coords)
  
#' select spot with cluster annotation
norm$barcode <- gsub("1", i, norm$barcode)
norm <- inner_join(visium.clusters, norm)

spot.pie.plot.list <- list()
# Basic piechart for all spots
for (s in seq_along(1:length(norm$cluster))) {
  # Create Data
  plot.data <-data.frame(
    group = names(norm[s, 3:14]),
    value = (norm[s, 3:14]) %>% unlist() %>% unname()
  )
 p <-ggplot(plot.data, aes(x="", y=value, fill=group)) +
   geom_bar(stat="identity") +
   scale_fill_manual(values = cell_col[cell_type_names]) +
   coord_polar("y", start=0) +
   theme_void() +
   theme(legend.position = "none")
 spot.pie.plot.list[[s]] <- p
}

#' prepare layout matrix for spatial plotting
plot.coords <- data.frame(x=-norm$y, y=norm$x, spot=as.numeric(rownames(norm)))
plot.matrix <- xtabs(spot ~ x + y, data = plot.coords)
plot.matrix[plot.matrix==0]<-NA
rownames(plot.matrix) <-NULL
colnames(plot.matrix) <-NULL
attr(plot.matrix, "class") <- NULL
attr(plot.matrix, "call") <- NULL


png(file = "SpatialPie_cell_type_occur.png", 
    width = ncol(plot.matrix)*100, height = nrow(plot.matrix)*170, res=800) 
gridExtra::grid.arrange(grobs=spot.pie.plot.list,
                              layout_matrix = cbind(plot.matrix))
dev.off()
}

