source("/home/minerva/VisiumST/gitHub/library.R")
source("/home/minerva/VisiumST/gitHub/directory.R")

set.seed(1220)
dir.create(paste(dec.dir, "Isthmus", sep ="/")) 
dir.create(paste(dec.dir, "Ampulla", sep ="/")) 
dir.create(paste(dec.dir, "Fimbria", sep ="/")) 
dir.create(paste(dec.dir, "FT", sep ="/"))
ref.name <- "scRNA_Hammoud_2022"

## load and process scRNAseq reference data
healthy.tubes <- readRDS("/home/minerva/VisiumST/scRNAseq_Normal_FT_Hammoud_2022/healthy.tubes.rds")
fimbria.markers <- FindMarkers(healthy.tubes, assay = "RNA", ident.1 = c("Fimbria1", "Fimbria2", "Fimbria3"), group.by = "Dataset", verbose = T, min.pct = 0.01, min.diff.pct = 0.2) %>% rownames_to_column("gene_id") %>%
  left_join(y = unique(annotations[, c("gene_id", "gene_name", "description")]), by = c("gene_id" = "gene_id"))
write.table(fimbria.markers, file = paste(data.dir, "STU/DEA", ref.name, "fimbria.markers.csv", sep = "/"), sep =",", quote = F, row.names = F, col.names = T)

ciliated.epithelial <- readRDS("/home/minerva/VisiumST/scRNAseq_Normal_FT_Hammoud_2022/ciliated.epithelial.rds")
ciliated.fimbria.markers <- FindMarkers(ciliated.epithelial, assay = "RNA", ident.1 = c("Fimbria1", "Fimbria2", "Fimbria3"), group.by = "Dataset", verbose = T, min.pct = 0.01, min.diff.pct = 0.2) %>% rownames_to_column("gene_id") %>%
  left_join(y = unique(annotations[, c("gene_id", "gene_name", "description")]), by = c("gene_id" = "gene_id"))
write.table(ciliated.fimbria.markers, file = paste(data.dir, "STU/DEA", ref.name, "ciliated.fimbria.markers.csv", sep = "/"), sep =",", quote = F, row.names = F, col.names = T)

non.ciliated.epithelial <- readRDS("/home/minerva/VisiumST/scRNAseq_Normal_FT_Hammoud_2022/non.ciliated.epithelial.rds")
non.ciliated.fimbria.markers <- FindMarkers(non.ciliated.epithelial, assay = "RNA", ident.1 = c("Fimbria1", "Fimbria2", "Fimbria3"), group.by = "Dataset", verbose = T, min.pct = 0.01, min.diff.pct = 0.2) %>% rownames_to_column("gene_id") %>%
  left_join(y = unique(annotations[, c("gene_id", "gene_name", "description")]), by = c("gene_id" = "gene_id"))
write.table(non.ciliated.fimbria.markers, file = paste(data.dir, "STU/DEA", ref.name, "non.ciliated.fimbria.markers.csv", sep = "/"), sep =",", quote = F, row.names = F, col.names = T)

stroma <- readRDS("/home/minerva/VisiumST/scRNAseq_Normal_FT_Hammoud_2022/stroma.rds")
stroma.fimbria.markers <- FindMarkers(stroma, assay = "RNA", ident.1 = c("Fimbria1", "Fimbria2", "Fimbria3"), group.by = "Dataset", verbose = T, min.pct = 0.01, min.diff.pct = 0.2) %>% rownames_to_column("gene_id") %>%
  left_join(y = unique(annotations[, c("gene_id", "gene_name", "description")]), by = c("gene_id" = "gene_id"))
write.table(stroma.fimbria.markers, file = paste(data.dir, "STU/DEA", ref.name, "stroma.fimbria.markers.csv", sep = "/"), sep =",", quote = F, row.names = F, col.names = T)

Idents(healthy.tubes)<-"Dataset"
Isthmus.idents <- c("Isthmus1", "Isthmus2", "Isthmus3")
Isthmus <- subset(healthy.tubes, idents = Isthmus.idents )
Mix <- SubsetSTData(data, idents = 3) 
stroma <- SubsetSTData(data, idents = c(1, 4, 5))

## prepare counts from scRNA data
sc_counts <- healthy.tubes@assays[["RNA"]]@counts
# Prepare nUMI object
sc_umis <- healthy.tubes@meta.data[,c(19,20)]
sc_umis <- setNames(sc_umis[[2]], rownames(sc_umis))
# prepare celltype matrix
Idents(healthy.tubes)<-"cell_type"
cell_types<-healthy.tubes@active.ident
cell_types <- as.factor(cell_types) 
# Create single cell reference object
SCreference <- Reference(sc_counts, cell_types, sc_umis)

#### Load merged visium data directly from Space Ranger output directory
# Read in the tissue_positions_list.csv file 
# Optional: To match the orientation of the tissue on the slide with spacexr result plots, we will flip and mirror the coordinates:
merged_coords <- read.csv(paste("/home/minerva/VisiumST/scRNAseq_Normal_FT_Hammoud_2022/outs", "aggr_tissue_positions_list.csv", sep = "/"), header=F)
colnames(merged_coords)[1]<-"barcodes"
colnames(merged_coords)[2]<-"in_tissue"
colnames(merged_coords)[3]<-"x"
colnames(merged_coords)[4]<-"y"
colnames(merged_coords)[5]<-"pxl_row_in_fullres"
colnames(merged_coords)[6]<-"pxl_col_in_fullres"
write.table(merged_coords, paste("/home/minerva/VisiumST/scRNAseq_Normal_FT_Hammoud_2022/outs/spatial", "tissue_positions_list.csv", sep = '/'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',') 

merged<-read.VisiumSpatialRNA("/home/minerva/VisiumST/scRNAseq_Normal_FT_Hammoud_2022/outs")
features.merged <- read.table("/home/minerva/VisiumST/scRNAseq_Normal_FT_Hammoud_2022/outs/filtered_feature_bc_matrix/features.tsv", sep="\t") %>% dplyr::rename("id"=1, "name"=2, "type"=3)
merged@counts@Dimnames[1]<-list(features.merged$id)

## Differential regulate one gene
marker<-read.csv("/home/minerva/VisiumST/gitHub/Data/Seurat_Cell_Markers.csv")
gene_list <- data.frame(marker$Hammound_2022[1:40] )%>% dplyr::rename("name"=1)
gene_list <- left_join(gene_list, features.merged)

## extract coordinate from merged data 
  ##  Create and run RCTD (robust cell type decomposition)
  Merged_RCTD <- create.RCTD(merged, SCreference, max_cores = 4)
  Merged_RCTD <- run.RCTD(Merged_RCTD, doublet_mode = "full")
  saveRDS(Merged_RCTD, file.path(paste(dec.dir, "Merged", "Merged_RCTD_visium_full.rds", sep ="/")))
  
  ## Exploring the full mode results
  # Create variables from the myRCTD object to plot results
  barcodes <- colnames(Merged_RCTD@spatialRNA@counts) # list of spatial barcodes
  weights <- Merged_RCTD@results$weights # Weights for each cell type per barcode
  weights_table <- as.data.frame(weights[,cell_type_names]) %>% rownames_to_column() %>% dplyr::rename("barcode"=1)
  write.table(weights_table, paste0(paste(dec.dir, "Merged", sep ="/"), "/", 'weight_table.csv'), quote=FALSE, row.names=F, sep=',')
  
  # Normalize per spot weights so cell type probabilities sum to 1 for each spot
  norm_weights <- normalize_weights(weights) 
  cell_type_names<-colnames(norm_weights) # List of cell types
  # dim(norm_weights) # matrix of 2,264 spots and 21 cell types
  # For each spot barcode (row), you can see the normalized weight for each cell type (column)
  print(head(norm_weights[,cell_type_names]))
  Merged_RCTD@results$norm_weights <- norm_weights
  dim(norm_weights) # matrix of the number of spots and 12 cell types
  # For each spot barcode (row), you can see the normalized weight for each cell type (column)
  norm_weights_table <- as.data.frame(norm_weights[,cell_type_names]) %>% rownames_to_column() %>% dplyr::rename("barcode"=1)
  write.table(norm_weights_table, paste0(paste(dec.dir, "Merged", sep ="/"), "/", 'norm_weight_table.csv'), quote=FALSE, row.names=F, sep=',')
  
  # aggregate pixel count of each cell type without filtering based on the cell type weight
  agg_cell_type <- aggregate_cell_types(Merged_RCTD, barcodes = barcodes, doublet_mode = F) 
  agg_cell_type_table <- as.data.frame(agg_cell_type) %>% rownames_to_column() %>% dplyr::rename("cell_type"=1, "count"=2)
  write.table(agg_cell_type_table, paste0(paste(dec.dir, "Merged", sep ="/"), "/", 'agg_cell_type_table.csv'), quote=FALSE, row.names=F, sep=',')

  # count cell type based on cell type weight threshold 1
  count_cell_type <- count_cell_types(Merged_RCTD, barcodes = barcodes, cell_types = cell_type_names, cell_type_threshold = 11, weight_threshold = 1, doublet_mode = F)
  count_cell_type_table <- as.data.frame(count_cell_type) %>% rownames_to_column() %>% dplyr::rename("cell_type"=1, "count"=2)
  write.table(count_cell_type_table, paste0(paste(dec.dir, "Merged", sep ="/"), "/", 'count_cell_type_table.csv'), quote=FALSE, row.names=F, sep=',')
  
  # Plot cell type probabilities (normalized) per spot (red = 1, blue = 0 probability)
  # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 'results/cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, paste(dec.dir, "Merged", sep ="/"), norm_weights, Merged_RCTD@spatialRNA) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0(paste(dec.dir, "Merged", sep ="/"), "/", "cell_type_occur.pdf"), height=8, width=8, units='in', dpi=300)
  
  # Plots the confident weights for each cell type as in full_mode (saved as 'results/cell_type_weights_unthreshold.pdf')
  plot_weights(cell_type_names, Merged_RCTD@spatialRNA, paste(dec.dir, "Merged", sep ="/"), norm_weights) 
  # Plots all weights for each cell type as in full_mode. (saved as 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, Merged_RCTD@spatialRNA, paste(dec.dir, "Merged", sep ="/"), norm_weights) 
  
  
  ## Running CSIDE
  # create a binary variable for specific cell type enrichment using median normalized weight as cutoff 
  ciliated_weight <- Merged_RCTD@results$norm_weights[, 'ciliated epithelial cell']
  ciliated_med <- median(ciliated_weight)
  ciliated_weight[ciliated_weight>=ciliated_med] <- 1 #above or equal median called as 1
  ciliated_weight[ciliated_weight<ciliated_med] <- 0 #below median called as 0
  plot_puck_continuous(Merged_RCTD@spatialRNA, rownames(Merged_RCTD@spatialRNA@coords), ciliated_weight)
  ggsave(paste(dec.dir, "Merged", "ciliated_weight_binary.pdf", sep ="/"), height=5, width=5, units='in', dpi=500)
  
  aggregate_cell_types(Merged_RCTD, barcodes = barcodes)
  Merged_RCTD <- run.CSIDE.single(Merged_RCTD, ciliated_weight, cell_types = cell_type_names, doublet_mode = F, fdr = 0.01, cell_type_threshold = 50)
  print(Merged_RCTD@de_results$sig_gene_list)  
  
  
  # create a binary variable for specific cell type enrichment using median weight as cutoff 
  ciliated_weight <- Merged_RCTD@results$weights[, 'ciliated epithelial cell']
  ciliated_med <- median(ciliated_weight)
  ciliated_weight[ciliated_weight>=ciliated_med] <- 1 #above or equal median called as 1
  ciliated_weight[ciliated_weight<ciliated_med] <- 0 #below median called as 0
  plot_puck_continuous(Merged_RCTD@spatialRNA, rownames(Merged_RCTD@spatialRNA@coords), ciliated_weight)
  ggsave(paste(dec.dir, "Merged", "ciliated_binary_weight.pdf", sep ="/"), height=5, width=5, units='in', dpi=500)
  
  #DEgenes
  Merged_RCTD@config$max_cores <- 1
  cell_type_names <- c('secretory cell', 'myofibroblast cell')
  rownames(Merged_RCTD@originalSpatialRNA@counts) # at least 1000 genes 
  aggregate_cell_types(Merged_RCTD, barcodes = barcodes)
  Merged_RCTD <- run.CSIDE.single(Merged_RCTD, ciliated_weight, cell_types = cell_type_names, doublet_mode = F, fdr = 0.01)
  print(Merged_RCTD@de_results$sig_gene_list)  

  ## Differentially down regulate one gene
  marker<-read.csv("/home/minerva/VisiumST/gitHub/Data/Seurat_Cell_Markers.csv")
  gene_list <- data.frame(marker$Hammound_2022[1:40] )%>% dplyr::rename("name"=1)
  gene_list <- left_join(gene_list, features)
  myRCTD@originalSpatialRNA@counts <- myRCTD@spatialRNA@counts[gene_list$id,]
  class_num <- rep(0, length(barcodes)); names(class_num) <- barcodes
  class_num[region_middle] <- 1; class_num[region_right] <- 2 
  plot_class(myRCTD@spatialRNA, barcodes, factor(class_num),  title ='plot of regions') # plot regions
  
  
#### prepare visium data from different clusters
FTE <- readRDS("/home/minerva/VisiumST/gitHub/Data/SCT_RDS/SCT.Merged.v1.FTE.rds")
Stromal <- readRDS("/home/minerva/VisiumST/gitHub/Data/SCT_RDS/SCT.Merged.v1.Stromal.rds")
Mixed <- readRDS("/home/minerva/VisiumST/gitHub/Data/SCT_RDS/SCT.Merged.v1.Mixed.rds")



#### extract coordinate, counts and nUMI from FTE clusters
Idents(FTE)<-"orig.ident" 
names <- levels(Idents(FTE))

FTE_coords <- list()
for (id in names) {
  vis_coords <- FTE@images[[id]]@coordinates[,c(2,3)]
  colnames(vis_coords)[1]<-"x"
  colnames(vis_coords)[2]<-"y"
  FTE_coords[[id]]<-vis_coords
}

  FTE_coords <- bind_rows(FTE_coords)
  write.table(FTE_coords, paste(dec.dir, "FTE", 'FTE_tissue_positions.csv', sep ="/"), quote=FALSE, row.names=T, sep=',')
  
  FTE_counts <- FTE@assays[["Spatial"]]@counts
  FTE_counts@Dimnames[[1]] <- features.merged$id
  
  FTE_uUMI <- FTE@meta.data[2]
  FTE_uUMI <- setNames(FTE_uUMI[[1]], rownames(FTE_uUMI))
  
  FTE_vis <- SpatialRNA(FTE_coords, FTE_counts, FTE_uUMI)
  
  ##  Create and run RCTD (robust cell type decomposition)
  FTE_RCTD <- create.RCTD(FTE_vis, SCreference, max_cores = 4)
  FTE_RCTD <- run.RCTD(FTE_RCTD, doublet_mode = "full")
  saveRDS(FTE_RCTD, file.path(paste(dec.dir, "FTE", "FTE_RCTD_visium_full.rds", sep ="/")))
  
  ## Exploring the full mode results
  # Create variables from the myRCTD object to plot results
  barcodes <- colnames(FTE_RCTD@spatialRNA@counts) # list of spatial barcodes
  weights <- FTE_RCTD@results$weights # Weights for each cell type per barcode
  
  # Normalize per spot weights so cell type probabilities sum to 1 for each spot
  norm_weights <- normalize_weights(weights) 
  cell_type_names<-colnames(norm_weights) # List of cell types
  # dim(norm_weights) # matrix of 2,264 spots and 21 cell types
  # For each spot barcode (row), you can see the normalized weight for each cell type (column)
  print(head(norm_weights[,cell_type_names]))
  
  # Plot cell type probabilities (normalized) per spot (red = 1, blue = 0 probability)
  # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 'results/cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, paste(dec.dir, "FTE", sep ="/"), norm_weights, FTE_RCTD@spatialRNA) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste(dec.dir, "FTE", "cell_type_occur.pdf", sep ="/"), height=8, width=8, units='in', dpi=300)
  
  # Plots the confident weights for each cell type as in full_mode (saved as 'results/cell_type_weights_unthreshold.pdf')
  plot_weights(cell_type_names, FTE_RCTD@spatialRNA, paste(dec.dir, "FTE", sep ="/"), norm_weights) 
  # Plots all weights for each cell type as in full_mode. (saved as 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, FTE_RCTD@spatialRNA, paste(dec.dir, "FTE", sep ="/"), norm_weights) 
  
  
## Running CSIDE
  # create a binary variable for specific cell type enrichment using median weight as cutoff 
  ciliated_weight <- FTE_RCTD@results$weights[, 'ciliated epithelial cell']
  ciliated_med <- median(ciliated_weight)
  ciliated_weight[ciliated_weight>=ciliated_med] <- 1 #above or equal median called as 1
  ciliated_weight[ciliated_weight<ciliated_med] <- 0 #below median called as 0
  plot_puck_continuous(FTE_RCTD@spatialRNA, rownames(FTE_RCTD@spatialRNA@coords), ciliated_weight)
  ggsave(paste(dec.dir, "FTE", "ciliated_binary_weight.pdf", sep ="/"), height=5, width=5, units='in', dpi=500)
  
  #DEgenes
  FTE_RCTD@config$max_cores <- 1
  cell_type_names <- c('secretory cell', 'myofibroblast cell')
  rownames(FTE_RCTD@originalSpatialRNA@counts) # at least 1000 genes 
  aggregate_cell_types(FTE_RCTD, barcodes = barcodes)
  FTE_RCTD <- run.CSIDE.single(FTE_RCTD, ciliated_weight, cell_types = cell_type_names, doublet_mode = F, fdr = 0.01)
  print(FTE_RCTD@de_results$sig_gene_list)
  
  
#### extract coordinate, counts and nUMI from stromal clusters
  Idents(Stromal)<-"orig.ident" 
  names <- levels(Idents(Stromal))
  
  Stromal_coords <- list()
  for (id in names) {
    vis_coords <- Stromal@images[[id]]@coordinates[,c(2,3)]
    colnames(vis_coords)[1]<-"x"
    colnames(vis_coords)[2]<-"y"
    Stromal_coords[[id]]<-vis_coords
  }
  
  Stromal_coords <- bind_rows(Stromal_coords)
  write.table(Stromal_coords, paste(dec.dir, "Stromal", 'Stromal_tissue_positions.csv', sep ="/"), quote=FALSE, row.names=T, sep=',')
  
  Stromal_counts <- Stromal@assays[["Spatial"]]@counts
  Stromal_counts@Dimnames[[1]] <- features.merged$id
  
  Stromal_uUMI <- Stromal@meta.data[2]
  Stromal_uUMI <- setNames(Stromal_uUMI[[1]], rownames(Stromal_uUMI))
  
  Stromal_vis <- SpatialRNA(Stromal_coords, Stromal_counts, Stromal_uUMI)
  
  ##  Create and run RCTD (robust cell type decomposition)
  Stromal_RCTD <- create.RCTD(Stromal_vis, SCreference, max_cores = 4)
  Stromal_RCTD <- run.RCTD(Stromal_RCTD, doublet_mode = "full")
  saveRDS(Stromal_RCTD, file.path(paste(dec.dir, "Stromal", "Stromal_RCTD_visium_full.rds", sep ="/")))
  
  ## Exploring the full mode results
  # Create variables from the myRCTD object to plot results
  barcodes <- colnames(Stromal_RCTD@spatialRNA@counts) # list of spatial barcodes
  weights <- Stromal_RCTD@results$weights # Weights for each cell type per barcode
  
  # Normalize per spot weights so cell type probabilities sum to 1 for each spot
  norm_weights <- normalize_weights(weights) 
  cell_type_names<-colnames(norm_weights) # List of cell types
  # dim(norm_weights) # matrix of 2,264 spots and 21 cell types
  # For each spot barcode (row), you can see the normalized weight for each cell type (column)
  print(head(norm_weights[,cell_type_names]))
  
  # Plot cell type probabilities (normalized) per spot (red = 1, blue = 0 probability)
  # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 'results/cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, paste(dec.dir, "Stromal", sep ="/"), norm_weights, Stromal_RCTD@spatialRNA) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste(dec.dir, "Stromal", "cell_type_occur.pdf", sep ="/"), height=8, width=8, units='in', dpi=300)
  
  # Plots the confident weights for each cell type as in full_mode (saved as 'results/cell_type_weights_unthreshold.pdf')
  plot_weights(cell_type_names, Stromal_RCTD@spatialRNA, paste(dec.dir, "Stromal", sep ="/"), norm_weights) 
  # Plots all weights for each cell type as in full_mode. (saved as 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, Stromal_RCTD@spatialRNA, paste(dec.dir, "Stromal", sep ="/"), norm_weights) 
  
  
  ## Running CSIDE
  # create a binary variable for specific cell type enrichment using median weight as cutoff 
  smm_weight <- Stromal_RCTD@results$weights[, 'smooth muscle cell']
  smm_med <- median(smm_weight)
  smm_weight[smm_weight>=smm_med] <- 1 #above or equal median called as 1
  smm_weight[smm_weight<smm_med] <- 0 #below median called as 0
  plot_puck_continuous(Stromal_RCTD@spatialRNA, rownames(Stromal_RCTD@spatialRNA@coords), ciliated_weight)
  ggsave(paste(dec.dir, "Stromal", "stromal_binary_weight.pdf", sep ="/"), height=5, width=5, units='in', dpi=500)
  
  #DEgenes
  Stromal_RCTD@config$max_cores <- 1
  cell_type_names <- c('secretory cell', 'myofibroblast cell')
  Stromal_RCTD <- run.CSIDE.single(Stromal_RCTD, ciliated_weight, cell_types = cell_type_names, doublet_mode = F, cell_type_threshold = 10, fdr = 0.05,  weight_threshold = 0.1)
  print(Stromal_RCTD@de_results$sig_gene_list)