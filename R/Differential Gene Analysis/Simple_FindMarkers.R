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

#' select features
vfeatures <- VariableFeatures(data)
gene_attr <- 
  data.frame(nUMI = Matrix::rowSums(data@assays[["SCT"]]@counts), 
             nSpots = Matrix::rowSums(data@assays[["SCT"]]@counts > 0))

hfeatures <- row.names(gene_attr[which(gene_attr$nUMI>500),]) 
vhfeatures <- unique(vfeatures[vfeatures%in%hfeatures])
feature.list <- list(vfeatures, vhfeatures)
names(feature.list) <- c("Variable_features","Variable_a500UMI_features")
test <- c("MAST", "wilcox")


#' Find all markers of all graph-based clusters in PCA-Harmony
for(i in seq_along(feature.list)) {
  for(t in seq_along(test)) {
features <- feature.list[[i]]
out.dir <- paste(deg.dir, file.name, cluster.ident, 
                 names(feature.list[i]), test[t], sep = "/") 
dir.create(out.dir, recursive = T)
setwd(out.dir)
Idents(data) <- "Visium_clusters"

#' Find differential expressed markers between 2 clusters
#' FTE clusters 
#' min.diff.pct at 0.25 has no gene found, set 10 instead

for (x in seq_along(levels(data))) {
  gene.list <- list()
  for (y in seq_along(1:11)) {
    if (y != x) {
      markers <- FindMarkers(data, 
                             assay = "SCT",
                             features = features,
                             ident.1 = levels(data)[[x]], 
                             ident.2 = levels(data)[[y]],
                             min.pct = -Inf,
                             min.diff.pct = -Inf,
                             logfc.threshold = -Inf,
                             test.use = test[t]) %>% 
        rownames_to_column() %>% 
        dplyr::rename("gene" = 1) %>% 
        mutate(rank=-log10(p_val)*sign(avg_log2FC)) 
        colnames(markers)[2] <- paste0("p_val_", levels(data)[[y]])
        colnames(markers)[3] <- paste0("avg_log2FC_", levels(data)[[y]])
        colnames(markers)[6] <- paste0("p_val_adj_", levels(data)[[y]])
        colnames(markers)[7] <- paste0("rank_", levels(data)[[y]])
      gene.list[[y]] <- markers[c(1,3,2,6,7)]
    } else {
      next
    }
  }
  gene.list[sapply(gene.list, is.null)] <- NULL
  save(gene.list, file = paste0("Cluster.",levels(data)[[x]],".genelist.rds"))
  all.markers <- Reduce(function(x, y) full_join(x, y, by = "gene"), gene.list)
  write.table(all.markers, 
            file = paste0("Cluster.",levels(data)[[x]],".findmarkers.csv"),
            sep =",", quote = F, row.names = F, col.names = T)
}
  }
}
