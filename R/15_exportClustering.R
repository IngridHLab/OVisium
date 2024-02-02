#' optimal resolution should be based on the clustering and the morphology 
#' will be 0.4 which seperate the FTE and Stromal

library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

#' Export uwot projection and clusters annotation for loupe browser
#' Prepare empty lists
clust.file <- list.files(path = rds.dir, pattern = "OVisium.*.merged_clust.rds")

for (rds in clust.file) { 
  data <- readRDS(paste(rds.dir, rds, sep = "/"))
  file.name <- gsub(".clust.rds", "", rds)
  vars <- "sample"
  reduction <- "harmony_SCT"
  reduction2d <- c("umap_harmony_SCT", "uwot_harmony_SCT")
  
for (res in seq(0.1, 1, by = 0.1)) {
  harmony.ident <- paste0("harmony_SCT_res_", res)
  Idents(data) <- harmony.ident
  out.dir <- paste(anno.dir, file.name, sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)

  #' split the Seurat object
  all.corrected.cluster <- list()
  corrected <- SplitObject(data, split.by = "sample")
  
  #' Edit the barcodes into a format that is compatible with the Loupe Browser
  for (i in seq_along(corrected)) {
    
    #' Extract sample id  
      id <- unique(corrected[[i]]$sample)
  
    #' Extract barcodes
      barcode <- rownames(Embeddings(corrected[[i]], reduction = reduction))
      barcode <- gsub(paste0("-1", "_", i), "", barcode)
      barcode <- paste0(barcode, "-", i)
      
    #' Export new clusters
      clusters <-  Idents(corrected[[i]])
      clusters.data <- cbind("Barcode" = barcode, 
                             data.frame("cluster" = clusters))
      colnames(clusters.data)[2]<-harmony.ident
      all.corrected.cluster[[id]] <- clusters.data
  }
  corrected.clusters <- bind_rows(all.corrected.cluster)
  write.table(corrected.clusters, file=paste0(harmony.ident, ".csv"), sep = ",",
              quote = F, row.names = F, col.names = T)
}
  
  #' Export 2-dimensional reductions 
  for (red in reduction2d) {
    all.corrected.umap <- list()
    for (i in seq_along(corrected)) {
      barcode <- rownames(Embeddings(corrected[[i]], reduction = red))
      barcode <- gsub(paste0("-1", "_", i), "", barcode)
      barcode <- paste0(barcode, "-", i)
      proj <- Embeddings(corrected[[i]], reduction = red)
      umap <- cbind("Barcode" = barcode, data.frame(proj))
      all.corrected.umap[[i]] <- umap
    }
    setwd(out.dir)
    corrected.umap <- bind_rows(all.corrected.umap)
    write.table(corrected.umap, 
                file=paste(red, "csv", sep = "."), 
                sep = ",", quote = F, row.names = F, col.names = T)
  }
}
