#' Determine the resolution during the FindClusters().
#' https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#' https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html

library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))


choose.clustRes <- function(data) {
  vars <- "sample"
  reductions <- c("pca_SCT", "harmony_SCT")

  #' Seurat plot clusters as tree
  pdf("clusterTrees.pdf", onefile = T, width = 6, height = 6)
  print(
    for (i in seq(0.1, 1, by = 0.1)) {
      resolution <- i
      for (red in reductions) {
        npcs <- length(data@reductions[[red]])
        red.ident <- paste(red, "res", resolution, sep = "_")
        #' run clustering on reduction data
        data <- 
          FindNeighbors(data, assay = "SCT", reduction = red, dims = 1:npcs) %>% 
          FindClusters(resolution = resolution) %>% 
          BuildClusterTree(assay = "SCT", reduction = red, dims = 1:npcs)
        
        #' back up the new clusters and save rds
          data[[red.ident]] <- Idents(data)
          
        #' Plot cluster tree
          PlotClusterTree(data, direction = "leftwards")  
          title(red.ident)
      }
    }
  )
dev.off()

  #' perform clustree analysis
  for (red in reductions) {
  red.prefix <- paste(red, "res_", sep = "_")

  p1 <- clustree(data, prefix = red.prefix)
  p2 <- clustree(data, node_colour = "sc3_stability", prefix = red.prefix)
  p3 <- clustree(data, node_colour = "C7", node_colour_aggr = "median", 
                 prefix = red.prefix)

  png(paste0("Clustree_", red, ".png"), width = 6000, height = 3000, res = 300) 
  print(
    plot_grid(p1, p2, p3, ncol = 3, labels = "AUTO", label_size = 12)
  )
  dev.off()
  }

  saveRDS(data, file = paste(rds.dir, "clusters.rds", sep ="/"))
}

dims.file <- list.files(path = rds.dir, 
                        pattern = "OVisium.*dims.rds")
for (rds in dims.file) { 
  data <- readRDS(paste(rds.dir, rds, sep = "/"))
  file.name <- gsub("_dims.rds", "", rds)
  data.type <- gsub("OVisium_", "", file.name)
  
  set.seed(1220)
  dir.create(paste(dims.dir, "Cluster", file.name, sep = "/"))
  setwd(paste(dims.dir, "Cluster", file.name, sep = "/"))
  choose.clustRes(data=data)
  
  setwd(rds.dir)
  file.rename("clusters.rds", paste0(file.name, "_clust", ".rds"))
  }
