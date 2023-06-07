source("/home/minerva/VisiumST/gitHub/library.R")
source("/home/minerva/VisiumST/gitHub/directory.R")

file <- list.files(path = path, pattern = "\\.v2.NMF.harmony.rds$")
file.name <- gsub(".NMF.harmony.rds", "", file)
data <- readRDS(paste(path, file, sep = "/"))

nfcs <- length(data@reductions[["NMF"]])
dims <- 1:nfcs
NMF.reduction <- paste("NMF", nfcs, sep = ".")
data@reductions[["NMF"]]<-data@misc$reductions.backup[[NMF.reduction]]

setwd(dims.dir)
pdf(paste(file.name, NMF.reduction, "clustertrees.pdf", sep = "."), onefile = T, width = 6, height = 6)
for (i in seq(0.1, 1, by = 0.1)) {
  
  set.seed(1220)
  resolution <- i
  ## run clustering on NMF data
  NMF.ident <- paste("NMF", npcs, "SNN", "res", resolution, sep = ".")
  data <- 
    FindNeighbors(data, assay = "SCT", reduction = "NMF", dims = dims) %>% 
    FindClusters(resolution = resolution) %>% 
    BuildClusterTree(assay = "SCT", reduction = "NMF", dims = dims)
  
  ## back up the new clusters and save rds
  data[[NMF.ident]] <- Idents(data)
  
  ## Plot cluster tree
  PlotClusterTree(data, direction = "leftwards")  
  title(NMF.ident)
}
dev.off()

setwd(dims.dir)
head(data@assays$SCT@var.features)
NMF.prefix <- paste("NMF", npcs, "SNN", "res.", sep = ".")
p1 <- clustree(data, prefix = NMF.prefix)
p2 <- clustree(data, node_colour = "sc3_stability", prefix = NMF.prefix)
p3 <- clustree(data, node_colour = "C7", node_colour_aggr = "median", prefix = NMF.prefix)


png(paste(file.name, "NMF.clustree.png", sep = "."), width = 6000, height = 3000, res = 300) 
plot_grid(p1, p2, p3, ncol = 3)
dev.off()


## Save data
setwd(path)
saveRDS(data, file = paste(file.name, "clusterings", "rds", sep = "."))