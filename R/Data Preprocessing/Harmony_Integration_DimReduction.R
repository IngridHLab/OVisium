#' Perform PCA with the optimal number of PCs and Harmony
#' https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
#' Here we use all variable features instead of most variable features across samples to integrate

library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

#' Run the PCA using the optimal number of PCs and Harmony using library as 
#' integration.
merged.file <- list.files(path = rds.dir, 
                            pattern = "OVisium_SCT_merged.rds")
reductions <- c("pca_SCT", "harmony_SCT")
nl <- NoLegend()

for (i in seq_along(merged.file)) { 
  data <- readRDS(paste(rds.dir, merged.file[i], sep = "/"))
  data@meta.data[["sample"]] <- as.factor(data@meta.data[["sample"]])
  file.name <- gsub(".rds", "", merged.file[i])
  data.type <- gsub("OVisium_", "", file.name)
  
  set.seed(1220)
  out.dir <- paste(dims.dir, toupper(reductions[2]), sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  png(file = paste0("Harmony_convergence_", data.type, ".png"), 
      width = 706, height = 504, res = 100) 
  
    data <- RunHarmony(data, reduction = reductions[1], assay.use = "SCT", 
                       group.by.vars = "sample", plot_convergence = TRUE, 
                       reduction.save = reductions[2])
  
  dev.off()
  
  
  for (red in reductions) {
  out.dir <- paste(dims.dir, toupper(red), sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  
  #' find underlying patterns of transcription
  p1<-DimPlot(data, reduction = red, pt.size = 0.8, 
              group.by = "sample", label = T, label.box = T, repel = F) + 
    xlab(paste0(red, "_1")) + 
    ylab(paste0(red, "_2")) + nl
  p2<-DimPlot(data, reduction = red, pt.size = 0.8, 
              group.by = "slide", label = T, label.box = T, repel = T) + 
    xlab(paste0(red, "_1")) + 
    ylab(paste0(red, "_2")) + nl
  p3<-DimPlot(data, reduction = red, pt.size = 0.8, 
              group.by = "run", label = T, label.box = T, repel = T) + 
    xlab(paste0(red, "_1")) + 
    ylab(paste0(red, "_2")) + nl
  p4<-DimPlot(data, reduction = red, pt.size = 0.8, 
              group.by = "origin", label = T, label.box = T, repel = T) + 
    xlab(paste0(red, "_1")) + 
    ylab(paste0(red, "_2")) + nl
  png(file = paste0("BatchCor_", data.type, ".png"), 
      width = 3000, height = 2800, res = 300) 
  print(
    ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, labels = c("A","B","C","D"))
  )
  dev.off()
  
  #' 2 dimensional reduction for easy visualization
  #' tSNE, UMAP or UWOT on PCA data
  npcs <- length(data@reductions[[red]])
  data <- RunTSNE(data, reduction = red, dims = 1:npcs, 
                  tsne.method = "Rtsne", seed.use = 1, 
                  reduction.name = paste0("tsne_", red))  

  data <- RunUMAP(data, reduction = red, dims = 1:npcs, 
                  umap.method = "umap-learn", metric = "correlation", 
                  min.dist = 0.1, n.neighbors = 30, 
                  reduction.name = paste0("umap_", red))
  
  data <- RunUMAP(data, reduction = red, dims = 1:npcs, 
                  n.components = 3, umap.method = "umap-learn", 
                  metric = "correlation", min.dist = 0.1, n.neighbors = 30, 
                  reduction.name = paste0("umap3d_", red))

  data <- RunUMAP(data, reduction = red, dims = 1:npcs, 
                  umap.method = "uwot", metric = "cosine", 
                  min.dist = 0.01, n.neighbors = 30, 
                  reduction.name = paste0("uwot_", red)) 
  
  data <- RunUMAP(data, reduction = red, dims = 1:npcs, 
                  n.components = 3, umap.method = "uwot", metric = "cosine", 
                  min.dist = 0.01, n.neighbors = 30, 
                  reduction.name = paste0("uwot3d_", red)) 
  
  #' plot 2-dimensional reduction with different annotation in the metadata
  #' BRCA mutation  
  reductions2d<-c("tsne", "umap", "uwot")
  for (red2d in reductions2d) {
    out.dir <- paste(dims.dir, toupper(red), toupper(red2d),  sep = "/")
    dir.create(out.dir, recursive = T)
    setwd(out.dir)
    png(file = paste0("Mutation_", data.type, ".png"), 
        width = 1500, height = 1500, res = 300) 
    print(
      DimPlot(data, group.by = "mutation", label = TRUE, label.box = T, 
              label.size=3, repel = T, reduction = paste0(red2d, "_", red)) + 
        xlab(paste0(red2d, "_1")) + 
        ylab(paste0(red2d, "_2")) + 
        nl
    )
    dev.off() 
    
    #' Tissue origin  
    png(file = paste0("TissueOrigin_", data.type, ".png"), 
        width = 1500, height = 1500, res = 300) 
    print(
      DimPlot(data, group.by = "origin", label = TRUE, label.box = T, 
              label.size=3, repel = T, reduction = paste0(red2d, "_", red)) + 
        xlab(paste0(red2d, "_1")) + 
        ylab(paste0(red2d, "_2")) +  
        nl
    )
    dev.off()
    
    ## Visium slide id  
    png(file = paste0("VisiumSlide_", data.type, ".png"), 
        width = 1500, height = 1500, res = 300) 
    print(
      DimPlot(data, group.by = "slide", label = TRUE, label.box = T, 
              label.size=3, repel = T, reduction = paste0(red2d, "_", red)) + 
        xlab(paste0(red2d, "_1")) + 
        ylab(paste0(red2d, "_2")) +  
        nl
    )
    dev.off()
    
    #' Illumina sequencing run 
    png(file = paste0("SeqRun_", data.type, ".png"), 
        width = 1500, height = 1500, res = 300) 
    print(
      DimPlot(data, group.by = "run", label = TRUE, label.box = T, 
              label.size=3, repel = T, reduction = paste0(red2d, "_", red)) + 
        xlab(paste0(red2d, "_1")) + 
        ylab(paste0(red2d, "_2")) +  
        nl
    )
    dev.off() 
    
    #' sample id 
    png(file = paste0("VisiumSample_", data.type, ".png"), 
        width = 2100, height = 1800, res = 300) 
    print(
      DimPlot(data, group.by = "sample", label = F, 
              reduction = paste0(red2d, "_", red)) + 
        xlab(paste0(red2d, "_1")) + 
        ylab(paste0(red2d, "_2"))
    )
    dev.off() 
    
    #' library id 
    png(file = paste0("Library_", data.type, ".png"), 
        width = 2100, height = 1800, res = 300) 
    print(
      DimPlot(data, group.by = "library", label = F, 
              reduction = paste0(red2d, "_", red)) + 
        xlab(paste0(red2d, "_1")) + 
        ylab(paste0(red2d, "_2"))
    )
    dev.off() 
    
    #' Patient id  
    png(file = paste0("Patient_", data.type, ".png"), 
        width = 2100, height = 1800, res = 300) 
    print(
      DimPlot(data, group.by = "patient", label = F, 
              reduction = paste0(red2d, "_", red)) + 
        xlab(paste0(red2d, "_1")) + 
        ylab(paste0(red2d, "_2"))
    )
    dev.off()
  }
  }
  saveRDS(data, file = paste(rds.dir, paste0(file.name, "_dims.rds"), 
                             sep = "/"))
  }
