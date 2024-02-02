#' Determine number of dims in the dataset.
#' https://github.com/hbctraining/scRNA-seq/blob/master/lessons/elbow_plot_metric.md

library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))



choose.pcs <- function(data, assay) {
  set.seed(1220)
  dir <- paste(dims.dir, paste0("PCA_", toupper(assay)), sep = "/")
  dir.create(dir, recursive = T)
  setwd(dir)
  
  #' perform PCA on high dimensions 
  #' 50 is enough and no need for 100 like the spacerange
  data <- RunPCA(data, assay = assay)
  data@misc$reductions.backup[["pca.50"]] <- data[["pca"]]

  #' JackStraw can not be performed on SCT data and it takes very long time. 
  #' In this Rscript, we only run the elbowplot to estimate the number of dims.
  png("pca_elbow.png", width = 1200, height = 1000, 
    res = 300)
  print(
    ElbowPlot(data, ndims = 50)
    )
  dev.off()

  #' Quantitative elbowplot:
  #' Determine percent of variation associated with each PC
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100

  #' Calculate cumulative sum of variation for each PC
  cumu <- cumsum(pct)

  #' Determine which PC exhibits cumulative variation greater than 90% and 
  #' variation associated with the PC as less than 5%
  co1 <- which(cumu > 90 & pct < 5)[1]

  #' Determine the difference between variation of a PC and subsequent PC 
  #' less than 0.1%
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
            decreasing = T)[1] + 1

  #' Minimum of the two calculation
  threshold <- min(co1, co2)
  
  #' Create a data frame with values
  plot_df <- data.frame(pct = pct, 
                      Cumu.Variation = cumu,
                     PCs = 1:length(pct))

  #' Elbow plot to visualize
  png("pca_quant_elbow.png", width = 1200, height = 1000, res = 300)
  print(
    ggplot(
      plot_df, 
      aes(Cumu.Variation, pct, label = PCs, color = PCs >= threshold)) + 
      geom_text(size = 3) + 
      geom_vline(xintercept = 90, color = "red") + 
      geom_hline(yintercept = min(pct[pct > 5]), color = "red") +
      theme_classic() + 
      theme(axis.text=element_text(color='black')) +
      theme(strip.background = element_blank()) +
      theme(legend.position = c(0.8, 0.8))
    )
  dev.off()

  #' check the first 20 PCs to confirm the signal below/above threshold
  npcs <- threshold-1 
  png("pca_heatmap.png", width = 4000, height = 5000, res=400)
  print(
    DimHeatmap(data, dims = 1:npcs, reduction = "pca", ncol = 4, 
               assays = assay) + 
      theme(plot.margin = margin(1,1,1,1, "cm"))
  )
  dev.off()
  
  return(npcs)
  }




