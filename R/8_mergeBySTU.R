#' Utilize the package developed by the spatialrearch group in sthlm for more:
#' STutility is an R-package with the goal of providing an easy-to-use visualization and analysis tool kit for spatial transcriptomics data. The package was developed by Ludvig Larsson and Joseph Bergenstråhle in the Genomics research group at the Royal Institute of Technology (KTH). The group is positioned at the Science for Life Laboratory (SciLifeLab) in Stockholm, Sweden. Please visit our website https://www.spatialresearch.org/ for more information about the group and our research.

library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

infoTable <- readr::read_csv(paste(rds.dir, "infoTable_20230707.csv", 
                                   sep = "/")) 
infoTable <- data.frame(infoTable)
#' Create Seurat object without standard filters
merged <- InputFromTable(infotable = infoTable, platform =  "Visium")
merged <- LoadImages(merged, time.resolve = T, xdim = 1500)

Idents(merged) <- "sample"
png(file = paste(qc.dir, "Spatial_H&E.png", sep = "/"), width = 6000, 
    height = 5000, res = 300) 
ImagePlot(merged, indices = 1:length(levels(merged)), ncols = 5, 
          method = "raster")
dev.off()

Idents(merged) <- "library"
Idents(merged) <- 
  factor(Idents(merged), 
         levels = c("PO2", "PO2_2re", "PO6", "PO7_A1", "PO7_2", "PO9_B1", 
                    "PO13", "PO13a", "PO20", "PO20a", "PO20_2", "PO28", "PO35", 
                    "PO37", "PO40", "PO40_2", "PO41", "PO45"))
merged[["library"]] <- Idents(merged)

saveRDS(merged, file = paste(rds.dir, "OVisium_merged.rds", sep = "/"))