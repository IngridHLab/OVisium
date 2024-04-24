#' Perform different gene expression analysis using Seurat functions:
#' FindAllMarkers(): Gene expression markers for all identity classes;
#' FindConservedMarkers(): Finds markers that are conserved between the groups;
#' FindMarker(): Gene expression markers of identity classes 

library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/17_geneAnnotations.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/get_marker.R", sep = "/"))

file.name <- "OVisium_SCT_merged" 
merged <- readRDS(paste(rds.dir, paste0(file.name, "_clust.rds"), sep = "/"))
cluster.ident <- "harmony_SCT_res_0.6"
Idents(merged) <- cluster.ident

#' Critical step to normalize the counts before DEA.
#' Given a merged object with multiple SCT models, this function uses minimum 
#' of the median UMI (calculated using the raw UMI counts) of individual 
#' objects to reverse the individual SCT regression model using minimum of 
#' median UMI of all as the sequencing depth covariate. The counts slot of the SCT 
#' assay is replaced with re-corrected counts and the data slot is replaced with 
#' log1p of re-corrected counts.
#' Exclude cluster 3
data <- SubsetSTData(merged, idents = c(0:3, 5:11))

#' Merge fimbrial and otherside tissues
Idents(data) <- "origin"
data <- RenameIdents(data, "Otherside" = "Fimbrial")
Idents(data) <- factor(Idents(data),
                       levels = c("Fimbrial", "Proximal"))
data[["Tissue_origin"]] <- Idents(data)

#' Create a new category of the clusters based on morphology
Idents(data) <-cluster.ident
data <- RenameIdents(data, "0" = "0_FTE")
data <- RenameIdents(data, "1" = "1_FTE")
data <- RenameIdents(data, "2" = "2_Stroma")
data <- RenameIdents(data, "3" = "3_Mix")
data <- RenameIdents(data, "5" = "5_Stroma")
data <- RenameIdents(data, "6" = "6_Stroma")
data <- RenameIdents(data, "7" = "7_Stroma")
data <- RenameIdents(data, "8" = "8_Stroma")
data <- RenameIdents(data, "9" = "9_Stroma")
data <- RenameIdents(data, "10" = "10_Stroma")
data <- RenameIdents(data, "11" = "11_Stroma")
Idents(data) <- factor(Idents(data),
                       levels = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                  "5_Stroma", "6_Stroma", "7_Stroma",
                                  "8_Stroma", "9_Stroma", "10_Stroma",
                                  "11_Stroma"))
data[["Visium_clusters"]] <- Idents(data)
setwd(rds.dir)
saveRDS(data, file = paste0(file.name, "_11_clust.rds"))

data <- PrepSCTFindMarkers(data, assay = "SCT", verbose = T)
setwd(rds.dir)
saveRDS(data, file = paste0(file.name, "_deg.rds"))
