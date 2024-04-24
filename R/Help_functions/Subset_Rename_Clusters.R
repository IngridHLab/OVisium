#' Exclude cluster 3
Idents(data) <- cluster.ident
data.sub <- SubsetSTData(data, idents = c(0:3, 5:11))

#' Merge fimbrial and otherside tissues
Idents(data.sub) <- "origin"
data.sub <- RenameIdents(data.sub, "Otherside" = "Fimbrial")
Idents(data.sub) <- factor(Idents(data.sub),
                       levels = c("Fimbrial", "Proximal"))
data.sub[["Tissue_origin"]] <- Idents(data.sub)

#' Create a new category of the clusters based on morphology
Idents(data.sub) <-cluster.ident
data.sub <- RenameIdents(data.sub, "0" = "0_FTE")
data.sub <- RenameIdents(data.sub, "1" = "1_FTE")
data.sub <- RenameIdents(data.sub, "2" = "2_Stroma")
data.sub <- RenameIdents(data.sub, "3" = "3_Mix")
data.sub <- RenameIdents(data.sub, "5" = "5_Stroma")
data.sub <- RenameIdents(data.sub, "6" = "6_Stroma")
data.sub <- RenameIdents(data.sub, "7" = "7_Stroma")
data.sub <- RenameIdents(data.sub, "8" = "8_Stroma")
data.sub <- RenameIdents(data.sub, "9" = "9_Stroma")
data.sub <- RenameIdents(data.sub, "10" = "10_Stroma")
data.sub <- RenameIdents(data.sub, "11" = "11_Stroma")
Idents(data.sub) <- factor(Idents(data.sub),
                       levels = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                                  "5_Stroma", "6_Stroma", "7_Stroma",
                                  "8_Stroma", "9_Stroma", "10_Stroma",
                                  "11_Stroma"))
data.sub[["Visium_clusters"]] <- Idents(data.sub)
