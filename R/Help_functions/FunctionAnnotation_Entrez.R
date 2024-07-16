library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

#' features to EntrezID from aggregated file
features.id <- readRDS(paste(deg.dir, "features_EntrezID.rds", sep = "/"))

#' MSigDB 2023 all gene subsets
#' manually download all gmt files for different gene sets 
release <- "msigdb_v2023.2.Hs_GMTs"
gmt.files <- list.files(paste(deg.dir, "MSigDB", release, sep = "/"), pattern = ".*entrez.gmt")
anno.files <- list()
for (i in seq_along(gmt.files)) {
  name <- str_replace_all(gmt.files[i], ".Hs.entrez.gmt", "")
  gmt <- 
    gmtPathways(paste(deg.dir, "MSigDB", release, gmt.files[i], sep = "/"))
  anno <- reshape2::melt(gmt) %>%
    dplyr::rename("gs_name" = 2, "entrez_gene" = 1) %>%
    dplyr::select(.,2,1) %>%
    na.omit()
  anno$gs_name<-str_remove(anno$gs_name, "^[A-Z]*_")
  anno.files[[name]] <- anno
}

#' Download cell markers:
#' http://yikedaxue.slwshop.cn/download.php

Cell_marker_Human <- read.xlsx(paste(deg.dir, "Cell_marker_Human.xlsx", sep = "/"))
Cell_marker_Seq <- read.xlsx(paste(deg.dir, "Cell_marker_Seq.xlsx", sep = "/"))

#' Human cell markers
cells <- Cell_marker_Human %>%
  dplyr::select(cell_name, GeneID) %>%
  na.omit()
anno.files[["Cell.markers.2.0"]] <- cells

#' Single cell markers
scells <- Cell_marker_Seq %>%
  dplyr::filter(species == "Human") %>%
  dplyr::select(cell_name, GeneID) %>%
  na.omit()
anno.files[["sCell.markers.2.0"]] <- scells

setwd(deg.dir)
saveRDS(anno.files, "FunctionAnnotation_entrez_files.rds")
names <- names(anno.files)
