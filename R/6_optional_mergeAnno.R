#' Reorder the csv files based on sample name
anno.list <- c("PO2.csv","PO2_2re.csv","PO6.csv","PO7_A1.csv","PO7_2.csv",
               "PO9_B1.csv","PO13.csv","PO13a.csv","PO20.csv","PO20a.csv",
               "PO20_2.csv","PO28.csv","PO35.csv","PO37.csv","PO40.csv",
               "PO40_2.csv","PO41.csv","PO45.csv")

#' Rename the extension of the barcodes to 1-18
#' Remove sample id in the morphology description
#' Merge all 18 samples annotation
Anno <- list()
for (i in seq_along(anno.list)) {
  sample.id <- gsub(".csv", "", anno.list[i])
  anno <- read_csv(file = paste(anno.dir, anno.list[i], sep = "/"), col_names = T)
  anno[["Barcode"]] <- gsub("1", i, anno[["Barcode"]])
  anno[["Graph-based"]] <- gsub(paste0(sample.id, "_"), "", anno[["Graph-based"]])
  Anno[[sample.id]] <- anno
}
morphology.clusters <- bind_rows(Anno)
base::colnames(morphology.clusters)[2] <- "morphology"

#' Write csv file for Loupe Browser
write.table(morphology.clusters, file="~/OVisium/Annotation/morphology.csv", 
            sep = ",", quote = F, row.names = F, col.names = T)