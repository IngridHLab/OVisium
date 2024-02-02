#' Export projection from Loupe Browser for captures areas that contain 2 samples:
#' Manually select the spots based on the H&E image in the Loupe Browser
#' The projection file contains barcodes and XY coordinates for each spot

home<-"/home/minerva/OVisium"
source(paste(home, "manuscript/gitHub/1_directory.R", sep = "/"))
source(paste(home, "manuscript/gitHub/2_library.R", sep = "/"))

PO7_A1_bar <- 
  read_tsv(gzfile(paste(work.dir, 
                        "Spaceranger/Outs/PO7_A1/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                        sep = "/"))) %>% mutate(1="Barcode")
PO9_B1_bar <- 
  read_tsv(gzfile(paste(work.dir, 
                        "Spaceranger/Outs/PO9_B1/filtered_feature_bc_matrix/barcodes.tsv.gz", 
                        sep = "/"))) %>% mutate(1="Barcode")

#' Import morphology annotation for the joined sections
PO7PO9_anno <- read_csv(paste(anno.dir, "Split_PO7PO9.csv", sep = "/"))
PO9PO7_anno <- read_csv(paste(anno.dir, "Split_PO9PO7.csv", sep = "/"))

#' Split annotation and define spots that are missing annotation to "Unknown"
PO7PO9_anno$Barcode <- gsub(".{2}$", "", PO7PO9_anno$Barcode)
PO9PO7_anno$Barcode <- gsub(".{2}$", "", PO9PO7_anno$Barcode)
PO7_A1_bar$Barcode <- gsub(".{3}$", "", PO7_A1_bar$Barcode)
PO9_B1_bar$Barcode <- gsub(".{3}$", "", PO9_B1_bar$Barcode)

PO7_A1_anno <- right_join(PO7PO9_anno, PO7_A1_bar[1])
PO9_B1_anno <- right_join(PO9PO7_anno, PO9_B1_bar[1])

PO7_A1_anno$Barcode <- paste0(PO7_A1_anno$Barcode, "-1")
PO9_B1_anno$Barcode <- paste0(PO9_B1_anno$Barcode, "-1")

PO7_A1_anno$`Graph-based` <- 
  gsub("PO7PO9", "PO7_A1", PO7_A1_anno$`Graph-based`) %>% 
  replace_na("PO7_A1_Unknown")
PO9_B1_anno$`Graph-based` <- 
  gsub("PO9PO7", "PO9_B1", PO9_B1_anno$`Graph-based`) %>% 
  replace_na("PO9_B1_Unknown")

write.table(PO7_A1_anno, file=paste(anno.dir, "PO7_A1.csv", sep = "/"), 
            sep = ",", quote = F, row.names = F, col.names = T)
write.table(PO9_B1_anno, file=paste(anno.dir, "PO9_B1.csv", sep = "/"), 
            sep = ",", quote = F, row.names = F, col.names = T)
