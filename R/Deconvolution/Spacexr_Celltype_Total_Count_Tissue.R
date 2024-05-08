library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

#' choose the sample 
#' set working directory to their spacexr output folder
# setwd()
count_cell_type_clusters_table <- readr::read_csv("count_cell_type_clusters_table.csv")

             
celltype.count<-count_cell_type_clusters_table %>% column_to_rownames("cluster")
library(tidyverse)  
celltype.count<-celltype.count %>%summarise_all(sum, na.rm = TRUE)
  celltype.count <- t(celltype.count)%>%as.data.frame()%>%rownames_to_column("cell_type")
  celltype.count$Tissue <- "Fimbrial"
  celltype.count$Tissue <- "Proximal"
  celltype.count$cell_type <- factor(celltype.count$cell_type, 
                                     levels = c("ciliated epithelial cell",
                                                "secretory cell",
                                                "fibroblast",
                                                "myofibroblast cell",
                                                "smooth muscle cell", 
                                                "pericyte",
                                                "blood vessel endothelial cell",
                                                "endothelial cell of lymphatic vessel",
                                                "B cell",
                                                "mature NK T cell",
                                                "mast cell",
                                                "macrophage"))
  
  celltype.count.list <- list()
  celltype.count.list[["PO40_P"]] <- celltype.count
  saveRDS(celltype.count.list, "CellType_count_list_fimbrial_Proximal.rds")
my_list_2 <- mapply('[<-', celltype.count.list, "Patient", value = names(celltype.count.list), SIMPLIFY = F)
cellcount.all <- dplyr::bind_rows(my_list_2) 

#' import percentage data 
  count <- 
    readr::read_csv("~/OVisium/Deconvolution_analysis/Spacexr/Version_1_spaceranger/AGG_BRCA_18/Status_celltype_count_table.csv")
  count$cell_type <- factor(count$cell_type, 
                                     levels = c("ciliated epithelial cell",
                                                "secretory cell",
                                                "fibroblast",
                                                "myofibroblast cell",
                                                "smooth muscle cell", 
                                                "pericyte",
                                                "blood vessel endothelial cell",
                                                "endothelial cell of lymphatic vessel",
                                                "B cell",
                                                "mature NK T cell",
                                                "mast cell",
                                                "macrophage"))
  #'
  #' plot propotion plot
  cell_col <- 
    c("#A6CEE3","#1F78B4","#B8E986", "#7ED321", "#417505", "#FFFF99", "#E31A1C", 
               "#FB9A99", "#CAB2D6", "#6A3D9A", "#FDBF6F", "#B15928")
  names(cell_col) <- c("ciliated epithelial cell",
                       "secretory cell",
                       "fibroblast",
                       "myofibroblast cell",
                       "smooth muscle cell", 
                       "pericyte",
                       "blood vessel endothelial cell",
                       "endothelial cell of lymphatic vessel",
                       "B cell",
                       "mature NK T cell",
                       "mast cell",
                       "macrophage")
               
png(file = "Total_celltype_count_tissue.png", 
      width = 5000, height = 5000, res = 500)
 print(
  ggplot(cellcount.all, aes(x=Tissue, y = V1, fill = cell_type)) +
    geom_bar(position = "fill", stat = "identity") + 
    scale_fill_manual(values = cell_col) +
    scale_y_continuous(labels=scales::percent) + 
    theme(axis.text = element_text(size=20),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("Status") + 
    ylab("Percentage") +
    ggtitle("4 BRCA1 Patient") 
 )
dev.off()

#' boxplot for each cell type to compare fimbrial and proximal
cellcount.all %<>% group_by(Patient) %>% mutate(per = prop.table(V1)*100)

png(file = "Total_celltype_count_per_tissue_split.png", 
    width = 6000, height = 5000, res = 500)
print(
  ggplot(cellcount.all, aes(x=cell_type, y = per, fill = Tissue)) +
    geom_boxplot(alpha = 1) + 
    theme(axis.text = element_text(size=10),
          axis.text.x = element_blank()) +
    xlab("Cell Type") + 
    ylab("Percentage") +
    ggtitle("4 BRCA1 Patients") +
    facet_wrap(~cell_type, scales = "free")
)
dev.off()

png(file = "Total_celltype_count_tissue_split.png", 
    width = 6000, height = 5000, res = 500)
print(
  ggplot(cellcount.all, aes(x=cell_type, y = V1, fill = Tissue)) +
    geom_boxplot(alpha = 1) + 
    theme(axis.text = element_text(size=10),
          axis.text.x = element_blank()) +
    xlab("Cell Type") + 
    ylab("Count") +
    ggtitle("4 BRCA1 Patients") +
    facet_wrap(~cell_type, scales = "free")
)
dev.off()
