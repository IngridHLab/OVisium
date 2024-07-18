library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

#' 
count_cell_type_clusters_table <- 
  readr::read_csv("~/OVisium/Deconvolution/Spacexr/Version_1_spaceranger/AGG_BRCA_18/count_cell_type_clusters_table.csv")

             
celltype.count<-count_cell_type_clusters_table %>% column_to_rownames("cluster")
library(tidyverse)  
celltype.count<-celltype.count %>%summarise_all(sum, na.rm = TRUE)
  celltype.count <- t(celltype.count)%>%as.data.frame()%>%rownames_to_column("cell_type")
  celltype.count$BRCA <- "BRCA1/2"
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
               
png(file = "Total_celltype_count_aggregated.png", 
      width = 5000, height = 5000, res = 500)
 print(
   ggplot(count, aes(x=Status, y = Percentage, fill = cell_type, label = Percentage)) + 
    geom_bar(position = "fill", stat = "identity") + 
    scale_fill_manual(values = cell_col) +
     scale_y_continuous(labels = scales::percent) + 
    theme(axis.text = element_text(size=15),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
     geom_text(position = position_fill(vjust=0.5), size = 5) +
    xlab("Status") + 
    ylab("Percentage") +
    ggtitle("Aggregated") 
 )
dev.off()

