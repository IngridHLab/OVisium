library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

#' choose the sample 
#' set working directory to their spacexr output folder

sample.name <- "POxx_x"
setwd(paste(home, "OVisium/Deconvolution_analysis/Spacexr/Version_1_spaceranger/sample.name", sep = "))
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
  celltype.count.list[[sample.name]] <- celltype.count
  saveRDS(celltype.count.list, "CellType_count_list_fimbrial_Proximal.rds")

#' read rds from Deconvolution folder "Spacexr/Version_1_spaceranger"
my_list_2 <- mapply('[<-', celltype.count.list, "Patient", value = names(celltype.count.list), SIMPLIFY = F)
cellcount.all <- dplyr::bind_rows(my_list_2) 

  cellcount.all$cell_type <- factor(cellcount.all$cell_type, 
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

  #' plot proportion plot
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
    ggtitle("4 BRCA1 Patients") 
 )
dev.off()

#' boxplot for each cell type to compare fimbrial and proximal
#' perform pairwise t.test between patient and tissue on each cell type
#' check homoscedasticity - Are the variances homogenous
fimbrial <- subset(cellcount.all, Tissue == "Fimbrial")
proximal <- subset(cellcount.all, Tissue == "Proximal")
var.test(fimbrial$per, proximal$per) 
#' p-value = 0.76 they are not heteroscedasticity 

#' combine PO20_F and PO20_F_a
PO20_F_combined <- rbindlist(celltype.count.list[c("PO20_F", "PO20_F_a")])
PO20_F_combined %<>%
  group_by(cell_type, Tissue) %>%
  summarise(
    n = n(),
    mean = mean(V1),
    sd = sd(V1)
  ) %>%
  ungroup()

celltype.count.list[["PO20_F_combined"]] <-  PO20_F_combined[c(1,4,2)] %>% rename_(., "V1" = 2) %>% as.data.frame()
my_list_3 <- mapply('[<-', celltype.count.list[c(1:4,10,7:9)], "Patient", value = names(celltype.count.list[c(1:4,7:10)]), SIMPLIFY = F)
cellcount.all <- dplyr::bind_rows(my_list_3) 

cellcount.all %<>% group_by(Patient) %>% mutate(per = prop.table(V1)*100)
cellcount.all$Patient <- gsub("_.*", "", cellcount.all$Patient)



png(file = "Total_celltype_count_per_tissue_split.png", 
    width = 6000, height = 5000, res = 500)
print(
  ggboxplot(cellcount.all, x="Tissue", y = "per", fill = "cell_type", 
            palette = cell_col) +
    stat_compare_means(paired = T) +
    stat_compare_means( aes(label = ..p.signif..), 
                        label.x = 1.5, label.y = 40) +
    geom_line(aes(group = Patient)) +
    facet_wrap(~cell_type, scales = "free")
)
dev.off()

png(file = "Total_celltype_count_tissue_split.png", 
    width = 6000, height = 5000, res = 500)
print(
  ggboxplot(cellcount.all, x="Tissue", y = "V1", fill = "cell_type", 
            palette = cell_col) +
    stat_compare_means(paired = T) +
    stat_compare_means( aes(label = ..p.signif..), 
                        label.x = 1.5, label.y = 40) +
    geom_line(aes(group = Patient)) +
    facet_wrap(~cell_type, scales = "free")
)
dev.off()
