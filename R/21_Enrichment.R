#' Prepare standard gene list for compareCluster 
#' perform geneset enrichment analysis
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

#' import all annotation files
setwd(deg.dir)
features.id <- readRDS("features_EntrezID.rds")
anno.files <- readRDS("FunctionAnnotation_entrez_files.rds")
names <- names(anno.files)

#' Enrichment analysis between clusters
 plot_compareCluster <- function(file, anno, dir, ...) {
 	
   setwd(dir)
   name <- str_replace_all(file, ".markers.filt.csv", "")
   
   universe <- features.id$ENTREZID
   markers <- readr::read_csv(file, col_names = T)
   
   #' Convert gene name to entrezid
   common <- inner_join(markers[,c("cluster", "gene", "avg_log2FC")], 
                        features.id, by = c("gene" = "SYMBOL"))
   
   #' Check genes that don't have EntrezID
   warning(paste(setdiff(markers$gene, common$gene), 
                 "EntrezID is missing; ", sep = " "))
   
  
   common$level[common$avg_log2FC >= log2(1.5)] <- "Up" # FC > 1.5
   common$level[common$avg_log2FC <= log2(0.75)] <- "Down" # FC < 0.75
   common <- na.omit(common)

   g1 <- c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
           "5_Stroma", "6_Stroma", "7_Stroma",
           "8_Stroma", "9_Stroma", "10_Stroma",
           "11_Stroma")
   g2 <- c("Up", "Down", "No Change")
   
   out.dir <- paste(dir, "Enrichment", sep = "/")
   dir.create(out.dir, recursive = T)
   setwd(out.dir)
   result <- compareCluster(ENTREZID~cluster+level, 
                            data         = common, 
                            universe     = universe,
                            ...)
   if(length(result) == 0) {
     warning(paste0(name, " has no result from compareCluster!"))
     next
   } else {
   result@compareClusterResult$cluster <- 
     factor(result@compareClusterResult$cluster, levels = g1)
   result@compareClusterResult$level <- 
     factor(result@compareClusterResult$level, levels = g2 )
   saveRDS(result, 
           file = paste("CompareCluster", anno, name, "result.rds", sep = "_"))
   
   if (length(unique(result@compareClusterResult$Description))>100){
     png(file = paste("CompareCluster", anno, name, "dotplot.png", sep = "_"), 
         width = 2400, height = 3500, res=300)
     print(
       dotplot(result, x="cluster", showCategory = 10, includeAll = F) + 
         facet_grid(~level) + 
         theme(axis.text.x = element_text(hjust = 1, angle = 60), 
               axis.text.y = element_text(size = 7)) +
         ggtitle(paste("Enrichment of", anno, "for", name, sep = " "))
     )
     dev.off()
   } else {
     png(file = paste("CompareCluster", anno, name, "dotplot.png", sep = "_"), 
         width = 2400, height = 3000, res=300)
     print(
       dotplot(result, x="cluster", showCategory = 10, includeAll = F) + 
         facet_grid(~level) + 
         theme(axis.text.x = element_text(hjust = 1, angle = 60), 
               axis.text.y = element_text(size = 8)) +
         ggtitle(paste("Enrichment of", anno, "for", name, sep = " "))
     )
     dev.off()
     }
   }
 }
 
#' Enrichment analysis on each cluster on positive regulated genes 
 plot_enrichBar <- function(file, anno, assay, dir, ...) {
 	
   setwd(dir)
   name <- str_replace_all(file, ".markers.filt.top100.csv", "")
   
   universe <- features.id$ENTREZID
   top.markers <- readr::read_csv(file, col_names = TRUE)
   
   #' Convert gene name to entrezid
   common <- inner_join(top.markers[,c("cluster", "gene")], 
                        features.id, by = c("gene" = "SYMBOL"))
   
   #' Check genes that don't have EntrezID
   warning(paste(setdiff(top.markers$gene, common$gene), 
                 "skipped; ", sep = " "))
   
   #' reorder the clusters
   common$cluster <- 
     factor(common$cluster, 
            levels = c("0_FTE","1_FTE", "3_Mix", "2_Stroma", 
                       "5_Stroma", "6_Stroma", "7_Stroma",
                       "8_Stroma", "9_Stroma", "10_Stroma",
                       "11_Stroma"))
   
   result <- list()
   out.dir <- paste(dir, "Enrichment", sep = "/")
   dir.create(out.dir, recursive = T)
   setwd(out.dir)
   for (c in 1:length(levels(common$cluster))) {
    
     c.name <- levels(common$cluster)[c]
     common.cluster <- subset(common, cluster == c.name)
     
     #' Enrichment analysis on individual cluster
     if (assay == "enrichGO") {
       eer <- enrichGO(gene          = common.cluster$ENTREZID,
                       universe      = universe, 
                       OrgDb         = org.Hs.eg.db,
                       ...)
     } else {
       eer <- enricher(gene          = common.cluster$ENTREZID,
                       universe      = universe,
                       ...)
     }
     
     if (dim(eer)[1] == 0 || length(eer) == 0) {
       warning(paste0(c.name, " has no over-representation on any ontology!"))
       next     
       } else {
       result[[c.name]] <- eer
     }
   }
   
   if (length(result) == 0) {
     warning(paste0(name, " has no over-representation result for any of the 
                    clusters!"))
     next
   } else {
   saveRDS(result, file = paste("PosEnrich", anno, name, "result.rds", sep = "_"))
   
   plot_list <- list()
   for (i in names(result)) {
       barplot(result[[i]], order = TRUE, showCategory = 10)  +
       theme(axis.text = element_text(size = 4)) +
       ggtitle(i) +
       theme(plot.title = element_text(hjust = 0.5)) -> p
     plot_list[[i]] <- p
   }
   
   png(file = paste("PosEnrich", anno, name, "barplot.png", sep = "_"), 
         width = 8000, height = 6000, res=300)
   print(
     ggarrange(plotlist =plot_list, ncol = 3, nrow = 4)
     )
   dev.off()
   }
 }
   
set.seed(1220)
file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"
feature.ident <- c("Variable_features")
#' Other features : 
#' "Variable_a500UMI_features","Scaled_features","Scaled_a500UMI_features","All_features"
test.ident <- c("MAST", "roc")

#' this is used for patient analysis
ident.1 = "Fimbrial" 
ident.2 = "Otherside"

for (i in seq_along(feature.ident)) {
  for (t in seq_along(test.ident)) {
    for (p in c(5,6)) {
      patient <- paste("Patient", p, sep = "_")
      assay <- paste(patient, ident.1, ident.2, sep = "_")
      input.dir <- paste(deg.dir, 
                         file.name, 
                         cluster.ident, 
                         feature.ident[[i]], 
                         test.ident[[t]], 
                         assay,
                         sep = "/")
      deg.files <- list.files(path = input.dir, pattern = "*.filt.csv")
      
    for (x in seq_along(anno.files)) {
      for (y in seq_along(deg.files)) {
          if (y > length(deg.files)) break 
          try( plot_compareCluster(file      		= deg.files[[y]], 
          			                   anno      		= names[x], 
                         	        fun       		= "enricher",
                         	        dir       		= input.dir, 
                         	        TERM2GENE 		= anno.files[[x]],
                      		        pvalueCutoff 	= 0.01,
                      		        qvalueCutoff 	= 0.05)
                      )
          y <- y + 1
      }
    }
  }
  }
}




plot_enrichBar(file         = deg.files[[y]], 
               anno         = names[x],
               assay        = "enricher", 
               dir          = input.dir,
               TERM2GENE    = anno.files[[x]],
               pvalueCutoff = 0.01,
               qvalueCutoff = 0.05)
