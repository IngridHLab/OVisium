#' Prepare standard gene list for compareCluster 
#' perform geneset enrichment analysis
library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

#' import all annotation files
setwd(deg.dir)
features.id <- readRDS("features_EntrezID.rds")
anno.files <- readRDS("FunctionAnnotation_entrez_files.rds")
names <- names(anno.files)

#' Enrichment analysis between clusters
plot_compareCluster_patient <- function(file, anno, dir, ...) {

setwd(dir)
name <- gsub(".list.*", "", file)

universe <- features.id$ENTREZID
markers.list <- readRDS(file)
markers<-bind_rows(markers.list)

#' Convert gene name to entrezid
common <- inner_join(markers, features.id, by = c("gene" = "SYMBOL"), 
                     relationship = "many-to-many")

#' Check genes that don't have EntrezID
warning(paste(setdiff(markers$gene, common$gene), 
              "skipped; ", sep = " "))

common$patient <- factor(common$patient, levels = c("Patient_1", "Patient_3",
                                                    "Patient_6", "Patient_10"))


common$level[common$avg_log2FC >= log2(1.5)] <- "Up" # FC > 1.5
common$level[common$avg_log2FC <= log2(0.75)] <- "Down" # FC < 0.75
common <- na.omit(common)

g1 <- c("Patient_1", "Patient_3",
        "Patient_6", "Patient_10")
g2 <- c("Up", "Down", "No Change")


out.dir <- paste(dir, "Enrichment", anno, sep = "/")
dir.create(out.dir, recursive = T)
sub.dir.names <- c("rds", "dotplot")
for (n in sub.dir.names) {
  sub.dir <- paste(out.dir, n, sep = "/")
  dir.create(sub.dir)
}

setwd(out.dir)
result <- compareCluster(ENTREZID~patient+level, 
                         data         = common, 
                         universe     = universe,
                         ...)
if(length(result) == 0) {
  warning(paste0(name, " has no result from compareCluster!"))
  next
} else {
  result@compareClusterResult$patient <- 
    factor(result@compareClusterResult$patient, levels = g1)
  result@compareClusterResult$level <- 
    factor(result@compareClusterResult$level, levels = g2 )
  saveRDS(result, 
          file = paste("rds/CompareCluster", anno, name, "result.rds", sep = "_"))
  
  if (length(unique(result@compareClusterResult$Description))>100){
    png(file = paste("dotplot/CompareCluster", anno, name, "dotplot.png", sep = "_"), 
        width = 2400, height = 3500, res=300)
    print(
      dotplot(result, x="patient", showCategory = 10, includeAll = F) + 
        facet_grid(~level) + 
        theme(axis.text.x = element_text(hjust = 1, angle = 60), 
              axis.text.y = element_text(size = 7)) +
        ggtitle(paste("Enrichment of", anno, "for", name, sep = " "))
    )
    dev.off()
  } else {
    png(file = paste("dotplot/CompareCluster", anno, name, "dotplot.png", sep = "_"), 
        width = 2400, height = 3000, res=300)
    print(
      dotplot(result, x="patient", showCategory = 10, includeAll = F) + 
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
plot_enrichBar_patient <- function(file, anno, assay, dir, ...) {
  
  setwd(dir)
  name <- gsub(".list.*", "", file)
  
  universe <- features.id$ENTREZID
  markers.list <- readRDS(file)
  markers<-bind_rows(markers.list)
  
  #' Convert gene name to entrezid
  common <- inner_join(markers, features.id, by = c("gene" = "SYMBOL"), 
                       relationship = "many-to-many")
  
  #' Check genes that don't have EntrezID
  warning(paste(setdiff(markers$gene, common$gene), 
                "skipped; ", sep = " "))
  
  common$patient <- factor(common$patient, levels = c("Patient_1", "Patient_3",
                                                      "Patient_6", "Patient_10"))

  out.dir <- paste(dir, "Enrichment", anno, sep = "/")
  dir.create(out.dir, recursive = T)
  sub.dir.names <- c("rds", "barplot")
  for (n in sub.dir.names) {
    sub.dir <- paste(out.dir, n, sep = "/")
    dir.create(sub.dir)
  }
  
  setwd(out.dir)
  result <- list()
  
  for (c in 1:length(levels(common$patient))) {
    
    c.name <- levels(common$patient)[c]
    common.patient <- subset(common, patient == c.name)
    gene <- common.patient$ENTREZID
    
    #' Enrichment analysis on individual cluster
    if (assay == "enrichGO") {
      eer <- enrichGO(gene          = gene,
                      universe      = universe, 
                      OrgDb         = org.Hs.eg.db,
                      ...)
    } else {
      eer <- enricher(gene          = gene,
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
    saveRDS(result, file = paste("rds/PosEnrich", anno, name, "result.rds", sep = "_"))
    
    plot_list <- list()
    for (i in names(result)) {
      barplot(result[[i]], order = TRUE, showCategory = 10)  +
        theme(axis.text = element_text(size = 4)) +
        ggtitle(i) +
        theme(plot.title = element_text(hjust = 0.5)) -> p
      plot_list[[i]] <- p
    }
    
    png(file = paste("barplot/PosEnrich", anno, name, "barplot.png", sep = "_"), 
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
 feature.ident <- c("Variable_features_filt")
 #' Other features : 
 #' "Variable_a500UMI_features","Scaled_features","Scaled_a500UMI_features","All_features"
 #' Different list based on the filtering
 #' Test methods to identify DE genes
 test.ident <- c("MAST")
 markers.list.dir <- "FindMarkers/2024-05-16_patient_pair"
 
 for (i in seq_along(feature.ident)) {
   for (t in seq_along(test.ident)) {
     input.dir <- paste(deg.dir, 
                        file.name, 
                        cluster.ident, 
                        feature.ident[[i]], 
                        test.ident[[t]],
                        markers.list.dir,
                        sep = "/")
     
     deg.files <- list.files(path = input.dir, pattern = "*DEA.list.rds")
      
    for (x in seq_along(anno.files)) {
      for (y in seq_along(deg.files)) {
          if (y > length(deg.files)) break 
          try( plot_compareCluster_patient(file      		= deg.files[[y]], 
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


 for (i in seq_along(feature.ident)) {
   for (t in seq_along(test.ident)) {
     input.dir <- paste(deg.dir, 
                        file.name, 
                        cluster.ident, 
                        feature.ident[[i]], 
                        test.ident[[t]],
                        markers.list.dir,
                        sep = "/")
     
     deg.files <- list.files(path = input.dir, pattern = "*DEA.list.rds")
     
     for (x in seq_along(anno.files)) {
       for (y in seq_along(deg.files)) {
         if (y > length(deg.files)) break 
         try( plot_enrichBar_patient(file         = deg.files[[y]], 
                             anno         = names[x],
                             assay        = "enricher", 
                             dir          = input.dir,
                             TERM2GENE    = anno.files[[x]],
                             pvalueCutoff = 0.01,
                             qvalueCutoff = 0.05)
         )
         y <- y + 1
       }
     }
   }
 }





