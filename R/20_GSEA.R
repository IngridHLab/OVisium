#' Functional annotation using all positive regulated genes
#' Use all differential expressed genes for functional annotation of 
#' each cluster
#' as well as GSEA Analysis on lastest version of hallmark geneset
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

#' import all annotation files
setwd(deg.dir)
features.id <- readRDS("features_EntrezID.rds")
anno.files <- readRDS("FunctionAnnotation_entrez_files.rds")
names <- names(anno.files)

#' GSEA function on positive regulated genes 
plot_GSEA <- function(file, anno, assay, dir, ...) {
  
  setwd(dir)
  name <- str_replace_all(file, ".markers.filt.csv", "")
  
  #' Normally the universe should be the list of genes that were actually assayed 
  #' in your transcriptome analysis. 
  #' Otherwise all the GO genes (all genes with GO terms)
  universe <- features.id$ENTREZID
  markers <- readr::read_csv(file, col_names = TRUE)

  #' Convert gene name to entrezid
  common <- inner_join(markers[,c("cluster", "gene", "avg_log2FC")], 
                       features.id, by = c("gene" = "SYMBOL"))
  
  #' Check genes that don't have EntrezID
  warning(paste(setdiff(markers$gene, common$gene), 
                "skipped; ", sep = " "))
  
  #' reorder the clusters
  #' change the clusters to 0_FTE and 10_Stroma for the paired analysis
  common$cluster <- factor(common$cluster)

  result <- list()
  out.dir <- paste(dir, "GSEA", sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  for (c in seq_along(levels(common$cluster))) {
  #' Prepare a gene list for GESA 
    c.name <- levels(common$cluster)[c]
    common.cluster <- subset(common, cluster == c.name) 
    gene <- common.cluster$avg_log2FC
    names(gene) <- common.cluster$ENTREZID
    gene <- sort(gene, decreasing = T)
    
  #' Perform GSEA on cluster
     if (assay == "gseGO") {
         gse <- gseGO(geneList = gene,
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",
                      seed = TRUE,
                      ...)
    } else if (assay == "GSEA") {
        gse <- GSEA(geneList = gene,
                    seed = TRUE,
                    ...)
    }
    
    if (dim(gse)[1] == 0 || length(gse) == 0) {
      warning(paste0(c.name, " has no gene set enrichment!"))
      next     
    } else {
      result[[c.name]] <- gse
      }
    }
  
  if (length(result) == 0) {
    warning(paste0(name, " has no gene set enrichment result for any of the 
                    clusters!"))
    next
  } else {
    saveRDS(result, file = paste("GSEA", anno, name, "result.rds", sep = "_"))
    
    plot_list_1 = list()
    for (i in names(result)) {
      n <- nrow(result[[i]]) 
      if (n == 1) {
         p <- gseaplot2(result[[i]], 
                    geneSetID = n, 
                    pvalue_table = F,
                    title = i, 
                    subplots = 1:3) 
         p <- ggarrange(plotlist = p, ncol = 1, align = "v", 
                        heights = c(4, 1, 2), common.legend = TRUE)
        } else if (n > 5) {
            p <- gseaplot2(result[[i]], 
                           geneSetID = 1:5,  
                           pvalue_table = F, 
                           title = i,
                           subplots = 1:3)
            p <- ggarrange(plotlist = p, ncol = 1, align = "v", 
                           heights = c(4, 1, 2))
          } else {
              p <- gseaplot2(result[[i]], 
                             geneSetID = 1:n,  
                             pvalue_table = F,
                             title = i,
                             subplots = 1:3)
              p <- ggarrange(plotlist = p, ncol = 1, align = "v", 
                             heights = c(4, 1, 2))
          } 
      plot_list_1[[i]] <- p
  }
    png(file = paste("GSEA", anno, name, "runplot.png", sep = "_"), 
      width = 5000, height = 6000, res = 300)
    print(
      ggarrange(plotlist=plot_list_1, ncol = 3, nrow = 4)
    )
    dev.off()
  
    plot_list_2 = list()
    for (i in names(result)) {
      n <- nrow(result[[i]])
      result[[i]] <- setReadable(result[[i]], "org.Hs.eg.db", 
                               keyType = "ENTREZID")
      p <- cnetplot(result[[i]], 
                    categorySize       ="pvalue", 
                    foldChange         = gene, 
                    cex_label_category = 1, 
                    cex_label_gene     = 0.5, 
                    legend_n           = 2) + 
        ggtitle(i) + 
        theme(plot.title = element_text(size = 20))
      plot_list_2[[i]] <- p
    }
    
    png(file = paste("GSEA", anno, name, "cnet.png", sep = "_"), 
        width = 5000, height = 6000, res = 300)
    print(
      ggarrange(plotlist=plot_list_2, ncol = 3, nrow = 4)
    )
    dev.off()
  }
}

#' perform GSEA on all genesets across all clusters
#' pvalueCutoff = 0.05 as the standard
set.seed(1220)
file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"
feature.ident <- c("Variable_features")
#' Other features : 
#' "Variable_a500UMI_features","Scaled_features","Scaled_a500UMI_features","All_features"
test.ident <- c("MAST", "roc")

#' this is only used for patient analysis
ident.1 = "Otherside" 
ident.2 = "Proximal"


for (i in seq_along(feature.ident)) {
  for (t in seq_along(test.ident)) {
    for (p in c(6)) {
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
          try( plot_GSEA(file      = deg.files[[y]], 
                         anno      = names[[x]], 
                         assay     = "GSEA",
                         dir       = input.dir,
                         TERM2GENE = anno.files[[x]]))
          y <- y + 1
      }
    }
    }
  }
  }