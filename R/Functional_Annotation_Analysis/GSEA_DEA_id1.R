#' Functional annotation using all positive regulated genes
#' Use all differential expressed genes for functional annotation of 
#' each cluster
#' as well as GSEA Analysis on lastest version of hallmark geneset
library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

#' import all annotation files
setwd(deg.dir)
features.id <- readRDS("features_EntrezID.rds")
anno.files <- readRDS("FunctionAnnotation_entrez_files.rds")
names <- names(anno.files)

#' GSEA function on positive regulated genes 
plot_GSEA_id1 <- function(file, anno, assay, dir, ...) {
  
  setwd(dir)
  name <- gsub(".findmarkers.*", "", file)
  
  #' Normally the universe should be the list of genes that were actually assayed 
  #' in your transcriptome analysis. 
  #' Otherwise all the GO genes (all genes with GO terms)
  universe <- features.id$ENTREZID
  markers.list <- readRDS(file)

  markers<-bind_rows(markers.list)
  #' Convert gene name to entrezid
  common <- inner_join(markers, features.id, by = c("gene" = "SYMBOL"), 
                       relationship = "many-to-many")
  
  #' Check genes that don't have EntrezID
  warning(paste(setdiff(markers$gene, common$gene), 
                "skipped; ", sep = " "))
  
  common$id1 <- factor(common$id1)
  common$id2 <- factor(common$id2)
                          
  out.dir <- paste(dir, "GSEA", anno, sep = "/")
  dir.create(out.dir, recursive = T)
  sub.dir.names <- c("rds", "runplot", "cnet", "tsbar")
  for (n in sub.dir.names) {
    sub.dir <- paste(out.dir, n, sep = "/")
    dir.create(sub.dir)
  }
  
  setwd(out.dir)
  result <- list()
  
  for (c in seq_along(levels(common$id1))) {
  #' Prepare a gene list for GESA 
    c.name <- levels(common$id1)[c]
    common.id1 <- subset(common, id1 == c.name) 
    gene <- common.id1$avg_log2FC # or avg_log2FC
    names(gene) <- common.id1$ENTREZID
    gene <- sort(gene, decreasing = T)
    
  #' Perform GSEA on cluster
    out.name <- paste(anno, name, sep = "_")
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
    
    saveRDS(result, file = paste0("rds/", out.name, "_result.rds"))
    
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
    png(file = paste0("runplot/", out.name, "_runplot.png"), 
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
    
    png(file = paste0("cnet/", out.name, "_cnet.png"), 
        width = 5000, height = 6000, res = 300)
    print(
      ggarrange(plotlist=plot_list_2, ncol = 3, nrow = 4)
    )
    dev.off()
    
    plot_list_3 = list()
    for (i in names(result)) {
      n <- nrow(result[[i]]) 
      result[[i]] <- mutate(result[[i]]@result, ordering = abs(NES)) %>%
        arrange(desc(ordering))
      
      if (n > 10) { 
        result[[i]] <- group_by(result[[i]], sign(NES)) %>%
        slice(1:10)
      } else { 
        result[[i]] <- group_by(result[[i]], sign(NES)) %>%
        slice(1:n)
      }
      
      p <- ggplot(result[[i]],
                  aes(NES, fct_reorder(Description, NES), 
                      fill = qvalue < 1e-04), showCategory=(20)) + 
        geom_bar(stat = "identity", color="black", width = 0.9) +  # to adjust width = 0.3
        scale_fill_manual(values = c("#FC2D00","#008EFC"),
                          labels=c("TRUE"= "FDR q-value < 1e-04",
                                   "FALSE"= "FDR q-value > 1e-04")) +
        labs(x="Normalized Enrichment Score (NES)",
             y=NULL,
             title = i,
             subtitle = "",
             caption  = "", 
             fill = "") + 
        theme_classic() +
        theme(legend.position = "bottom")
      plot_list_3[[i]] <- p
    }
    
    png(file = paste0("tsbar/", out.name, "_tsbar.png"), 
        width = 5000, height = 9000, res = 300)
    print(
      ggarrange(plotlist=plot_list_3, ncol = 2, nrow = 6)
    )
    dev.off()
    
  }
}

#' perform GSEA on all genesets across all clusters
#' pvalueCutoff = 0.05 as the standard
set.seed(1220)
file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"
feature.ident <- c("Variable_features_filt")
#' Other features : 
#' "Variable_a500UMI_features","Scaled_features","Scaled_a500UMI_features","All_features"
#' Different list based on the filtering
#' Test methods to identify DE genes
test.ident <- c("MAST")
markers.list.dir <- "FindMarkers/Combined_filtered_markers_v1/rds"

for (i in seq_along(feature.ident)) {
  for (t in seq_along(test.ident)) {
      input.dir <- paste(deg.dir, 
                         file.name, 
                         cluster.ident, 
                         feature.ident[[i]], 
                         test.ident[[t]],
                         markers.list.dir,
                         sep = "/")
      
      deg.files <- list.files(path = input.dir, pattern = "*rds")
    
    for (x in seq_along(anno.files)) {
      for (y in seq_along(deg.files)) {
          if (y > length(deg.files)) break 
           try(plot_GSEA_id1(file      = deg.files[[y]], 
                         anno      = names[[x]], 
                         assay     = "GSEA",
                         dir       = input.dir,
                         TERM2GENE = anno.files[[x]])
               )
          y <- y + 1
      }
    }
    }
  }
