#' Visualization of differential gene expression on reference cell markers
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
my_merge <- function(df1, df2){                                
  merge(df1, df2, by = "gene")
}

theme_splittext = function (...) 
{
  function(label, ...) {
    splitlab = paste(strwrap(label), collapse="\n")
    textGrob(splitlab, 0.5, 0.5,  ...) 
  }
}

file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"
feature.ident <- "Variable_features_filt"
test.ident <- "MAST"
markers.list.dir <- "FindMarkers/Combined_filtered_markers_v1/rds"

gsea.dir <- 
  paste(deg.dir, 
        file.name, 
        cluster.ident, 
        feature.ident, 
        test.ident,
        markers.list.dir,
        "GSEA", sep = "/")

#' Reference cell markers from single cell data
ref.markers<-read.csv(paste(deg.dir, "single_cell_markers.csv", sep = "/"))
sc.markers <- ref.markers$Hammound_2022 %>% data.frame() %>% dplyr::rename("SYMBOL"=1)
load(paste(rds.dir, "Variable_features_filt_SCT_log2counts+1_harmony.RData",sep = "/"))
Idents(data.sub.filt) <- factor(Idents(data.sub.filt),
                       levels = c("11_Stroma", "10_Stroma", "9_Stroma", 
                                  "8_Stroma", "7_Stroma", "6_Stroma", 
                                  "5_Stroma","2_Stroma","3_Mix","1_FTE", 
                                  "0_FTE"))
data.sub.filt[["Visium_clusters_rev"]] <- Idents(data.sub.filt)
features.id <- readRDS(paste(deg.dir, "features_EntrezID.rds", sep = "/"))


#' go features 0 vs 1
feature.list <- list()
hm <- 
  readRDS(paste(gsea.dir, "h.all.v2023.2_FTE_result.rds", sep = "/"))
cluster <- "0_FTE"
core <- hm[[cluster]]@result[["core_enrichment"]] %>% 
  str_split(., pattern = "/") 

  gene <- core[[2]] %>% data.frame() %>% dplyr::rename("ENTREZID"=1) %>% 
    inner_join(., features.id) %>% anti_join(., sc.markers, by="SYMBOL")
  feature.list[["TNFA_SIGNALING (FTE)"]] <- as.character(gene$SYMBOL)
  
  #' refined based on dotplot
  feature.list[["TNFA_SIGNALING (FTE)"]] <- c("TNFAIP2","TNFAIP3","CCL2","CLCF1",
                                              "CD83","IRF1","ZC3H12A","SAT1",
                                              "GADD45A")



#' hallmark features 10 vs 11
hm <- 
  readRDS(paste(gsea.dir, "h.all.v2023.2_Stroma_result.rds", sep = "/"))
cluster <- "10_Stroma"
core <- hm[[cluster]]@result[["core_enrichment"]] %>% 
  str_split(., pattern = "/") 
for (i in c(1,2,4,5,6)) {
  gene <- core[[i]] %>% data.frame() %>% dplyr::rename("ENTREZID"=1) %>% 
    inner_join(., features.id) %>% anti_join(., sc.markers, by="SYMBOL")
  feature.list[[hm[[cluster]]@result[["ID"]][i]]] <- as.character(gene$SYMBOL)
}

un <- feature.list %>% unlist()
unfeature.list <- Map(`[`, feature.list, relist(!duplicated(un), 
                                                skeleton = feature.list)) 
unfeature.list <- unfeature.list[c(1,2,3,4,5,6)]

names(unfeature.list)[1] <- "TNFa (FTE)"
names(unfeature.list)[2] <- "TNFa (Stroma)"
names(unfeature.list)[3] <- "HYPO"
names(unfeature.list)[4] <- "P53"
names(unfeature.list)[5] <- "APOPT"
names(unfeature.list)[6] <- "EMT"



#' Dotplot of hallmark features from both FTE and Stroma
out.dir <- paste(gsea.dir, "new_markers", sep = "/")
dir.create(out.dir, recursive = T)
setwd(out.dir)
features <- unfeature.list
png(file = paste0("dotplot_new_markers_by_cluster_origin.png"), 
    width = 7500, height = 3000, res = 400)
DotPlot(data.sub.filt, assay = "SCT", 
        features = features, 
        group.by = "Visium_clusters_rev", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 7,
        split.by = "Tissue_origin") +
  theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8),
        strip.text.x = element_text(angle = 0, size = 11)) +
  guides(size=guide_legend(title="Fraction of spots", 
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("") +
  theme(legend.position = "bottom")
dev.off()

png(file = paste0("dotplot_new_markers_by_cluster.png"), 
    width = 9000, height = 2800, res = 400)
DotPlot(data.sub.filt, assay = "SCT", 
        features = features, 
        group.by = "Visium_clusters_rev", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 7) + 
  theme(axis.text.x=element_text(angle = 45, hjust = 0.9, vjust = 0.8),
        strip.text.x = element_text(angle = 0, size = 15)) +
  guides(size=guide_legend(title="Fraction of spots", 
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("") +
  theme(legend.position = "bottom")
dev.off()
