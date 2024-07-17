#' Visualization of differential gene expression on reference cell markers
library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

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
setwd(gsea.dir)
feature.list <- list()
hm <- 
  readRDS(paste(gsea.dir, "h.all.v2023.2_FTE_result.rds", sep = "/"))
cluster <- "0_FTE"
core <- hm[[cluster]]@result[["core_enrichment"]] %>% 
  str_split(., pattern = "/") 

for (i in 1:3) {
  gene <- core[[i]] %>% data.frame() %>% dplyr::rename("ENTREZID"=1) %>% 
    inner_join(., features.id) %>% anti_join(., sc.markers, by="SYMBOL")
  feature.list[[paste(cluster, hm[[cluster]]@result[["ID"]][i], sep = "_")]] <- as.character(gene$SYMBOL)
}
  
  #' refined based on dotplot
  feature.list[["TNFA_SIGNALING (FTE)"]] <- c("TNFAIP2","TNFAIP3","CCL2","CLCF1",
                                              "CD83","IRF1","ZC3H12A","SAT1",
                                              "GADD45A")


#' hallmark features in cluster 10 
hm <- 
  readRDS(paste(gsea.dir, "h.all.v2023.2_Stroma_result.rds", sep = "/"))
cluster <- "10_Stroma"
core <- hm[[cluster]]@result[["core_enrichment"]] %>% 
  str_split(., pattern = "/") 
for (i in c(1:12)) {
  gene <- core[[i]] %>% data.frame() %>% dplyr::rename("ENTREZID"=1) %>% 
    inner_join(., features.id) %>% anti_join(., sc.markers, by="SYMBOL")
  feature.list[[paste(cluster, hm[[cluster]]@result[["ID"]][i], sep = "_")]] <- as.character(gene$SYMBOL)
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
features <- feature.list[[1]]
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

#' save feature list and export as table for loupe browser since these two clusters
#' share TNFa signaling
saveRDS(feature.list, "h.all.0_FTE&10_Stroma_genelist.rds")
feature.list.table <-melt(feature.list) %>% 
  dplyr::rename("Name"=1, "List" = 2) %>% 
  dplyr::relocate(Name, .after = List)
write.table(feature.list.table, 
            file = "GSEA_hallmark_0FTE_10Stroma.csv",
            sep =",", quote = F, row.names = F, col.names = T)

#' Plot individual gene set that enriched in our data
for (i in 1:length(feature.list)) {
  features <- feature.list[[i]]
  gene.set.id <- names(feature.list[i])
  
  #' split by Tissue_origin
  png(file = paste0("dotplot_", gene.set.id, "_by_cluster_origin.png"), 
      width = 5000, height = 3000, res = 400)
  p <- DotPlot(data.sub.filt, assay = "SCT", 
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
  print(p)
  dev.off()
  
  #' No split and just cluster
  png(file = paste0("dotplot_", gene.set.id, "_by_cluster.png"), 
      width = 5000, height = 2800, res = 400)
  p<-DotPlot(data.sub.filt, assay = "SCT", 
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
  print(p)
  dev.off()
}

#' hallmark features for all FTE and stroma clusters when they compared on their own
all.hm.gene.list <- list()
hm.FTE <- 
  readRDS(paste(gsea.dir, "h.all.v2023.2_FTE_result.rds", sep = "/"))
hm.stroma <- 
  readRDS(paste(gsea.dir, "h.all.v2023.2_Stroma_result.rds", sep = "/"))
all.hm <<- c(hm.FTE, hm.stroma)
for (i in 1:length(all.hm)) {
  cluster <- names(all.hm[i])
  core <- all.hm[[cluster]]@result[["core_enrichment"]] %>% str_split(., pattern = "/") 
  for (i in 1:length(core)) {
    gene <- core[[i]] %>% data.frame() %>% dplyr::rename("ENTREZID"=1) %>% inner_join(., features.id) 
    all.hm.gene.list[[paste(cluster, all.hm[[cluster]]@result[["ID"]][i], sep = "_")]] <- as.character(gene$SYMBOL)
  }
}

#' dotplot for all enriched gene set from all clusters 
for (i in 1:length(all.hm.gene.list)) {
  features <- all.hm.gene.list[[i]]
  gene.set.id <- names(all.hm.gene.list[i])
  
  #' split by Tissue_origin
  png(file = paste0("dotplot_", gene.set.id, "_by_cluster_origin.png"), 
      width = 5000, height = 3000, res = 400)
  p <- DotPlot(data.sub.filt, assay = "SCT", 
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
  print(p)
  dev.off()
  
  #' No split and just cluster
  png(file = paste0("dotplot_", gene.set.id, "_by_cluster.png"), 
      width = 5000, height = 2800, res = 400)
  p<-DotPlot(data.sub.filt, assay = "SCT", 
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
  print(p)
  dev.off()
}


#' immune sCell markers features 10 vs 11 
feature.list <- list()
scmark <- 
  readRDS(paste(gsea.dir, "sCell.markers.2.0_all_result.rds", sep = "/"))
cluster <- "11_Stroma"
core <- scmark[[cluster]]@result[["core_enrichment"]] %>% 
  str_split(., pattern = "/") 
for (i in 1:length(core)) {
  gene <- core[[i]] %>% data.frame() %>% dplyr::rename("ENTREZID"=1) %>% 
    inner_join(., features.id) %>% anti_join(., sc.markers, by="SYMBOL")
  feature.list[[scmark[[cluster]]@result[["ID"]][i]]] <- as.character(gene$SYMBOL)
}

#' Exclude microglial cell when analyse the fimbrial result
feature.list <- feature.list[c("Macrophage", "Monocyte", "Dendritic cell",
                               "T cell", "Natural killer cell", "B cell")] 
#' combine T and NK cells
feature.list[["T/NK cell"]]  <- c(feature.list[["T cell"]], feature.list[["Natural killer cell"]]) 

feature.list <- feature.list[c(1,2,3,7,6)]
un <- feature.list %>% unlist() 
unfeature.list <- Map(`[`, feature.list, relist(!duplicated(un), 
                                                skeleton = feature.list))

#' Dotplot of new markers of merged data by cellmarker
out.dir <- paste(gsea.dir, "new_markers", sep = "/")
dir.create(out.dir, recursive = T)
setwd(out.dir)
features <- unfeature.list
png(file = paste0("dotplot_scMark_immune_by_cluster_origin.png"), 
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

png(file = paste0("dotplot_scMark_immune_by_cluster.png"), 
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



#' FTE markers between 0 and 1 clusters 
out.dir <- paste(deg.dir, file.name, cluster.ident, feature.ident, test.ident, sep = "/")
dir.create(out.dir, recursive = T)
setwd(out.dir)
FTE.markers <- 
  readr::read_csv("Clusters.0vs1.markers.filt.top100.csv")
sc.markers <- sc.markers %>% dplyr::rename("gene"=1)
new.FTE.markers <- anti_join(FTE.markers, sc.markers)
new.FTE.markers$pct.diff<-abs(new.FTE.markers$pct.1-new.FTE.markers$pct.2)
new.FTE.markers<-new.FTE.markers %>% arrange(desc(pct.diff))

cilia <- new.FTE.markers[which(str_detect(new.FTE.markers$description, 
                                          pattern = "cilia*") == T),]$gene
memb <- new.FTE.markers[which(str_detect(new.FTE.markers$description, 
                                        pattern = "membrane") == T),]$gene
bind <- new.FTE.markers[which(str_detect(new.FTE.markers$description, 
                                        pattern = "binding protein") == T),]$gene
grow <- new.FTE.markers[which(str_detect(new.FTE.markers$description, 
                                         pattern = "growth|cyclin|proliferation") == T),]$gene
trans <- new.FTE.markers[which(str_detect(new.FTE.markers$description, 
                                          pattern = "transcript|translat") == T),]$gene
FTE.features = list("Cilia" = cilia, 
                "Transmemb" = memb, 
                "BindingPro" = bind,
                "CellCycleGrowth" = grow, 
                "GeneExpression" = trans)

un <- FTE.features %>% unlist() 
unFTE.features <- Map(`[`, FTE.features, relist(!duplicated(un), 
                                                skeleton = FTE.features))

#' Dotplot of reference markers of merged data
features <- unFTE.features
png(file = paste("dotplot_0vs1_new_by_cluster_origin.png", sep = "."), 
    width = 5000, height = 3500, res = 400)
DotPlot(data, assay = "SCT", 
        features = features, 
        group.by = "Visium_clusters", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 10,
        split.by = "Tissue_origin") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(angle = 0, size = 11)) +
  guides(size=guide_legend(title="Fraction of spots", 
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("") + 
  theme(legend.position = "bottom")
dev.off()

png(file = paste("dotplot_0vs1_new_by_cluster.png", sep = "."), 
    width = 5000, height = 2200, res = 400)
DotPlot(data, assay = "SCT", 
        features = features, 
        group.by = "Visium_clusters", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 10) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(angle = 0, size = 11)) +
  guides(size=guide_legend(title="Fraction of spots", 
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("") + 
  theme(legend.position = "bottom")
dev.off()

#' twoside barplot
#' https://www.genepattern.org/modules/docs/GSEA/17#gsc.tab=0
#' qvalue cutoff 0.25 for large genesets

result.hm <- mutate(hallmark[[cluster]]@result, ordering = abs(NES)) %>%
  arrange(desc(ordering))

## Warning: The `add` argument of `group_by()` is deprecated as of dplyr 1.0.0.
## Please use the `.add` argument instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_warnings()` to see where this warning was generated.
setwd(out.dir)
png(file = paste0("tsbarplot_gsea_hallmark_", cluster, ".png"), 
    width = 2500, height = 2000, res = 400)
ggplot(result.hm,
       aes(NES, fct_reorder(Description, NES), 
           fill = qvalue < 1e-04), showCategory=(20)) + 
  geom_bar(stat = "identity", color="black", width = 0.25) +  # to adjust width = 0.3
  scale_fill_manual(values = c("#FC2D00","#008EFC"),
                    labels=c("TRUE"= "FDR q-value < 1e-04",
                             "FALSE"= "FDR q-value > 1e-04")) +
  labs(x="Normalized Enrichment Score (NES)",
       y=NULL,
       title = "Cluster 0 vs 1",
       subtitle = "",
       caption  = "", 
       fill = "") + 
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()

#' single cell markers
result.sc <- mutate(scmark[[cluster]]@result, ordering = abs(NES)) %>%
  arrange(desc(ordering))

n <- 10
result_bar <- group_by(result, sign(NES)) %>%
  slice(1:n)

png(file = paste0("tsbarplot_gsea_scell_", cluster, ".png"), 
    width = 2500, height = 2000, res = 400)
ggplot(result.sc,
       aes(NES, fct_reorder(Description, NES), 
           fill = qvalue < 1e-04), showCategory=(20)) + 
  geom_bar(stat = "identity", color="black", width = 0.9) +  # to adjust width = 0.3
  scale_fill_manual(values = c("#FC2D00","#008EFC"),
                    labels=c("TRUE"= "FDR q-value < 1e-04",
                             "FALSE"= "FDR q-value > 1e-04")) +
  labs(x="Normalized Enrichment Score (NES)",
       y=NULL,
       title = "Cluster 10 vs 11",
       subtitle = "",
       caption  = "", 
       fill = "") + 
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()
