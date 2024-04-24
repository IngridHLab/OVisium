library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))

my_merge <- function(df1, df2){                                
  merge(df1, df2, by = "gene")
}
spat.dir <- paste(work.dir, "Spatgenes", sep = "/")
file.list <- list.files(path=paste(spat.dir, "spat_corr_list"), pattern = ".csv")
spatgene.list <- list()
for (i in seq_along(file.list)) {
  spat <- readr::read_csv(file = paste(spat.dir, "spat_corr_list", file.list[i], sep = "/"), col_names = T)
  sample.id <- gsub(".csv", "", file.list[i])
  spatgene.list[[sample.id]] <- spat
}
spatgenes.common <- Reduce(my_merge, spatgene.list) 
colnames(spatgenes.common) <- c("gene", 1:18)
spatgenes.common$average <-rowMeans(spatgenes.common[2:19])
spatgenes.common <- spatgenes.common[order(-spatgenes.common$average),] 
write.table(spatgenes.common, file=paste0(work.dir, "/Spatgenes/", "common_spatgenes.csv"), sep =",", quote = F, row.names = F, col.names = T)

#' Visualize the top genes in dotplot, vlnplot and feature plot
spatgenes.common.filt <- spatgenes.common %>% 
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene))
FTE.features <- spatgenes.common.filt$gene[c(12,19,20,25,28,30,36,39,44,45,48)]
Stroma.features <- base::setdiff(spatgenes.common.filt$gene, FTE.features)[1:39]
features <- list("FTE" = FTE.features, "Stroma" = Stroma.features)

setwd(spat.dir)
png(file = paste("dotplot_spatgene50_by_cluster.png", sep = "."), 
    width = 5000, height = 1800, res = 350)
DotPlot(data, assay = "SCT", 
        features = features, 
        group.by = "Visium_clusters", 
        cluster.idents = F, 
        cols="RdYlBu", 
        dot.scale = 7) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(angle = 0, size = 14)) +
  guides(size=guide_legend(title="Fraction of spots", 
                           override.aes=list(shape=21, colour="grey", 
                                             fill="grey"))) +
  xlab("") + 
  ylab("") 
dev.off()

png(file = paste("vlnplot_spatgene50_by_cluster.png", sep = "."), 
    width = 8000, height = 7000, res = 400)
VlnPlot(data, assay = "SCT", 
        group.by = "Visium_clusters", 
        features = unlist(features),
        ncol = 8,
        pt.size = 0) & 
  xlab("") & 
  ylab("")
dev.off()

Idents(data) <- "Visium_clusters"
#' FeaturePlot of reference markers of merged data
png(file = paste("featureplot_spatgene50_uwot.png", sep = "."), 
    width = 8000, height = 7000, res = 300)
FeaturePlot(data, 
            features = unlist(features), 
            reduction = "uwot_harmony_SCT", 
            slot = "data", 
            ncol= 8, 
            order = F, 
            label = T, 
            label.size = 2, 
            label.color = "black", 
            repel = T,
            combine = T) & 
  NoAxes() &
  NoLegend() &
  labs(color = "Data")
dev.off()
