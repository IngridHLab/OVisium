#' some scripts for making different heatmap

#' Original Seurat Heatmap on all clusters
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
clust.col <- 
  c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#F0027F","#BF5B17",
             "#666666","#1B9E77","#D95F02","#386CB0","black")
             
names(clust.col) <-c("0_FTE","1_FTE", "3_Mix", "2_Stroma",
                     "5_Stroma", "6_Stroma", "7_Stroma",
                     "8_Stroma", "9_Stroma", "10_Stroma",
                     "11_Stroma")

p<-DoHeatmap(pat.data, assay = "SCT", 
             features = HRD_gene_signature,
             group.by = "origin", #'or "Tissue_origin"
             group.colors = clust.col, 
             size = 8, 
             angle = 70, 
             disp.max = 3) +
  theme(text = element_text(size = 18)) +
  theme(axis.text = element_text(size = 10)) + 
  guides(color="none") +
  labs(fill = "SCT:Residual") +
  scale_fill_gradientn(colors = rev(mapal))
p_b <- ggplot_build(p)
p_b$layout$panel_scales_y[[1]]$range$range
png(file = paste(patient, "heatmap_HRD_signature_144_origin.png", sep = "_"),
    width = 8000 , height = 9000, res = 400) 
print( 
  p + scale_y_discrete(limits =  rev(p_b$layout$panel_scales_y[[1]]$range$range))
)
dev.off()


#' multibar heatmap
cols.use <- list(Visium_clusters=c("#7FC97F","#BEAED4","#FDC086","#FFFF99",
                                            "#F0027F","#BF5B17","#666666",
                                            "#1B9E77","#D95F02","#386CB0",
                                            "black", "white"))
                                            names(cols.use[['Visium_clusters']]) <-c("0_FTE","1_FTE", "3_Mix", "2_Stroma",
                                                                                     "5_Stroma", "6_Stroma", "7_Stroma",
                                                                                     "8_Stroma", "9_Stroma", "10_Stroma",
                                                                                     "11_Stroma", "empty")
                                            
                                            
                                            p <- DoMultiBarHeatmap(pat.data,
                                                                   features = HRD_gene_signature, 
                                                                   group.by = "origin",
                                                                   cols.use = cols.use,
                                                                   additional.group.by = c( "Visium_clusters", "Phase"),
                                                                   additional.group.sort.by = c( "Visium_clusters"),
                                                                   disp.max = 3) +
                                              theme(text = element_text(size = 18)) +
                                              theme(axis.text = element_text(size = 10)) + 
                                              labs(fill = "SCT:Residual") +
                                              scale_fill_gradientn(colors = rev(mapal))
                                            p_b <- ggplot_build(p)
                                            p_b$layout$panel_scales_y[[1]]$range$range
                                            
                                            png(file = paste(patient, "mbHeatmap_HRD_signature_144.png", sep = "_"), 
                                                width = 8000 , height = 9000, res = 400) 
                                            print( 
                                              p + scale_y_discrete(limits = rev(p_b$layout$panel_scales_y[[1]]$range$range)) 
                                            )
                                            dev.off()
                                            
                                            }

#' original pheatmap (classic)
pat.m <- as.matrix(GetAssayData(pat.data, slot = "scale.data"))
pat.m.HRD <- subset(pat.m, rownames(pat.m) %in% HRD_gene_signature)
quantile(pat.m.HRD, c(0.1, 0.95))

#' prepare annotations for pheatmap
annotation_col = data.frame(Cluster = pat.data$Visium_clusters,
                            Origin = pat.data$Tissue_origin,
                            variant = factor(pat.data$variant),
                            age = factor(pat.data$age))

annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
)
rownames(annotation_row) = paste("Gene", 1:20, sep = "")


ann_colors = list(
  Cluster = clust.col,
  Origin = c(Fimbrial = "#1B9E77", Proximal = "#D95F02"),
  variant = c(c.3048_3052dup = "#7570B3", c.68_69del = "#E7298A", c.1082_1092del = "#66A61E"),
  age = c("35+"="#8B7500", "40+"="#FFD700", "45+"="#FFBBFF", "50+"="#9400D3")
)
png(file = paste(patient, "pHeatmap_HRD_signature_144.png", sep = "_"), 
    width = 8000 , height = 9000, res = 400) 

ComplexHeatmap::pheatmap(subset.pat.m, annotation_col = annotation_col,
                         annotation_colors = ann_colors, show_colnames = FALSE, treeheight_col = 0)
dev.off()
