#' Spatial visualization of clustering of all samples
#' https://ludvigla.github.io/STUtility_web_site/Maize_data.html#Clustering
#' Decide which data and resolution before ploting
library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
dir.create(paste(dims.dir, "Cluster", sep = "/"), recursive = T)

#' The optimal resolution based on morphology and preliminary DEG is 0.6 from 
#' the combine variable feature data
file.name <- list.files(rds.dir, pattern = "OVisium.*._merged_clust.rds")
data.type <- gsub("_clust.rds", "", file.name)
data <- readRDS(paste(rds.dir, file.name, sep = "/"))
  
#' plot each patient 
npcs <- length(data@reductions[["pca_SCT"]])
res <- 0.6
vars <- "sample"
Idents(data) <- vars
id <- levels(data)
data[["sample"]] <- Idents(data)
                                 

#' prepare color code
clust.col <- 
  c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#7570B3","#F0027F","#BF5B17",
             "#666666","#1B9E77","#D95F02","#386CB0","black")

cscale <- c("darkblue", "blue", "lightblue", "white", 
                      "lightgray", "mistyrose", "red", "darkred", "black")

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))  
                      
#' plot harmony integrated PCs embedding values spatially of each sample
for (i in seq_along(id)) {
  #' Plot PCA embedding value spatially
  out.dir <- paste(dims.dir, "PCA_SCT", data.type, sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  png(file = paste0("pca_spatial_", id[[i]], ".png"), 
      width = 5000, height = 5000, res = 400) 
  print(
    ST.DimPlot(data, 
               dims = 1:npcs,
               indices = i,
               grid.ncol = 5,
               reduction = "pca_SCT", 
               pt.size = 1, 
               center.zero = T, 
               cols = cscale, 
               show.sb = FALSE) &
      theme(legend.position = "None", plot.subtitle = element_blank())
  )
  dev.off() 
  
  #' plot Harmony embedding value spatially 
  out.dir <- paste(dims.dir, "HARMONY_SCT", data.type, sep = "/")
  dir.create(out.dir, recursive = T)
  setwd(out.dir)
  png(file = paste0("harmony_spatial_", id[[i]], ".png"), 
      width = 5000, height = 5000, res = 400) 
  print(
    ST.DimPlot(data, 
               dims = 1:npcs,
               indices = i,
               grid.ncol = 5,
               reduction = "harmony_SCT", 
               pt.size = 1, 
               center.zero = T, 
               cols = cscale, 
               show.sb = FALSE) &
      theme(legend.position = "None", plot.subtitle = element_blank())
  )
  dev.off() 
    
#' plot clusters spatially
    harmony.ident <- paste0("harmony_SCT_res_", res)
    Idents(data) <- harmony.ident
    out.dir <- paste(dims.dir, "Cluster", data.type, harmony.ident, sep = "/")
    dir.create(out.dir, recursive = T)
    setwd(out.dir)
    png(file = paste("clusters_spatial_", id[[i]], ".png", sep = ""), 
      width = 8000, height = 6000, res = 400) 
    p1 <- ST.FeaturePlot(data, 
                       features = harmony.ident, 
                       cols = clust.col,
                       pt.size = 1.2, 
                       split.labels = T, 
                       indices = i,
                       show.sb = FALSE, 
                       ncol = 3) & 
      theme(legend.position = "None", plot.title = element_blank()) + 
      theme(text = element_text(size = 20))
  
    p2 <- ST.FeaturePlot(data, 
                       features = harmony.ident,
                       cols = clust.col,
                       pt.size = 2, 
                       indices = i) + 
      labs(fill="Clusters") + 
      theme(text = element_text(size = 20), 
            plot.title = element_blank(), 
            legend.position = "bottom") +
      guides(fill = guide_legend(override.aes = list(size=10)))
    
    p3 <- FeatureOverlay(data, 
                         features = harmony.ident, 
                         sampleids = i,
                         cols = clust.col,
                         pt.size = 0, 
                         pt.alpha = 0, 
                         pt.border = F,
                         show.sb = FALSE) &
      theme(legend.position = "None", plot.title = element_blank(), 
            plot.subtitle = element_blank()) 
    
    p4 <- plot_grid(p2, p3, ncol = 1, rel_heights = c(1.2,1))
    print(    
      plot_grid(p1, p4, ncol = 2, rel_widths = c(1.8,1))
    )
    dev.off() 

    }
  
out.dir <- paste(dims.dir, "Cluster", data.type, harmony.ident, sep = "/")
setwd(out.dir)

#' plot all clusters spatially across all samples 
png(file = "clusters_spatial_all.png", width = 8500, height = 7000, res = 500) 
print(
  ST.FeaturePlot(data, 
                 features = harmony.ident, 
                 cols = clust.col, 
                 ncol = 5, 
                 show.sb = F) +
    labs(fill="Clusters") + 
    theme(text = element_text(size = 15)) +
    guides(fill = guide_legend(override.aes = list(size=10)))
)
dev.off()

#' plot all clusters spatially across all samples with H&E in the background
png(file = "clusters_spatialHE_all.png", width = 8500, height = 7000, res = 500) 
print(
  FeatureOverlay(data, 
                 features = harmony.ident, 
                 sampleids = 1:length(id),
                 pt.size = 1,  
                 pt.border = F,
                 pt.alpha = 0.5,
                 split.labels = T, 
                 show.sb = FALSE, 
                 ncol = 5) &
    theme(plot.title = element_blank()) &
    NoLegend()
  
)
dev.off() 

#' plot dimplot with harmony clustering annotation
png(file = "clusters_dimPlot_uwot.png", width = 1600, height = 1500, res = 300) 
print(
  DimPlot(data, 
          group.by = harmony.ident,
          reduction = "uwot_harmony_SCT", 
          cols = clust.col) + 
    xlab("uwot_harmony_SCT_1") + 
    ylab("uwot_harmony_SCT_2") 
)
dev.off()

png(file = "clusters_dimPlot_umap.png", width = 1600, height = 1500, res = 300) 
print(
  DimPlot(data, 
          group.by = harmony.ident,
          reduction = "umap_harmony_SCT", 
          cols = clust.col) + 
    xlab("umap_harmony_SCT_1") + 
    ylab("umap_harmony_SCT_2") 
)
dev.off()

#' One approach to visualize the result of dimensionality reduction is to use 
#' the first three dimensions and transform the values into RGB color space. 
#' This 3 dimensional umap can then be utilized for spatial visualization. 
out.dir <- paste(dims.dir, "HARMONY_SCT", data.type, sep = "/")
dir.create(out.dir, recursive = T)
setwd(out.dir)
png(file = "umap_3d_spatial_all.png", width = 8500, height = 7000, res = 500) 
print(
  ST.DimPlot(data, 
             dims = 1:3,
             indices = 1:length(id),
             grid.ncol = 5,
             reduction = "umap3d_harmony_SCT",
             blend = T,
             pt.size = 1,
             cols = cscale, 
             show.sb = FALSE) + 
    theme(text = element_text(size = 20))
)
dev.off() 

out.dir <- paste(dims.dir, "HARMONY_SCT", data.type, sep = "/")
dir.create(out.dir, recursive = T)
setwd(out.dir)
png(file = "uwot_3d_spatial_all.png", width = 8500, height = 7000, res = 500) 
print(
  ST.DimPlot(data, 
             dims = 1:3,
             indices = 1:length(id),
             grid.ncol = 5,
             reduction = "uwot3d_harmony_SCT",
             blend = T,
             pt.size = 1, 
             cols = cscale, 
             show.sb = FALSE) + 
    theme(text = element_text(size = 20))
)
dev.off() 

#' plot and visualize cell type deconvolution results
NumSpotPlot.list <- list()
ProportionPlot.list <- list()
Idents(data) <- harmony.ident 
Idents(data) <- factor(Idents(data))
data[["cluster"]] <- Idents(data)
names(clust.col) <- levels(data)
Idents(data) <- "sample"            


p1 <- ggplot(data@meta.data, aes(x = sample, fill = cluster)) + 
  geom_bar() +
  scale_fill_manual(values = clust.col) +
  labs(x="Sample",
       y="Number of Spots",
       title = "",
       subtitle = "",
       caption  = "") +
  theme(plot.caption = element_text(hjust = 0, face= "italic"), 
        plot.title.position = "plot",
        plot.caption.position =  "plot",
        legend.position = "none")
#' stat on number 
getOption("digits")
options(digits = 1)
d1<-data@meta.data %>% 
  dplyr::count(sample, cluster) %>% 
  group_by(cluster) %>%
  summarize(min=min(n),
            max=max(n),
            mean=mean(n),
            median=median(n),
            std=sd(n))    
tbl1 <- gridExtra::tableGrob(d1, rows=NULL)
               
#' plot percentage on figure as well (option)
data@meta.data %>% 
  dplyr::count(sample, cluster) %>%       
  group_by(sample) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(sample, pct, fill=cluster) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = clust.col) +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  labs(x="Sample",
       y="Spot Proportions",
       title = "",
       subtitle = "",
       caption  = "") +
  theme(plot.caption = element_text(hjust = 0, face= "italic"), 
        plot.title.position = "plot",
        plot.caption.position =  "plot", 
        legend.position = "right") 


  p2 <- ggplot(data@meta.data, aes(x = sample, fill = cluster)) + 
    geom_bar(position = "fill") +
    scale_fill_manual(values = clust.col) +
    labs(x="Sample",
         y="Spot Proportions",
         title = "",
         subtitle = "",
         caption  = "") +
    theme(plot.caption = element_text(hjust = 0, face= "italic"), 
          plot.title.position = "plot",
          plot.caption.position =  "plot", 
          legend.position = "right") 
#' stat on proportions
  d2<-data@meta.data %>% 
    dplyr::count(sample, cluster) %>%       
    group_by(sample) %>%
    mutate(pct= prop.table(n) * 100) %>%
    group_by(cluster) %>%
    summarize(min=min(pct),
              max=max(pct),
              mean=mean(pct),
              median=median(pct),
              std=sd(pct))
  tbl2 <- gridExtra::tableGrob(d2, rows=NULL)
  
  
  p3 <- ggplot(data@meta.data, aes(x = cluster, fill = sample)) + 
    geom_bar() +
    scale_fill_manual(values = col_vector) +
    labs(x="Cluster",
         y="Number of Spots",
         title = "",
         subtitle = "",
         caption  = "") +
    theme(plot.caption = element_text(hjust = 0, face= "italic"), 
          plot.title.position = "plot",
          plot.caption.position =  "plot",
          legend.position = "none")
  
  
  p4 <- ggplot(data@meta.data, aes(x = cluster, fill = sample)) + 
    geom_bar(position = "fill") +
    scale_fill_manual(values = col_vector) +
    labs(x="Cluster",
         y="Spot Proportions",
         title = "",
         subtitle = "",
         caption  = "") +
    theme(plot.caption = element_text(hjust = 0, face= "italic"), 
          plot.title.position = "plot",
          plot.caption.position =  "plot", 
          legend.position = "right")
              
  
png(file = "clusters_sample_bar.png", width = 6000, height = 8000, res = 500)
    ggarrange(p1,NULL,p2,NULL, NULL,NULL,p3,NULL,p4,tbl1, NULL, tbl2, ncol = 3, nrow = 4, 
              labels = c("A","","B","","","","C","","D","","",""), 
              font.label = list(size = 20), 
              widths = c(1, 0.05, 1), 
              heights = c(1,0.05,1,0.6))
    dev.off()
