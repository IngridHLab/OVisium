source("/home/minerva/VisiumST/gitHub/library.R")
source("/home/minerva/VisiumST/gitHub/directory.R")

path <- path
file <- list.files(path = path, pattern = "\\.dims.rds$")
file.name <- gsub(".rds", "", file)
data <- readRDS(paste(path, file, sep = "/"))

set.seed(1220)
nl <- NoLegend()
na <- NoAxes()
theme <- theme(plot.title = element_text(vjust = 0.8, 
                                         hjust = 0.5, size = rel(2)))

## Critical step to normalize the counts before DEA.
## Given a merged object with multiple SCT models, this function uses minimum 
## of the median UMI (calculated using the raw UMI counts) of individual 
## objects to reverse the individual SCT regression model using minimum of 
## median UMI of all as the sequencing depth covariate. The counts slot of the SCT 
## assay is replaced with re-corrected counts and the data slot is replaced with 
## log1p of re-corrected counts.
data <- PrepSCTFindMarkers(data, assay = "SCT")

setwd(deg.dir)
## DEG from all graph-based clusters
npcs <- 14
dims <- 1:npcs 
resolution <- 0.3
harmony.ident <- paste("PCA", npcs, "Harmony", "SNN", "res", resolution, sep = ".")
Idents(data) <- data[[harmony.ident]]


all.markers <- FindAllMarkers(data, assay = "SCT", 
                              logfc.threshold = 0.25, min.pct = 0.1)
write.table(all.markers, 
            file = paste(file.name, harmony.ident, "all.deg.markers.csv", sep = "."), 
            sep =",", quote = F, row.names = F, col.names = T)


## Identification of Spatially Variable Features. Search for features 
## exhibiting spatial patterning in the absence of pre-annotation. 
## The default method markvariogram, moransi. Other tools SpatialID and Splotch 
## are not supported yet. nfeature, Max number of features as the top spatially 
## variable; If only compute on 1:1000 features: 
## features = VariableFeatures(sampleID)[1:1000] otherwise all features; 
## This is a heavy step and don't run it on markdown mode!
data <- FindSpatiallyVariableFeatures(
  data, assay = "SCT", slot = "scale.data", 
  features = VariableFeatures(data)[1:1000], 
  selection.method = "markvariogram") 

## The output is simply a list without pvalue. 
all.spatial.markers <- tibble(
  SpatiallyVariableFeatures(data, selection.method = "markvariogram")) %>% 
  dplyr::rename("gene" = 1)

write.table(all.spatial.markers, 
            file = paste(file.name, harmony.ident, "all.spatial.markers.csv", sep = "."), 
            sep =",", quote = F, row.names = F, col.names = T)

## Save data
setwd(path)
saveRDS(data, file = paste(file.name, "DEA.rds", sep = "."))

## We can plot the top features from each clusters into a heatmap for easy
## comparison: filter away MT- and RPL genes, filter all the markers which have 
## more than 1.5 fold up/down-regulation and padjust below 0.05.  
top.markers <- all.markers %>% 
  group_by(cluster) %>%
  dplyr::filter(!grepl("^MT-", gene), 
                !grepl("^RP[SL]", gene),
                !grepl("^MTRNR", gene),
                !grepl("^LINC", gene),
                abs(avg_log2FC) > 0.58,
                p_val_adj < 0.05) %>% 
  top_n(n = 10, wt = abs(avg_log2FC))

write.table(top.markers, 
            file = paste(file.name, harmony.ident, "top.deg.csv", sep = "."), 
            sep =",", quote = F, row.names = F, col.names = T)


## Select 2 features with most abs fold changes from each cluster and remove 
## any feature with MT or pseudogene:
most <- top.markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 2, order_by = avg_log2FC)
most.features <- most[which(duplicated(most$gene) == FALSE),]$gene 
png(paste(file.name, harmony.ident, "most.png", sep = "."), 
    width = 1500, height = 1500)
print(
  FeaturePlot(data, features = most.features, keep.scale = "all") & nl & na
)
dev.off()
write.table(most, 
            file = paste(file.name, harmony.ident, "most.deg.csv", sep = "."), 
            sep =",", quote = F, row.names = F, col.names = T)



## Lets cross the spatial markers with the DEG list
spatial.de.markers <- inner_join(all.markers, all.spatial.markers)
top.spatial.de.markers <- spatial.de.markers %>% 
                          group_by(cluster) %>%
                          dplyr::filter(!grepl('MT-', gene), 
                                        !grepl('RPL', gene),
                                        !grepl('MTRNR2L12', gene),
                                        abs(avg_log2FC) > 0.58,
                                        p_val_adj < 0.05) %>% 
                          top_n(n = 10, wt = abs(avg_log2FC))

write.table(top.spatial.de.markers, 
            file = paste(file.name, harmony.ident, "top.spatial.deg.csv", sep = "."), 
            sep =",", quote = F, row.names = F, col.names = T)

most.spatial <- top.spatial.de.markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 2, order_by = avg_log2FC) 
most.spatial.features <- most.spatial[which(duplicated(most.spatial$gene) == FALSE),]$gene 
png(paste(file.name, harmony.ident, "most.spatial.png", sep = "."), 
    width = 1500, height = 1500)
print(
  FeaturePlot(data, features = most.spatial.features, keep.scale = "all") & nl & na
)
dev.off()
write.table(most.spatial, 
            file = paste(file.name, harmony.ident, "most.spatial.deg.csv", sep = "."), 
            sep =",", quote = F, row.names = F, col.names = T)


## Visualize overlap genes using venndiagram
s1 <-  all.markers %>% 
  dplyr::filter(!grepl('MT-', gene), 
                !grepl('RPL', gene),
                !grepl('MTRNR2L12', gene))
s2 <-  all.spatial.markers %>% 
  dplyr::filter(!grepl('MT-', gene), 
                !grepl('RPL', gene),
                !grepl('MTRNR2L12', gene))

studies = list( S1=s1$gene %>% unique(),
                S2=s2$gene)
ol = calculate.overlap(x = studies)
ol_size=sapply(ol, length)

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = studies,
  category.names = c("DEG across all cluster" , "Spatial variable genes"),
  filename = paste(file.name, harmony.ident, "spatial.deg.venn.tiff"),
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  # Numbers
  cex = 0.6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  #rotation = 1
)

## alternative using vennDetail package for individual clusters
venn <- venndetail(list(DEG = all.markers$gene, 
                        Spatial = all.spatial.markers$gene))
getFeature(venn, subset = "Shared", rlist = list(dataframe1, dataframe2))





## Plot multibarheatmap provided by Scillus
## https://github.com/satijalab/seurat/issues/2201

## scale all features for better heatmap
data <- ScaleData(data, features = rownames(data))

## Create some new idents for annotation of heatmap
Idents(data) <- "loupeA.ident"
data[["Morphology"]] <- gsub(pattern = "_\\d+$", "", data@active.ident)
Idents(data) <- "Morphology" 
Idents(data) <- factor(Idents(data), 
                       levels = c("FTE","Mix", "FTE Stroma","Stroma"))
data[["Morphology"]] <- Idents(data)

Idents(data) <- "orig.ident"
data[["Sample"]] <- Idents(data) 
data[["Mutation"]] <- gsub(pattern = "PO9_A1", "BRCA2", data@active.ident)
data[["Mutation"]] <- gsub(pattern = "PO9_B1", "BRCA2", data$Mutation)
data[["Mutation"]] <- gsub(pattern = "PO3.+$", "BRCA2", data$Mutation)
data[["Mutation"]] <- gsub(pattern = "PO42", "BRCA2", data$Mutation)
data[["Mutation"]] <- gsub(pattern = "^PO.+$", "BRCA1", data$Mutation)
Idents(data) <- "Mutation"
Idents(data) <- factor(Idents(data),
                       levels = c("BRCA1", "BRCA2"))
data[["Mutation"]] <- Idents(data)

Idents(data) <- "orig.ident"
data[["Age"]] <- gsub(pattern = "PO2$", "35+", data@active.ident)
data[["Age"]] <- gsub(pattern = "PO2_2re", "35+", data$Age)
data[["Age"]] <- gsub(pattern = "PO28$", "35+", data$Age)
data[["Age"]] <- gsub(pattern = "PO45$", "35+", data$Age)
data[["Age"]] <- gsub(pattern = "PO3.+$", "50+", data$Age)
data[["Age"]] <- gsub(pattern = "PO42$", "50+", data$Age)
data[["Age"]] <- gsub(pattern = "^PO.+$", "40+", data$Age)
Idents(data) <- "Age"
Idents(data) <- factor(Idents(data),
                       levels = c("35+", "40+", "50+"))
data[["Age"]] <- Idents(data)

## Plot regular seurat heatmap
morph.cols <- c("#8B7500", "#FFD700", "#FFBBFF", "#9400D3")
sample.cols <- scCustomize_Palette(num_groups = 20, ggplot_default_colors = T)
mut.cols <- c("#63B8FF", "#B3EE3A")
age.cols <- c("#B03060", "#FFA07A", "#FFE4B5")

pdf(paste(file.name, harmony.ident, "deg.hm.pdf", sep = "."), 
    width = 15, height = 15, onefile = F) 
plot_heatmap(
  data,
  markers = unique(top.markers$gene),
  sort_var = c(harmony.ident, 
               "Morphology"),
  anno_var = c(harmony.ident,
               "Morphology",
               "percent.Mito", 
               "nCount_Spatial", 
               "Mutation", 
               "Age",
               "Sample"),
  anno_colors = list("Paired",
                     morph.cols,
                     "Greens", 
                     c("blue","white","red"), 
                     mut.cols,
                     age.cols,
                     sample.cols),
  hm_colors = c("purple","black","yellow"))
dev.off()

## Seurat heatmap
Idents(data) <- harmony.ident
png(paste(file.name, harmony.ident, "deg.hm.seurat.png", sep = "."), 
    width = 1200, height = 1200)
print(
DoHeatmap(data, features = unique(top.markers$gene), size = 3)
)
dev.off()

## spatial visualization
SpatialFeaturePlot(data, features = top.spatial.markers, 
                   ncol = 10, alpha = c(0.1, 1))

## Save data
setwd(path)
saveRDS(data, file = paste(file.name, harmony.ident, "rds", sep = "."))

