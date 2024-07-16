library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))

file.name <- "OVisium_SCT_merged"
cluster.ident <- "harmony_SCT_res_0.6"
#' Extract QC result from the merged data
qc <- readRDS(paste(rds.dir, "OVisium_SCT_merged.rds", sep = "/"))

#' raw RNA data
x<-qc@meta.data$nCount_RNA
x<-qc@meta.data$nFeature_RNA
x <-qc@meta.data$Mito.percent
x <-qc@meta.data$Ribo.percent
x <-log10(data$nFeature_RNA) / log10(data$nCount_RNA)
x <-qc@meta.data$log10FeaturesPerUMI %>% na.omit()
gene_attr <- data.frame(nUMI = Matrix::rowSums(data@assays$RNA@counts), 
                            nSpots = Matrix::rowSums(data@assays$RNA@counts > 0))
x<-gene_attr$nUMI
x<-gene_attr$nSpots

#' SCT data
c.features <- row.names(data@assays[["SCT"]]@scale.data) 
data.subset<-subset(data, feature=c.features) #' use to check vfeature 
x <-data@meta.data$nFeature_SCT
x <-data@meta.data$nFeature_SCT
x <-data@meta.data$Mito.percent
x <-data@meta.data$Ribo.percent
x <-log10(data$nFeature_SCT) / log10(data$nCount_SCT) %>% na.omit()
gene_attr_sct <- data.frame(nUMI = Matrix::rowSums(data@assays$SCT@counts), 
                            nSpots = Matrix::rowSums(data@assays$SCT@counts > 0))
x<-gene_attr_sct$nUMI
x<-gene_attr_sct$nSpots
row.names(data@assays[["SCT"]]@scale.data)

#' calculate statistic value
median(x)
mad(x)
IQR(x)
sd(x)


#' Extract info from the DEA result
out.dir<-paste(deg.dir, file.name, cluster.ident, "Variable_features", "MAST", 
               sep = "/")
setwd(out.dir)
all.markers.filt <- readr::read_csv("Clusters.markers.filt.csv")

pos.FTE.markers.filt <- all.markers.filt %>% 
  dplyr::filter(avg_log2FC > 1,
                cluster == "0_FTE" | cluster == "1_FTE")
neg.FTE.markers.filt <- all.markers.filt %>% 
  dplyr::filter(avg_log2FC <= 1,
                cluster == "0_FTE" | cluster == "1_FTE")

pos.Stroma.markers.filt <- all.markers.filt %>% 
  dplyr::filter(avg_log2FC > 1,
                cluster != "0_FTE" & cluster != "1_FTE" & cluster != "3_Mix")
neg.Stroma.markers.filt <- all.markers.filt %>% 
  dplyr::filter(avg_log2FC <= 1,
                cluster != "0_FTE" & cluster != "1_FTE" & cluster != "3_Mix")


top <- pos.markers.filt %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 100)

#' DEA between two FTE clusters
FTE0vs1.markers.filt <- readr::read_csv("Clusters.0vs1.markers.filt.csv")
FTE.0.pos.markers.filt <- FTE0vs1.markers.filt %>% 
  dplyr::filter(avg_log2FC > log2(2))

FTE.1.markers.filt <- FTE0vs1.markers.filt %>% 
  dplyr::filter(avg_log2FC < log2(0.75))

#' DEA between two Stomal clusters
S10vs11.markers.filt <- readr::read_csv("Cluster.10vs11.markers.filt.csv")
S.10.pos.markers.filt <- S10vs11.markers.filt %>% 
  dplyr::filter(avg_log2FC > log2(2))

S.11.markers.filt <- S10vs11.markers.filt %>% 
  dplyr::filter(avg_log2FC < log2(0.5))
