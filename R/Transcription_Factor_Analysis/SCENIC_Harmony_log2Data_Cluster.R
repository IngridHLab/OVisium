#' Initial Required
devtools::install_github("aertslab/AUCell")
#' RcisTarget identifies transcription factor binding motifs (TFBS) over-represented on a gene list:
devtools::install_github("aertslab/RcisTarget")
devtools::install_github("aertslab/GENIE3")
install.packages("~/Downloads/hdf5r_1.3.8.tar.gz", repos = NULL, type = "source")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC")

#' dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                          "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather",
                          "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                          "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather")
             # mc9nr: Motif collection version 9: 24k motifs)

for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}


library(fs)
home <- path_home()
source(paste(home, "OVisium/manuscript/gitHub/1_library.R", sep = "/"))
source(paste(home, "OVisium/manuscript/gitHub/2_directory.R", sep = "/"))
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCopeLoomR)
library(SCENIC)
setwd("~/OVisium/DE_functional_analysis/SCENIC")

## Get data from sce object
load(paste0(rds.dir, "/Variable_features_filt_SCT_log2counts+1_harmony.RData"))
exprMat <- as.matrix(GetAssayData(data.sub.filt, slot = "counts", assay = "SCT"))  
cellInfo <- data.frame(CellType=Idents(data.sub.filt))


# cellInfo$nGene <- colSums(exprMat>0)
head(cellInfo)
cellInfo <- data.frame(cellInfo)
cbind(table(cellInfo$CellType))
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")


# Color to assign to the variables (same format as for NMF::aheatmap)
clust.col <- 
  c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#F0027F","#BF5B17",
             "#666666","#1B9E77","#D95F02","#386CB0","black")
             
names(clust.col) <-c("0_FTE","1_FTE", "3_Mix", "2_Stroma",
                     "5_Stroma", "6_Stroma", "7_Stroma",
                     "8_Stroma", "9_Stroma", "10_Stroma",
                     "11_Stroma")
colVars <-list(CellType = clust.col)
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

### Initialize settings
library(SCENIC)
library(RcisTarget)
org= "hgnc"
dbdir <- "cisTarget_databases"
myDatasetTitle <- "SCENIC Harmony SCTcounts"
data("defaultDbNames")
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org= org, 
                                  dbDir= dbdir, 
                                  dbs= dbs, 
                                  datasetTitle= myDatasetTitle, 
                                  nCores=16) 
motifAnnotations_hgnc <-motifAnnotations
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)

## interesting genes
interestingGenes <- c("PAX8", "SOX17", "AGR2", "SOX2")
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
#[1]  2386 29791

#' Run correlation
runCorrelation(exprMat_filtered, scenicOptions)

#' Run Genie3 to infer potential TF targets
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 

### Identify cell states: 
# 3. Score GRN (regulons) in the cells (with AUCell) 
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

###### Save the analysis
saveRDS(scenicOptions, file="scenicOptions_240516.Rds")

## Binarize the network activity (regulon *on/off*)
#' Note: rbokeh was removed from BioConductor and AUCell was not allowed to depend on it anymore, so the plotTsne_AUCellApp viewer was removed too
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)

#' The follow step is skipped
# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions@settings$devType="png"

#' 4.2 Binarize the network activity (regulon on/off). 
#' scenicOptions@settings$devType="png" if save as png instead of pdf
#' This step is useful to reduce batch effect such as number of genes and total counts 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions, exprMat = exprMat_log)
#' Binary regulon activity: 64 TF regulons x 29785 cells.
#' (94 regulons including 'extended' versions)
#' 57 regulons are active in more than 1% (297.85) cells.
#' more about data binarization: https://medium.com/@brijesh_soni/topic-13-binarization-ca1059d8c1ce

## 4.3 Clustering / dimensionality reduction on the regulon activity (optional)
# 4.3.1 set number of PCs
nPcs <- c(5, 15, 50)

# 4.3.2 Calculates the t-SNE based on the regulon activity
#' Perplexity for the t-SNE (can be several values, see examples)
scenicOptions@settings$seed <- 1220 
fileNames <- tsneAUC(scenicOptions, aucType = "AUC", nPcs = nPcs, perpl = c(5,15,50))
perpl<-getSettings(scenicOptions,"aucell/tsne_perpl")


#' plot tSNE
par(mfrow=c(3,3))
plotTsne_compareSettings(fileNames[1], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames[2], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames[3], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames[4], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames[5], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames[6], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames[7], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames[8], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames[9], scenicOptions, showLegend = F, varName = "CellType", cex = .5)

# 4.3.3 Calculates the t-SNE using only "high-confidence" regulon
fileNames.OHC <- tsneAUC(scenicOptions, aucType = "AUC", nPcs = nPcs, perpl = c(5,15,50), onlyHighConf = T, filePrefix = "int/tSNE_OHC")
fileNames.OHC 

#' plot tSNE
par(mfrow=c(3,3))
plotTsne_compareSettings(fileNames.OHC[1], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames.OHC[2], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames.OHC[3], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames.OHC[4], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames.OHC[5], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames.OHC[6], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames.OHC[7], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames.OHC[8], scenicOptions, showLegend = F, varName = "CellType", cex = .5)
plotTsne_compareSettings(fileNames.OHC[9], scenicOptions, showLegend = F, varName = "CellType", cex = .5)


# 4.3.4 save default setting for t-SNE
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 15
scenicOptions@settings$defaultTsne$perpl <- 50
saveRDS(scenicOptions, file="scenicOptions_240607.Rds")


## 5, Exploring and interpreting the results
# 5.1 Identify cell states and their regulators
# 5.1.1 Projection the AUC and TF expression on t-SNEs
library(SCENIC)
loadInt(scenicOptions) # import all rds files in the int folder
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
exprMat <- readRDS("~/OVisium/DE_functional_analysis/SCENIC/exprMat.rds")

tsneFileName(scenicOptions)
par(mfrow=c(4,4))
plotTsne_compareSettings(tsneFileName(scenicOptions), scenicOptions, showLegend=T, varName = "CellType", cex = .5) 

#' show TF results
#' https://github.com/aertslab/SCENIC/issues/72
#' regarding what is 'extended' and low-confident motif
#' bracket show how many targeted genes 
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))][c("ETS2 (131g)", "FOS (12g)", "JUN (28g)", "CREM (56g)", "NFIL3 (39g)", "CEBPB (13g)", "KLF6 (16g)", "ATF3 (40g)", "CEBPD (16g)", "BCL3 (30g)", "IRF1 (32g)", "RELB (82g)", "NFKB2 (53g)")], plots="Expression")

#' or
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("ETS2", "FOS", "JUN", "CREM", "NFIL3", "CEBPB", "KLF6", "ATF3", "CEBPD", "BCL3", "IRF1", "RELB", "NFKB2")],], plots="Expression")

# 5.1.2 Density plot to detect most likely stable states (high density area)
library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col = brewer.pal(9, "YlOrBr"), axes = F)
contour(dens2d, add = T, nlevels = 5, drawlabels = F)


# 5.2 Regulon targets
#' Output file : Step2_regulonTargetsInfo

#' subset regulonTargetsInfo before exporting it as HTML
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo") 
tableSubset <- regulonTargetsInfo[TF %in% c("SOX17", "SOX4") & highConfAnnot == T]
viewMotifs(tableSubset, options=list(pageLength=10))

#' Targets for interested TFs
regulons.all <- loadInt(scenicOptions, "regulons")
regulons.all[c("ETS2", "FOS")]

#' Only regulons with 10 genes or more are scored with AUCell
regulons.aucell <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons.aucell))))

# 5.3 Regulon motifs
#' The full list of TF motifs: SCENIC/output/Step2_MotifEnrichment.tsv
#' Showing motifs for interested TF
motifEnrichment_selfMotifs_wGenes <-loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset.motifs <- motifEnrichment_selfMotifs_wGenes[highlightedTFs %in% c("SOX17", "SOX4")] #JUN SOX4 FOSB
viewMotifs(tableSubset.motifs)

# 5.4 Combine SCENIC analysis with other analysis tools
#' Regulators for known cell types or clusters
# 5.4.1 ComplexHeatmap to show average regulon activity by cell clusters
library(ComplexHeatmap)
cellInfo <- readRDS("~/OVisium/DE_functional_analysis/SCENIC/cellInfo.rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                    function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_scaled <- t(scale(t(regulonActivity_byCellType), scale = T, center = T))
ComplexHeatmap::Heatmap(
  regulonActivity_byCellType_scaled, name = "Regulon Activity"
)

#' Top postive Regulators for each cell type: Relative Activity
topRegulators <- reshape2::melt(regulonActivity_byCellType_scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity > 0),]
viewTable(topRegulators)

#' Regulators targets and motifs
#' subset regulonTargetsInfo before exporting it as HTML
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo") 
regulonActivity_byCellType_name <- str_remove(rownames(regulonActivity_byCellType_scaled), "_extended")
regulonActivity_byCellType_name <- gsub(" .*", "", regulonActivity_byCellType_name)
tableSubset <- regulonTargetsInfo[TF %in% regulonActivity_byCellType_name]
viewMotifs(tableSubset, options=list(pageLength=50))


# 5.4.2 Binarized data
minPerc <- .4 # percentage of cells in a cluster with the regulon is active
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo) %in% colnames(binaryRegulonActivity)),, drop = F]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells),
                                                     cellInfo_binarizedCells$CellType),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop = F]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc) > 0),]
ComplexHeatmap::Heatmap(
  binaryActPerc_subset, name = "Binary Regulon Activity", col = c("white","pink","red")
)


#' >0.4 Regulators for each cell type: Relative Activity
topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity > minPerc),]
viewTable(topRegulators)

# 5.4.3 Cell-typen specific regulators based on regulon specificity score (RSS)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC))]
rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellInfo[colnames(regulonAUC), "CellType"])
rss_for_Plot <- plotRSS(rss)
plotly::ggplotly(rss_for_Plot$plot)

#' plot for individual cluster or cell type
par(mfrow=c(3,4))
plot_list <- list()
for (i in 1:11) {
   cluster.name <- levels(cellInfo$CellType)[i]
   p <- plotRSS_oneSet(rss, setName = cluster.name) 
   plot_list[[i]] <- p
}
gridExtra::grid.arrange(grobs = plot_list) 


sessionInfo() %>% capture.output(file="session_info.txt")
