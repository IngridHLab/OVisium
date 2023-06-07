source("/home/minerva/VisiumST/gitHub/library.R")
source("/home/minerva/VisiumST/gitHub/directory.R")

file <- list.files(path = path, pattern = "\\Merged.sct.v2.rds$")
file.name <- gsub(".rds", "", file)
data <- readRDS(paste(path, file, sep = "/"))


## Determine number of dims in the dataset.
data <- se.merged
file.name <- "STU.sct.merged"
#ribo_genes <- grep(pattern = "^RP[SL]", x = rownames(data), value = TRUE)
#mito_ribo_genes <- grep(pattern = "^MRP[SL]", x = rownames(data), value = TRUE)
#mito_genes <- grep(pattern = "^MT-", x = rownames(data), value = TRUE)
#pseu_genes <- grep(pattern = "^MTRNR", x = rownames(data), value = TRUE)
#linc_genes <- grep(pattern = "^LINC", x = rownames(data), value = TRUE)
#VariableFeatures(data) <- setdiff(VariableFeatures(data), c(mito_ribo_genes, mito_genes, pseu_genes, linc_genes))
#VariableFeatures(data) <- setdiff(VariableFeatures(data), c(mito_ribo_genes, mito_genes))
set.seed(1220)
data <- RunPCA(data, assay = "SCT")
data@misc$reductions.backup[["pca.50"]] <- data[["pca"]]

## Note, JackStraw can not be performed on SCT data and it takes very long time.
## In this Rscript, we only run the elbowplot to estimate the number of dims.
setwd(dims.dir)
png(paste(file.name, "elbow.png", sep = "."), width = 1200, height = 1000, 
    res = 300)
print(
  ElbowPlot(data, ndims = 50)
  )
dev.off()

## Quantitative elbowplot:
## Determine percent of variation associated with each PC
pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100

## Calculate cumulative sum of variation for each PC
cumu <- cumsum(pct)

## Determine which PC exhibits cumulative variation greater than 90% and 
## variation associated with the PC as less than 5%
co1 <- which(cumu > 90 & pct < 5)[1]

## Determine the difference between variation of a PC and subsequent PC less 
## than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), 
            decreasing = T)[1] + 1

## Minimum of the two calculation
threshold <- min(co1, co2)

## Create a data frame with values
plot_df <- data.frame(Variation = pct, 
                      Cumu.Variation = cumu,
                      PCs = 1:length(pct))

## Elbow plot to visualize
png(paste(file.name, "quant.elbow.png", sep = "."), 
    width = 1200, height = 1000, res = 300)
print(
  ggplot(
    plot_df, 
    aes(Cumu.Variation, Variation, label = PCs, color = PCs <= threshold)) + 
    geom_text(size = 3) + 
    geom_vline(xintercept = 90, color = "red") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "red") +
    theme_classic() + 
    theme(axis.text=element_text(color='black')) +
    theme(strip.background = element_blank()) +
    theme(legend.position = c(0.8, 0.8))
)
dev.off()

## check the first 30 PCs to confirm the signal below/above threshold
png(paste(file.name, "pca.heatmap.png", sep = "."), 
    width = 4000, height = 5000, res=400)
print(
  DimHeatmap(data, dims = 1:18, reduction = "pca", ncol = 3)
)
dev.off()
