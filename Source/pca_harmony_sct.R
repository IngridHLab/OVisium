source("/home/minerva/VisiumST/gitHub/library.R")
source("/home/minerva/VisiumST/gitHub/directory.R")

file <- list.files(path = path, pattern = "\\Merged.sct.v2.rds$")
file.name <- gsub(".rds", "", file)
data <- readRDS(paste(path, file, sep = "/"))
data <- subset(data, subset = loupeA.ident != "Unknown")

set.seed(1220)
## Set number of dimension based on JackStraw and Elbow results;
npcs <- threshold # or threshold obtained from choose_pca_dims.R
dims <- 1:npcs

pca.reduction <- paste("PCA", npcs, sep = ".")
## Rerun PCA on SCTransformed data:
data <- RunPCA(data, npcs = npcs, assay = "SCT", 
               features = VariableFeatures(data))
data@misc$reductions.backup[[pca.reduction]] <- data[["pca"]]

## Run harmony and backup reduction
harmony.reduction <- paste("PCA", npcs, "Harmony.by.seq.id", sep = ".")
png(paste(file.name, "pca", npcs, "harmony.iteration.png", sep = "."), 
    width = 706, height = 504, res = 100)
print(
  data <- RunHarmony(data, reduction = "pca", assay.use = "SCT", 
                     group.by.vars = "seq_id", plot_convergence = TRUE)
)
dev.off()
data@misc$reductions.backup[[harmony.reduction]] <- data@reductions$harmony

setwd(path)
saveRDS(data, file = paste(file.name, "pca.harmony.rds", sep = "."))