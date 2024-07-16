#' Prepare gene annotation using emsembl id
library(fs)
home <- path_home()
source(paste(home, "OVisium/R/Help_functions/Load_Library.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Create_Directory.R", sep = "/"))
source(paste(home, "OVisium/R/Help_functions/Get_Gene_Annotation.R", sep = "/"))

#' Extract all 
agg.file <- "AGG_BRCA_18_20231018"
features <- 
  read_table(paste(sr.dir, "Outs", agg.file, "outs", 
                   "filtered_feature_bc_matrix", "features.tsv.gz", sep = "/"), 
             col_names = F)[c(1,2)] %>% 
  dplyr::rename(ENSEMBL = 1, SYMBOL = 2)

#' check unique emsembl id with same gene name, total 10 in our dataset
  dup <- features[which(duplicated(features$SYMBOL) == T),]

#' To convert gene id to entrez id for GSEA: Using the `bitr()` in the
#' `clusterProfiler` package to convert gene id or ensembl id to NCBI
#' Entrez id based on org.Hs.eg.db:
  ens.to.ent <- bitr(features$ENSEMBL, 
                   fromType = "ENSEMBL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

#' Warning message: In bitr(features$ENSEMBL, fromType = "ENSEMBL", 
#' toType = "ENTREZID",  : 34.82% of input gene IDs are fail to map...

#' Using `mapIds()` in `AnnotationDbi` package which contain no duplicates:
  ens.to.ent.2 <- mapIds(org.Hs.eg.db, 
                         keys=features$ENSEMBL, 
                         keytype = "ENSEMBL", 
                         column = "ENTREZID") %>% 
    as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    dplyr::rename(ENSEMBL = 1, ENTREZID = 2) %>% 
    na.omit()

#' Re-run `bitr()` on the remaining through gene name and obtained additional 
#' 1456 genes.
features.ent <- inner_join(features, ens.to.ent.2) 
diff <- anti_join(features, ens.to.ent.2)
sym.to.ent <- bitr(diff$SYMBOL, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
features.ent.diff <- inner_join(diff, sym.to.ent) 

#' Another way to convert gene id is using `biomaRt` package based Ensembl data. 
#' Many AL\^ genes have different gene name. Note, sometime it can be issue to 
#' connect to the Ensembl server:
#' Version 110 is used here:
ensembl <- useMart("ensembl")
listMarts(ensembl) # to check biomart version
bm <- 
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), 
        mart = useDataset("hsapiens_gene_ensembl", ensembl)) 
ens.to.ent.3 <- inner_join(diff, bm[c(1,3)], 
                           by = c("ENSEMBL" = "ensembl_gene_id")) %>% 
  na.omit() %>% 
  dplyr::rename(ENTREZID = 3) %>%  
  mutate_if(is.numeric, as.character)
features.ent.diff.2 <- inner_join(diff, ens.to.ent.3) 


#' `Annotables` package with grch37 or grch38 which is similar to `biomaRt` 
#' but data is loaded locally so need to keep track on the version:
ens.to.ent.4 <- 
  inner_join(diff, grch38[c(1,2)], by = c("ENSEMBL" = "ensgene")) %>% 
  na.omit() %>% 
  dplyr::rename(ENTREZID = 3) %>% 
  mutate_if(is.numeric, as.character)
features.ent.diff.3 <- inner_join(diff, ens.to.ent.4)


#' Use `AnnotationHub` to obtain ensembl
diff.2 <- 
  inner_join(diff, annotations, by = c("ENSEMBL" = "gene_id"))[c(1,2)] %>% 
  na.omit() 
sym.to.ent.2 <- bitr(diff.2$SYMBOL, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
features.ent.diff.4 <- inner_join(diff.2, sym.to.ent.2) 

# Final merging retrieve additional 1844 Entrez ids:
features.ent.d.all <- bind_rows(features.ent.diff, features.ent.diff.2, 
                                features.ent.diff.3, features.ent.diff.4)
features.ent.d.all <- 
  features.ent.d.all[which(duplicated(features.ent.d.all$ENSEMBL) == F), ]
features.ent.all <- bind_rows(features.ent, features.ent.d.all)
setwd(deg.dir)
write.table(features.ent.all, file = "features_EntrezID.csv", row.names = F, 
            quote = F, sep = ",")
save.image("features.covert.entrez.RData")
saveRDS(features.ent.all, "features_EntrezID.rds")
