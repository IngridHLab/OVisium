set.seed(1220)
 file.name <- "OVisium_SCT_merged"
 cluster.ident <- "harmony_SCT_res_0.6"
 feature.ident <- c("Variable_features_filt")
 #' Other features : 
 #' "Variable_a500UMI_features","Scaled_features","Scaled_a500UMI_features","All_features"
 #' Different list based on the filtering
 #' Test methods to identify DE genes
 test.ident <- c("MAST")
 markers.list.dir <- "FindMarkers/Patient_Pair"
 
 for (i in seq_along(feature.ident)) {
   for (t in seq_along(test.ident)) {
     input.dir <- paste(deg.dir, 
                        file.name, 
                        cluster.ident, 
                        feature.ident[[i]], 
                        test.ident[[t]],
                        markers.list.dir,
                        sep = "/")
     
     deg.files <- list.files(path = input.dir, pattern = "*DEA.list.rds")
      
    for (x in seq_along(anno.files)) {
      for (y in seq_along(deg.files)) {
          if (y > length(deg.files)) break 
          try( plot_compareCluster_patient(file      		= deg.files[[y]], 
          			                   anno      		= names[x], 
                         	        fun       		= "enricher",
                         	        dir       		= input.dir, 
                         	        TERM2GENE 		= anno.files[[x]],
                      		        pvalueCutoff 	= 0.01,
                      		        qvalueCutoff 	= 0.05)
                      )
          y <- y + 1
      }
    }
  }
  }


 for (i in seq_along(feature.ident)) {
   for (t in seq_along(test.ident)) {
     input.dir <- paste(deg.dir, 
                        file.name, 
                        cluster.ident, 
                        feature.ident[[i]], 
                        test.ident[[t]],
                        markers.list.dir,
                        sep = "/")
     
     deg.files <- list.files(path = input.dir, pattern = "*DEA.list.rds")
     
     for (x in seq_along(anno.files)) {
       for (y in seq_along(deg.files)) {
         if (y > length(deg.files)) break 
         try( plot_enrichBar_patient(file         = deg.files[[y]], 
                             anno         = names[x],
                             assay        = "enricher", 
                             dir          = input.dir,
                             TERM2GENE    = anno.files[[x]],
                             pvalueCutoff = 0.01,
                             qvalueCutoff = 0.05)
         )
         y <- y + 1
       }
     }
   }
 }


