get_marker <- function(dataset, 
                       clusters, 
                       comparison,
                       ...) {
  de <- list()
  
  for (i in seq_along(1:length(clusters))) {
    
    d <- FindMarkers(dataset, 
                     subset.ident = clusters[i],
                     group.by = comparison, 
                     ...) %>%
      cbind(cluster = clusters[i], gene = row.names(.))
    de[[i]] <- d
  }
  bind_rows(de)
}

get_conserved_marker <- function(dataset, 
                                 clusters, 
                                 group,
                                 ...) {
  de <- list()
  
  for (i in seq_along(clusters)) {
    
    d <- FindConservedMarkers(dataset, 
                              ident.1 = clusters[i],
                              grouping.var = group, 
                              ...) %>% 
      cbind(cluster = clusters[i], gene = row.names(.))
    
    if (is.null(d) == FALSE) {
      d$avg_log2FC <- 
        apply(d[,names(d) %in% colnames(d)[grepl("log2",colnames(d))]], 
              1, mean, na.rm = T)
      de[[i]] <- 
        d[c("cluster", "gene", "avg_log2FC", "max_pval", "minimump_p_val")]
    } else { 
      next 
      }
  }
  bind_rows(de)
}
