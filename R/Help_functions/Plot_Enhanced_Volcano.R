#' Generate volcano plot for DEA result
png(file = paste0(patient, "_volcano_", n, ".png"), 
    width = 8000, height = 8000, res = 500)
EnhancedVolcano(Clusters_markers_filt,
                lab = Clusters_markers_filt$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "patient",
                subtitle = names(de[n]),
                pCutoff = 0.01,
                FCcutoff = 1,                              
                pointSize = 3.0,
                labSize = 6.0,
                legendLabels=c('Not sig.',
                               'avg_Log2(FC)',
                               'p_adj_val',
                               'p_adj_val & avg_Log2(FC)')) +
  labs(x = "avg_Log2(Fold Change)", y = "-log10 (Adjust P_value)") +
  facet_wrap(. ~ cluster, ncol = 3)
dev.off()
