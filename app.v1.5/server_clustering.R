print("start server _clustering observer")
#fdr_filter = 0.05
#fc_filter_low = 1.5

print(paste("server clustering MSI df dim: ",dim(rv$msi_df), collapse = ", "))

# Reshape DE statistics: select only FDR and logFC columns, split condition from statistic
# Pattern: {condition}.{statistic} (e.g., "Low_Medium.logFC", "Low_Medium.FDR")
# Be explicit: only select columns that END with ".FDR" or ".logFC"
de_stat_cols <- colnames(rv$msi_df)[endsWith(colnames(rv$msi_df), ".FDR") | endsWith(colnames(rv$msi_df), ".logFC")]

print("DE stat columns found:")
print(de_stat_cols)

msi_de_stats <- rv$msi_df %>% rename_with(~ sub("^.*\\.(FDR|logFC|logCPM|LR|PValue|)$", "\\1", .x))
msi_de_stats <- msi_de_stats %>%
  mutate(across(any_of(c("FDR", "logFC")), ~ as.numeric(.x)))
  

# Filter for DEGs
degs <<- msi_de_stats %>%
  filter(abs(logFC) >= fc_filter_low & FDR <= fdr_filter) %>%
  pull(primary_name)

print(paste("fc_filter_low: ", fc_filter_low))
print(paste("fdr_filter: ", fdr_filter))
print(paste("top msi (reshaped):"))
print(head(msi_de_stats %>% arrange(logFC)))

print(paste("filter msi (after filter):"))
filtered_degs <- msi_de_stats %>%
  filter(abs(logFC) >= fc_filter_low & FDR <= fdr_filter) %>%
  arrange(desc(abs(logFC)))
print(head(filtered_degs))

print(paste("degs after pull (before unique):"))
print(paste("  length:", length(degs)))
print(paste("  class:", class(degs)))
print(paste("  first few:", paste(head(degs, 5), collapse = ", ")))

degs <<- unique(degs)
print(paste("degs after unique:"))
print(paste("  length:", length(degs)))
print(paste("  first few:", paste(head(degs, 5), collapse = ", ")))

degs <<- filter_by_protein_type(degs, protein_type)
print(paste("degs after protein type filter:"))
print(paste("  length:", length(degs)))
print(paste("  first few:", paste(head(degs, 5), collapse = ", ")))
#print("DEGS")

rv$all_together_tcpm <- rv$msi_df %>% 
  filter(primary_name %in% degs)
  #separate_rows(Name, sep = ",") 


  
print(paste("server clustering all_together_tcpm df dim: ",dim(rv$all_together_tcpm), collapse = ", "))

rv$all_together_tcpm <- as.data.frame(rv$all_together_tcpm)
rownames(rv$all_together_tcpm) <- rv$all_together_tcpm$primary_name 



rv$all_together_tcpm <- rv$all_together_tcpm %>%
  select(-all_of(c("primary_name",annotations_columns))) 



rv$all_together_tcpm <- rv$all_together_tcpm %>% 
  select(-contains(c("logCPM", "logFC", "LR", "PValue", "PValue", "FDR", "genelabels", "change")))  %>%
  select(-group_id)

cols = colnames(rv$all_together_tcpm) 
rv$all_together_tcpm <- rv$all_together_tcpm %>% 
  select(all_of(c(cols[endsWith(cols, ".Low")],
                  cols[endsWith(cols, ".Medium")], 
                  cols[endsWith(cols, ".High")])), 
         everything()) 
print(colnames(rv$all_together_tcpm))

rv$all_together_tcpm <- rv$all_together_tcpm %>% filter_all(all_vars(!is.na(.)))

rv$all_cluster_heatmap <- tryCatch(
  {

    pheatmap(rv$all_together_tcpm,
             fontsize_row = 3,
             angle_col=90,
             scale="row",
             cluster_cols=F,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             main = paste0("Expression of all DEGs (",length(degs), " with FDR < " ,fdr_filter ," and |logFC| > ",fc_filter_low ," ) for each sample"),
             width=8,
             height=10,
             show_rownames=F)
  },
  error=function(cond) {
    empty_image()
  }
)

rv$show_clustering_mout = T
#rv$show_norm_and_de_mout = T
#shinyjs::show("menu_norm_and_de_panel")
print("server_clustering observer OK")



