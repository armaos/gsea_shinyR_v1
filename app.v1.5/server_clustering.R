print("start server _clustering observer")
#fdr_filter = 0.05
#fc_filter_low = 1.5

print(paste("server clustering MSI df dim: ",dim(rv$msi_df), collapse = ", "))

degs <<- subset(rv$msi_df %>%
                select(primary_name, contains("FDR"), contains("logFC"))  %>% 
                gather(k, v , -primary_name) %>% 
                separate(k, c("cond", "value_type"), sep = "[.]") %>% 
                spread(value_type, v) %>% 
                filter(abs(logFC) >= fc_filter_low & FDR <= fdr_filter))$primary_name

degs<<-unique(degs)
#print("DEGS")
print(paste(head(degs), collapse = "," ))

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



