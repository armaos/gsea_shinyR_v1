
# the annotation table 
setwd(wd)

#annotations_default = as.data.frame(read.delim2("Ka.3.1_annotation_table_test.tsv", header = T, sep = "\t", fill = T))
#annotations_default = as.data.frame(read.delim2("Ka.2_plus_3.2_annotation_table_test.tsv", header = T, sep = "\t", fill = T,  quote = ""))
annotations_default = as.data.frame(read.delim2("Ka.2_plus_3.2_annotation_table.tsv", header = T, sep = "\t", fill = T))


# the cluster
# clusters_deg <- left_join(
#   as.data.frame(read.table("clusters_DEGS/cluster_BD_vs_GD.csv", header = T, sep = ",", fill = T,  quote = "")) %>%
#     rename(cluster_BD_vs_GD = cluster) %>% select(kal_t, cluster_BD_vs_GD),
#   as.data.frame(read.table("clusters_DEGS/cluster_BL_vs_GL.csv", header = T, sep = ",", fill = T,  quote = "")) %>%
#     rename(cluster_BL_vs_GL = cluster) %>% select(kal_t, cluster_BL_vs_GL),
#   by = "kal_t"
#   ) %>%
#   left_join(
#     as.data.frame(read.table("clusters_DEGS/cluster_dol.csv", header = T, sep = ",", fill = T,  quote = "")) %>%
#       rename(cluster_D_vs_O_vs_L = cluster) %>% select(kal_t, cluster_D_vs_O_vs_L),
#     by = "kal_t"
#   ) %>%
#   left_join(
#     as.data.frame(read.table("clusters_DEGS/cluster_together.csv", header = T, sep = ",", fill = T,  quote = "")) %>%
#       rename(cluster_together = cluster) %>% select(kal_t, cluster_together),
#   by = "kal_t")

#annotations_default <- annotations_default %>% left_join(clusters_deg, by = c ("Name"= "kal_t"))


representatives = annotations_default %>% select(Name)
representatives = as.data.frame(read.delim2("represenatives.txt", header = T))
annotations_default <- annotations_default %>% 
  select(colnames(annotations_default)[!grepl("exp_", colnames(annotations_default))]) %>%
  select(-any_of(c("protein_product", "seqid", "start", "end", "strand"))) %>%
  mutate(Similar_to = ifelse(Similar_to == "-", Name,  paste(Similar_to, Name, sep = ","))) %>% 
  separate_rows(Similar_to, sep = "[,]") %>% 
  mutate(Name = ifelse(Similar_to == "-", Name, Similar_to)) %>%
  select(-Similar_to) %>%
  select_if(function(x) !(all(is.na(x)) | all(x=="-")))

annotations_columns = sort(colnames(annotations_default))
print("dim annotations_default")
print(paste(dim(annotations_default), collapse = " "))





# here AnnotationDbi::Ontology("GO:0006486") for go ontology grouping
##go2gene_all <- as.data.frame(read.table("kal_v3_all_GO.txt", header = T))
# go2gene_all <- annotations_default %>%
#   select(Name, Ontology_term , em_GOs, b2go_GO_id, ips_GO, arabThal_GO_ID) %>%
#   gather(key, val, c("Ontology_term" , "em_GOs", "b2go_GO_id", "ips_GO", "arabThal_GO_ID")) %>%
#   filter(val != "-") %>%
#   select(val, Name) %>%
#   separate_rows(val, sep = ",") %>%
#   unique
# go2gene_all %>% write.table("Ka.2_plus_3.1_GO_terms_manipulate.txt")
go2gene_all <- as.data.frame(read.table("Ka.2_plus_3.1_GO_terms_manipulate.txt", header = T))

go_terms =  rbind(
  as.data.frame(read.delim2("GOIDs_GOterms.txt", header = T, sep = "\t", stringsAsFactors = F, fill = T)) %>%
    mutate(ArabT = "0") %>%
    filter(!grepl("obsolete ", GO_term)),
  as.data.frame(read.delim2("GOIDs_GOterms_ArabT.txt", header = F, sep = " ", stringsAsFactors = F, fill = T)) %>%
    mutate(ArabT = "1") %>%
    rename(GO_ID = V2) %>%
    mutate(GO_namespace = ifelse(grepl("GO_MF_", V1), "molecular_function", ifelse(
      grepl("GO_BP_", V1), "biological_process", "cellular_component" ) ) ) %>%
    mutate(V1 = str_remove(V1, "GO_MF_AT_"),
           V1 = str_remove(V1, "GO_CC_AT_"),
           V1 = str_remove(V1, "GO_BP_AT_")
           ) %>%
    rename(GO_term = V1) %>%
    #mutate(GO_term = str_replace_all(GO_term, "-", " ")) %>%
    mutate(GO_term = str_replace_all(GO_term, "_", " ")) %>%
    mutate(GO_term = str_to_lower(GO_term))
)

demo_genes <- as.data.frame(read.table("genes_to_test.txt", header = T))

# these are the KEGG thresholds per node coming from the prepare_graph_nodes_table.py scirpt that reads all the nodes from the paths.
# koala_thresholds = as.data.frame(read.table("koala_genes_vector.txt"))
# koala_thresholds <- koala_thresholds %>%
#   group_by(V1 ) %>% 
#   summarise(V2 = max(V2))

#### read KO annotations (i have manipulated the data to load faster)
## ko_annots <- as.data.frame(read.table("Ka.2_plus_3.1_KO_terms.txt", header = T, sep = "\t"))
# ko_annots <- annotations_default %>% select(Name, koala_KEGG_ko , em_KEGG_ko) # this is better since i get the annotations directly from the table
# ko_annots = full_join(
#   ko_annots %>%
#     select(Name, koala_KEGG_ko) %>%
#     filter(koala_KEGG_ko != "-") %>%
#     separate_rows(koala_KEGG_ko, sep = "[,]") %>%
#     separate(koala_KEGG_ko, c("KO", "koala_threshold"), sep = " ") %>%
#     mutate(koala_threshold = str_remove(koala_threshold, "[(]"),
#            koala_threshold = str_remove(koala_threshold, "[)]"),
#            koala_threshold = str_remove(koala_threshold, "[%]"),
#            koala_threshold = as.numeric(koala_threshold)),
#   ko_annots %>%
#     select(Name, em_KEGG_ko) %>%
#     filter(em_KEGG_ko != "-") %>%
#     separate_rows(em_KEGG_ko, sep = "[,]") %>%
#     mutate(em_KEGG_ko = str_remove(em_KEGG_ko, "ko:"),
#            em = 1),
#   by = c("Name" = "Name", "KO" = "em_KEGG_ko")
# )
# ko_annots %>% write.table("Ka.2_plus_3.1_KO_terms_manipulate.txt", quote = F, row.names = F)

ko_annots <- as.data.frame(read.table("Ka.2_plus_3.1_KO_terms_manipulate.txt", header = T)) 
koala_thresholds <- suppressWarnings(ko_annots %>%
                                       group_by(KO) %>%
                                       summarise(koala_threshold = max(koala_threshold, na.rm = T),
                                                 em = max(em)) %>%
                                       mutate(koala_threshold = ifelse(is.infinite(koala_threshold), NA, koala_threshold)))


ko_symbols <- as.data.frame(read.delim2("KO_titles_new.txt", header = F, sep = ";", quote = "")) %>%
  mutate(V3 = str_remove(V3, " ")) %>%
  rename(Symbol = V2, Node_Description = V3)


ko_Pathway_terms <- as.data.frame(read.table("Ka_KEGGpathways.txt", header = F, sep = "\t", stringsAsFactors = F, 
                                             colClasses = c("character", "character")))
#### read kegg mapp annotations (i have manipulated the data to load faster
## ko_Pathway = "Ka.2_plus_3.1_KO_Pathways.txt"
## ko_Pathway <- as.data.frame(read.table(ko_Pathway, header = T))
# ko_Pathway <- annotations_default %>% select(Name, em_KEGG_Pathway , koala_KEGG_Pathway) # this is better since i get the annotations directly from the table
# ko_Pathway <- ko_Pathway %>%
#   filter(em_KEGG_Pathway != "-" | koala_KEGG_Pathway != "-" ) %>%
#   gather(key, val, c(em_KEGG_Pathway, koala_KEGG_Pathway)) %>%
#   separate_rows(val, sep = "[,]") %>%
#   filter(val != "-" & grepl("ko", val)) %>%
#   select(-key) %>%
#   unique
# ko_Pathway %>% mutate(val = str_replace(val ,"ko", "")) %>% write.table("Ka.2_plus_3.1_KO_Pathways_manipulate.txt", quote = F, row.names = F)
ko_Pathway <- as.data.frame(read.table("Ka.2_plus_3.1_KO_Pathways_manipulate.txt", header = T, colClasses = c("character", "character"))) %>%
  dplyr::select(val, Name)  %>% 
  mutate(val = str_replace(val ,"ko", "")) %>%
  filter(val %in% ko_Pathway_terms$V1)  

######


ko_nodes <- as.data.frame(read.table("KEGGmaps_and_nodes.txt", header = F, sep = "\t") ) %>% 
  rename(map = V1, nodes = V2) %>%
  mutate(map = str_replace(map, "ko", ""))

## aws variables
region = "us-east-1"
bucket_name= "forjazul.galaxy"
ora_tab_list <- NULL

custom_upset_terms_list <- c()

fdr_filter <- 0
fc_filter_low <- 0
fc_filter_high <- 0
degs <- c()
rv <- reactiveValues(gsea_plots = NULL,
                     IP = NULL, 
                     show_norm_and_de_mout = NULL,
                     show_clustering_mout = NULL,
                     show_gso_mout = NULL,
                     show_KeggScape_mout = NULL,
                     show_gsea_mout = NULL,
                     show_gsea_plots_mout = NULL,
                     show_ora_mout = NULL,
                     
                     runs_galaxy_df = data.frame(),
                     galaxy_outputs = data.frame(),
                     
                     condition = "",
                     conditions = c(),
                     
                     update_annotations = 0,
                     update_msi_df = 0,
                     update_msi_clustering = 0,
                     
                     norm_and_de_signal =0,
                     data = NULL, 
                     text = NULL, 
                     okplot = FALSE, 
                     #okplot_KEGGmap_ora = FALSE,
                     okplot_KEGGmap_gso = FALSE,
                     okplot_KEGGmap_gsea = FALSE,
                     okplot_gseaplot = FALSE,
                     map = "",
                     map_gsea = "",
                     go_term = "",
                     run_counter = 0, 
                     
                     gse = NULL,
                     gse_desc = NULL,
                     
                     # barplot1 = NULL,
                     # barplot2 = NULL,
                     # dotplot = NULL,
                     # emapplot = NULL,
                     # cnetplot = NULL,
                     # heatplot = NULL,
                     # treeplot = NULL,
                     # upsetplot = NULL,
                     
                     # dotplot_gsea = NULL,
                     # ridgeplot_gsea = NULL,
                     
                     filtered = data.frame(), 
                     #pheatmap = NULL,
                     #volcano0 = NULL,
                     pca_plot = NULL,
                     all_cluster_heatmap = NULL,
                     all_together_tcpm = c(),
                     
                     gsea_list = list(),
                     #gsea = NULL,
                     current_gsea = NULL,
                     #gseaplot = NULL,
                     
                     enriched_kegg_pathways = NULL,
                     #enriched_kegg_pathways_desc = NULL,
                     #enriched_kegg_pathways_id = NULL,
                     
                     input_data = NULL,
                     non_redundant_input_data_size = NULL,
                     input_genes_with_annotations = c(),
                     
                     kegg_ora_vector = NULL,
                     kegg_gsea_vector = NULL,
                     
                     
                     ora = data.frame(), 
                     msi_df = data.frame(), 
                     gsea_input = NULL,
                     gsea = data.frame(),
                     gsea_kegg = data.frame(),
                     go = NULL,
                     enricher_plots = NULL,
                     
                     fc = c(),
                     
                     dim_msi_df = NULL ,
                     
                     kegg_gso_col_low = "red",
                     kegg_gso_col_high = "green",
                     kegg_gsea_col_low =  "red",
                     kegg_gsea_col_high = "green",
                     
                     aws_files = c(),
                     custom_conditions_names = data.frame()
                     
)

print("Initialization complete")