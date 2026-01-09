cur_dir = getwd()
data_folder = "/mnt/storage/alexandros/bcro/data/central_efs/data/"
setwd(data_folder)

annotations_default = as.data.frame(read.delim2("v_union2_v3.2_annotation_table.txt", header = T, sep = "\t", fill = T))


# this file was calculated in the annot_table.r in bcro
txi_names = as.data.frame(read.delim2("all_Names_mapped.txt", header = T, sep = "\t", fill = T,  quote = ""))  %>% 
  dplyr::rename( GENEID = gene_id_internal)
mitochondrial <- subset(txi_names %>% filter(source == "merge_mito_tabular.txt"))$GENEID
chloroplast <- subset(txi_names %>% filter(source == "merge_chloro_tabular.txt"))$GENEID
nuclear <- subset(txi_names %>% filter(!source %in% c("merge_chloro_tabular.txt", "merge_mito_tabular.txt" )))$GENEID
#tx2gene <- txi_names %>% select(Name, GENEID) %>% dplyr::rename(TXNAME = Name)

annotations_default <- annotations_default %>% 
  select(colnames(annotations_default)[!grepl("exp_", colnames(annotations_default))]) %>%
  select(-any_of(c("protein_product", "seqid", "start", "end", "strand"))) 


annotations_columns = sort(colnames(annotations_default))


#####
go_file = "Ka.vunion2_GO_terms_manipulate.txt"
if(!file.exists(go_file)) {
  annotations_default %>%
    filter(Name %in% txi_names$GENEID) %>% 
    select(Name, em_GOs, b2go_GO_id, ips_GO_ids, arabThal_GO_ID, w2go_GO_id, cr2go_GO_id) %>%
    gather(key, val, c( "em_GOs", "b2go_GO_id", "ips_GO_ids", "arabThal_GO_ID", "w2go_GO_id", "cr2go_GO_id")) %>%
    filter(val != "-" & val != "" & !is.na(val)) %>%
    separate_rows(val, sep = ",") %>%
    group_by(key, val) %>%
    summarise(freq = n_distinct(Name)) %>%
    write.table("all_go_terms.txt", sep = "\t", row.names = F, quote = F)

  go2gene_all <- annotations_default %>%
        select(Name, em_GOs, b2go_GO_id, ips_GO_ids, arabThal_GO_ID, w2go_GO_id, cr2go_GO_id) %>%
        gather(key, val, c( "em_GOs", "b2go_GO_id", "ips_GO_ids", "arabThal_GO_ID", "w2go_GO_id", "cr2go_GO_id")) %>%
        filter(val != "-" & val != "" & !is.na(val)) %>%
        select(val, Name) %>%
        separate_rows(val, sep = ",") %>%
        unique
  go2gene_all <- go2gene_all %>% filter(val != "-" & val != "" & !is.na(val) & val != "NA")
  go2gene_all %>% write.table(go_file, sep = "\t", row.names = F, quote = F)
}else{
  go2gene_all <- as.data.frame(read.table(go_file, header = T, sep = "\t"))
}



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
### until here we need

demo_genes <- as.data.frame(read.table("genes_to_test.txt", header = T))





# this file i do in bcro big with my script prepare annot_table
ko_annots <- as.data.frame(read.table("vunion2_KO_terms.txt", header = T, sep = "\t")) 
# here i have updated this file to include also the KO with any KOala threshold , even below 50. check kfk_KO_terms_thresholds_shiny in the prepare annot table. the filter was initially vefore.
# but i need to include all for the plot graph of the explorer.
#So i do the following

ko_annots <- ko_annots %>% 
filter(em == 1 | (koala_threshold >= 50 & koala_Evalue <= 10e-4)) %>%
rename(Name = V1) %>% 
unique


ko_annots <- suppressWarnings(
  ko_annots %>%
    group_by(Name, KO) %>%
    summarise(em = max(em, na.rm=T),
              koala_threshold = max(koala_threshold, na.rm=T)
              ) %>%
    ungroup
)

koala_thresholds <- suppressWarnings(ko_annots %>%
                                        group_by(KO) %>%
                                        summarise(koala_threshold = max(koala_threshold, na.rm = T),
                                                  em = max(em, na.rm = T)) %>%
                                        mutate(koala_threshold = ifelse(is.infinite(koala_threshold), NA, koala_threshold)))


ko_symbols <- as.data.frame(read.delim2("KO_titles_new.txt", header = F, sep = ";", quote = "")) %>%
  mutate(V3 = str_remove(V3, " ")) %>%
  rename(Symbol = V2, Node_Description = V3) 


# this file comes from my MAC i do with the :
# the Descriptions of the KEGG pathways come from the /bootstrap_db/database/Ka_KEGGpathways.txt that is produced by the KEGGmap2jsexport_vUnion2.py
ko_Pathway_terms <- as.data.frame(read.table("Ka_KEGGpathways.txt", header = F, sep = "\t", stringsAsFactors = F, 
                                              colClasses = c("character", "character")))  %>%
  mutate(V2 = sub(" \\[.*", "", V2))
                      
                                              


# this file i do in bcro big with my script prepare_annot_table
ko_Pathway <- as.data.frame(read.table("vunion2_KO_Pathways.txt", header = T, colClasses = c("character", "character"))) %>%
  dplyr::select(val, Name)  %>% 
  mutate(val = str_replace(val ,"ko", "")) %>%
  filter(val %in% ko_Pathway_terms$V1)  


ko_nodes <- as.data.frame(read.table("KEGGmaps_and_nodes.txt", header = F, sep = "\t") ) %>% 
  rename(map = V1, nodes = V2) %>%
  mutate(map = str_replace(map, "ko", ""))


setwd(cur_dir)