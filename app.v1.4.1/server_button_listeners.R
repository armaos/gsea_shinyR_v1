# button to fill in the demo genes
observeEvent(input$demo, {
  updateTextAreaInput(session, "gso_text_input",
                      value = paste(head(unique(demo_genes$Name) , 40), collapse = "\n") )
  updateSelectInput(session, "choose_gso", selected = "Custom Set")
  
})



# this returns the filtered table from the tcpm input table
observeEvent(input$submit_custom_data_preprocessing, {
  print("eventReactive(submit_custom_data_preprocessing")
  shinyjs::show("overview_text")
  html("overview_text", "Reading Count Data...")
  
  print("Reading...")
  data0 <- switch(input$choose_custom_tcmp,
                  testdata = read.table(file = 'tcpm.csv', sep = ',', header = TRUE) %>% select(-"Undetermined_S0"),
  )
  
  
  colnames(data0)[colnames(data0)=="gene_name"]="GeneSymbol"
  
  #data<- head(data0,1000)
  data <- data0
  #print(paste("data", head(data, 1)))
  colnames(data)[colnames(data)=="Geneid"]="Geneid"
  
  source(file = file.path(version, "server_data_preprocessing.R"),
         local = TRUE,
         encoding = "UTF-8")
})



# run the dif exp analysis from Meik : this is When i select the Submit DE analysis
observeEvent(input$run_DE, {
  source(file = file.path(version, "server_run_DE_Meik.R"),
         local = TRUE,
         encoding = "UTF-8")

})


# submit contrast
observeEvent(input$submit_contrast, {
  #shinyjs::disable("submit_contrast")
  #shinyjs::disable("submit_norm_and_de")
  #print("submit_norm_and_de")
  shinyjs::show("submitted_contrast1")
  #shinyjs::show("submitted_contrast2")
  rv$condition  <- input$Condition
  print("observeEvent(input$submit_contrast")
  print(paste(rv$condition, collapse = ","))
  cond = rv$condition[1]
  
  rv$msi_df = as.data.frame(read.table(paste0("MSI_data/", cond, "_all_values.csv"), header = T , sep = "\t"))
  rv$msi_df <- rv$msi_df %>% select(-c(genelabels, change))
  colnames(rv$msi_df)[1] = "Name"
  rv$msi_df <- rv$msi_df %>% select(colnames(rv$msi_df)[!grepl("X", colnames(rv$msi_df))])

  for(col in c("logCPM", "logFC", "LR", "PValue", "PValue", "FDR", "genelabels", "change")){
    colnames(rv$msi_df) = rename_msi_df_columns(colnames(rv$msi_df), col, paste(cond, col, sep = "."))
  }
  
  
  if(length(rv$condition) > 1){
    for(cond in rv$condition[2:length(rv$condition)]){
      msi_tmp = as.data.frame(read.table(paste0("MSI_data/", cond, "_all_values.csv"), header = T , sep = "\t"))
      colnames(msi_tmp)[1] = "Name"
      rv$msi_df <- rv$msi_df %>% select(colnames(rv$msi_df)[!grepl("X", colnames(rv$msi_df))])
      for(col in c("logCPM", "logFC", "LR", "PValue", "PValue", "FDR", "genelabels", "change")){
        colnames(msi_tmp) = rename_msi_df_columns(colnames(msi_tmp), col, paste(cond, col, sep = "."))
      }
      rv$msi_df <- rv$msi_df %>%
        full_join(msi_tmp)
    }  
  }
  
  print(paste("MSI df dim: ",dim(rv$msi_df), collapse = ", "))
  

  # add DE values to the annotations table
  rv$update_annotations =  rv$update_annotations + 1
  rv$update_msi_df =  rv$update_msi_df + 1
  rv$update_msi_clustering =  rv$update_msi_clustering + 1
  
})




#if i define gene set by filters
observeEvent(
  eventExpr = input$submit_thresholds, {
    rv$submission_by = input$choose_gso
    print(paste("rv$submission_by", rv$submission_by))
    # this is the reading of MSI data
    fc_filter_low <<- input$slider_lfc[1] 
    fc_filter_high <<- input$slider_lfc[2] 
    
    #pval_filter = input$slider_pval[1] 
    #fdr_filter <<- get_fdr_threshold(input$slider_fdr[1])
    
    print(paste("submit thresholds ", fdr_filter))
    #### server clustering does the heatmap but also calculate the deg (are the Kals that are above the threshold in any of the conditions)
    source(file = file.path(version, "server_clustering.R"),
           local = TRUE,
           encoding = "UTF-8")

    #msi_df <- rv$msi_df %>%  
    
    non_redundant_input_data_size  <- rv$msi_df %>%  
      filter(Name %in% degs) %>%
      select(Name, contains("FDR"), contains("logFC"))
    
    rv$non_redundant_input_data_size = dim(non_redundant_input_data_size)[1]
    rv$input_data <- non_redundant_input_data_size %>%
      separate_rows(Name, sep = ",")
      
    
    # rv$input_data <- msi_df$Name
    # rv$fc = msi_df$logFC
    # names(rv$fc)= rv$input_data
    print(paste("filtered genes load succesfully by thresholds", length(rv$input_data$Name)))
  })


#if i define gene set by groups
observeEvent(
  eventExpr = input$submit_gene_select_by_group, {
    rv$submission_by = input$choose_gso
    print(paste("rv$submission_by", rv$submission_by))
    group = input$gene_select_by_group
    values = input$gene_select_by_group_values
    
    
    #msi_df <- rv$msi_df %>%  
    rv$input_data <- rv$msi_df %>% 
      select(Name, contains("logFC"), contains("FDR"), all_of(group)) %>%
      #select(all_of(c("Name", "logFC", group))) %>%
      filter(!is.na(!!as.symbol(group))) %>% 
      separate_rows(!!as.symbol(group), sep = "[,]") %>%  
      filter(!!as.symbol(group) %in% values ) %>%
      #distinct(Name, logFC) %>%
      unique() %>%
      #distinct(Name, contains("logFC")) %>%
      separate_rows(Name, sep = ",")
    

    
    print(paste("filtered genes load succesfully by group", length(rv$input_data$Name)))
  })


#if i define gene set by input
observeEvent(
  eventExpr = input$submit_custom_gene_set, {
    rv$submission_by = input$choose_gso
    print(paste("rv$submission_by", rv$submission_by))
    
    input_genes = switch(input$choose_custom_gso,
                         File = read_file(input$gso_file_input$datapath),
                         Text = input$gso_text_input
    )
    
    rv$input_genes = str_split(input_genes, "\n")[[1]] 
    
    rv$input_data <- rv$msi_df %>%  
      separate_rows(Name, sep = ",") %>%
      filter(Name %in% input_genes) %>%
      select(Name, contains("logFC"), contains("FDR")) 
    
    # rv$input_data <- input_genes
    # 
    # rv$fc = msi_df$logFC
    # names(rv$fc)= msi_df$Name
    print(paste("custom genes load succesfully by input names", length(rv$input_data$Name)))
  })




###

observeEvent(
  eventExpr = input$run_gsea, {
    if (rv$condition == ""){
      shinyjs::show("text_configure_gsea")
    }else{
      print("GSEA")
      
      shinyjs::hide("text_configure_gsea")
      
      
      html("text_gsea_plot_overview", "Processing...")
      shinyjs::show("text_gsea_plot_overview")
      
      html("text_gsea_plot", "Processing...")
      shinyjs::show("text_gsea_plot")
      
      html("text_gsea_table", "Processing...")
      shinyjs::show("text_gsea_table")
      
      html("text_KEGGscape_GSEA", "Processing...")
      shinyjs::show("text_KEGGscape_GSEA")
      
      ######
      
      # rv$ridgeplot_gsea <- NULL
      # rv$dotplot_gsea <- NULL
      
      #rv$gsea <- data.frame()
      rv$okplot_KEGGmap_gsea <- FALSE
      rv$map_gsea <- ""
      rv$okplot_gseaplot <-  FALSE
      gsea_cond = input$choose_gsea_condition
      rv$gsea_list[[gsea_cond]] = data.frame()
      
      # print("run gsea")
      # print(paste(dim(rv$msi_df), collapse = ", "))
      # print(paste(colnames(rv$msi_df), collapse = ", "))
      # print(paste(head(rv$msi_df,1), collapse = ", "))
      # 
      gsea_input = get_filtered_msi_df(rv$msi_df %>% select(Name, contains("logFC"), contains("FDR")), gsea_cond, 1, 0 ) %>%
        arrange(desc(logFC)) 
        
      gsea_input_vector = gsea_input$logFC
      names(gsea_input_vector) = gsea_input$Name
      print(paste("length gsea_input : ", length(gsea_input_vector)))

      
      
      ##### seems that this is not used at all
      # gsea_kegg_v = ko_annots %>% 
      #   inner_join(get_filtered_msi_df(rv$msi_df, gsea_cond, 1, 0 )  %>% 
      #                select(Name, logFC) %>% 
      #                separate_rows(Name,sep = ",")
      #              ) %>%
      #   #inner_join(msi_df %>% select(Name, logFC) %>% separate_rows(Name,sep = ",")) %>%
      #   group_by(KO) %>%
      #   summarise(logFC = max(logFC, na.rm = T)) %>%
      #   ungroup()
      # rv$kegg_gsea_vector = gsea_kegg_v$logFC
      # names(rv$kegg_gsea_vector) <- gsea_kegg_v$KO
      
     
      
      source(file = file.path(version, "server_create_gsea_tab_dynamic.R"),
             local = TRUE,
             encoding = "UTF-8")
      
      # first run the KEGG gsea just to have it stored so that the maps can be loaded
      # if(dim(rv$gsea_kegg)[1] == 0){
      #   rv$gsea_kegg = GSEA(gsea_input_vector, 
      #                       eps = 1e-40,
      #                       TERM2GENE = ko_Pathway, 
      #                       TERM2NAME = ko_Pathway_terms, 
      #   )
      #   
      #   updateSelectInput(session, "choose_map_gsea", choices = rv$gsea_kegg$Description[!is.na(rv$gsea_kegg$Description)])
      #   rv$map_gsea = ifelse(length(rv$gsea_kegg[, 1]) > 0 , rv$gsea_kegg[1, 1], rv$map_gsea)
      #   shinyjs::show("gsea_kegg_pathway_dropdown")
      #   print(paste("rv$map_gsea ", rv$map_gsea ))
      # }
      #print(paste("length rv$gsea_kegg" , length(rv$gsea_kegg)))
      
      if (input$enrichment_term_gsea %in% c("biological_process", "molecular_function", "cellular_component")){
        term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == input$enrichment_term_gsea & ArabT == "0"))$GO_ID)
        term2name = go_terms %>% filter(GO_namespace == input$enrichment_term_gsea & ArabT == "0")
    
        
      }else if (input$enrichment_term_gsea %in% c("biological_process_arth", "molecular_function_arth", "cellular_component_arth")){
        term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term_gsea, "_arth") & ArabT == "1"))$GO_ID)
        term2name = go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term_gsea, "_arth") & ArabT == "1")
        
        
      }else if(input$enrichment_term_gsea == "KEGG_Pathways"){
        term2gene = ko_Pathway
        term2name = ko_Pathway_terms
       
      }
      
      rv$gsea_list[[gsea_cond]] <- GSEA(gsea_input_vector, 
                                       eps = 1e-40,
                                       TERM2GENE = term2gene, 
                                       TERM2NAME = term2name)
      
      
      
      rv$current_gsea <- rv$gsea_list[[gsea_cond]] 
      
      print("doing dotplot_gsea")
      #rv$dotplot_gsea = dotplot(rv$gsea_list[[gsea_cond]])
      rv[[paste0("dotplot_gsea","_", gsea_cond)]] = dotplot(rv$gsea_list[[gsea_cond]])
      
      print("doing ridgeplot_gsea")
      #rv$ridgeplot_gsea = ridgeplot(rv$gsea_list[[gsea_cond]] )
      rv[[paste0("ridgeplot_gsea","_", gsea_cond)]] = ridgeplot(rv$gsea_list[[gsea_cond]] )
      
      valid_terms = sort(rv$gsea_list[[gsea_cond]]$Description[!is.na(rv$gsea_list[[gsea_cond]]$Description)])
      updateSelectInput(session, paste0("choose_GO_", gsea_cond ), choices = valid_terms)
      #updateCheckboxGroupInput(session, "custom_upset_terms", choices = rv$gsea_list[[gsea_cond]]$Description[!is.na(rv$gsea_list[[gsea_cond]]$Description)])
      #print("before custom_upset_terms_list")
      #print(length(custom_upset_terms_list))
      custom_upset_terms_list <<- sort(unique(c(custom_upset_terms_list, valid_terms)))
      #print("after custom_upset_terms_list")
      #print(length(custom_upset_terms_list))
      updatePickerInput(session, "custom_upset_terms", choices = custom_upset_terms_list)
      
      rv$show_gsea_plots_mout = T
      
      shinyjs::show("gsea_plots_tabsetPanel")
      shinyjs::show("gsea_plots_dynamic_tabsetPanel")
      
      shinyjs::show("custom_upset_box")
      
      shinyjs::hide("text_KEGGscape_GSEA")
      shinyjs::show("gsea_kegg_pathway_table")
      
      
      
      
      html("text_gsea_table", "GSEA: DONE")
      
      shinyjs::show("gsea_table_box")
      shinyjs::show(paste0("choose_GO_box_",gsea_cond ))
      
      #rv$okplot_KEGGmap_gsea <- TRUE
      #rv$okplot_gseaplot <- TRUE
      print("GSEA: DONE")
    }
  })



observeEvent(
  eventExpr = input$run_ora, {
    shinyjs::hide("box_ora_send_to_plot")
    shinyjs::hide("ora_table_box")
    print(" input$run_ora")
    
    
    # rv$barplot1 <- NULL
    # rv$barplot2 <- NULL
    # rv$dotplot <- NULL
    # rv$emapplot <- NULL
    # rv$cnetplot <- NULL
    # rv$heatplot <- NULL
    # rv$treeplot <- NULL
    # rv$upsetplot <- NULL
    
    rv$ora <- data.frame()
    rv$ora_desc <- NULL
    
    html("text_ora_table", "Processing...")
    shinyjs::show("text_ora_table")
    print("ORA with enricher")
    
    universe = switch(input$choose_universe,
                      all_t = annotations() %>%
                        distinct(Name),
                      all_annot_t = annotations() %>%
                        filter_at(annotations_columns[annotations_columns != "Name" & annotations_columns != "mcl_cluster"], any_vars(. != "-")) %>%
                        distinct(Name),
                      all_de = rv$msi_df %>% select(Name) %>% separate_rows(Name, sep = ",") 
                      
    )
    

    ora_input_data <- get_ora_vector()
    
    print(paste("in ORA, lenght of data", dim(ora_input_data)))
    
    # we should also add the input data to the universe
    universe <- rbind(universe, data.frame(Name = ora_input_data$Name)) %>% distinct(Name)
    print(paste("universe:", dim(universe)[1]))
    
    if (input$enrichment_term %in% c("biological_process", "molecular_function", "cellular_component")){
      term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == input$enrichment_term & ArabT == "0"))$GO_ID) 
      term2name = go_terms %>% filter(GO_namespace == input$enrichment_term & ArabT == "0")
      
    }else if (input$enrichment_term %in% c("biological_process_arth", "molecular_function_arth", "cellular_component_arth")){
      term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term, "_arth") & ArabT == "1"))$GO_ID)
      term2name = go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term, "_arth") & ArabT == "1")
      
    }else if(input$enrichment_term == "KEGG_Pathways"){
      term2gene = ko_Pathway
      term2name = ko_Pathway_terms
    }
    
    print(paste("dim(term2gene)", dim(term2gene)))
    print(paste("dim(term2name)", dim(term2name)))
    
    rv$ora <- enricher(ora_input_data$Name, 
                       universe = universe$Name,
                       TERM2GENE = term2gene,
                       TERM2NAME = term2name, 
                       minGSSize = input$ORA_minGSSize,
                       maxGSSize = input$ORA_maxGSSize)
    
    print(paste("ORA with enricher :DONE with n. input data:", length(ora_input_data$Name)))
    if(dim(data.frame(rv$ora))[1] > 0 ){
      shinyjs::show("box_ora_send_to_plot")
      shinyjs::show("ora_table_box")
      html("text_ora_table", "ORA completed!")
    }else{
      html("text_ora_table", paste0("No gene sets have size between ",input$ORA_minGSSize , " and ",  input$ORA_maxGSSize," 500 ..."))
    }
  })


observeEvent(
  eventExpr = input$submit_plot_ora, {
    #html("text_ora_plots", "Processing...!")
    html("text_ora_plots_dynamic", "Processing...!")
    ora_cond = input$choose_ora_condition
    print(paste("input$submit_plot_ora", dim(as.data.frame(rv$ora))[1]))
    
    if(length(rv$ora ) > 0){
      ora_input_data <- get_ora_vector()
      # ora_input_data = get_filtered_msi_df(rv$input_data, ora_cond, fdr_filter, fc_filter_low)
      fc_vector = ora_input_data$logFC
      names(fc_vector) = ora_input_data$Name
      
      source(file = file.path(version, "server_create_ora_tab_dynamic.R"),
             local = TRUE,
             encoding = "UTF-8")
      
      max_to_show = input$ORA_max_elements
      print(max_to_show)
      p.adjust_low = 0
      p.adjust_high = 1
      if(!is.null(input$ORA_p.adjust)){
        p.adjust_low = as.numeric(str_split(input$ORA_p.adjust, ' ')[[1]][1]) 
        p.adjust_low = ifelse(is.numeric(p.adjust_low), p.adjust_low, 0)
        
        p.adjust_high = as.numeric(str_split(input$ORA_p.adjust, ' ')[[1]][3])
        p.adjust_high = ifelse(is.numeric(p.adjust_high), p.adjust_high, 1)
        
      }
      
      print(paste(p.adjust_low , p.adjust_high))
      
      rv$ora_desc = head(
        subset(
          as.data.frame(rv$ora) %>%
            filter(p.adjust >= p.adjust_low & p.adjust <= p.adjust_high )
        )$Description,
        max_to_show )
      
      print(paste("length of ora desc" , length(rv$ora_desc)))
      if(length(rv$ora_desc) > 0 ){
        
        
        print(paste0("barplot1","_", ora_cond))
        # rv$barplot1 = barplot(rv$ora, showCategory=rv$ora_desc, main = "Barplot (count)" )
        rv[[paste0("barplot1","_", ora_cond)]] = barplot(rv$ora, showCategory=rv$ora_desc, main = "Barplot (count)" )
        
        print(paste0("barplot2","_", ora_cond))
        # rv$barplot2 = mutate(rv$ora,
        #                      qscore = -log(p.adjust, base=10)) %>%
        #   barplot(x="qscore", showCategory=rv$ora_desc, main ="Barplot (qscore)")
        rv[[paste0("barplot2","_", ora_cond)]] = mutate(rv$ora,
                                                 qscore = -log(p.adjust, base=10)) %>%
          barplot(x="qscore", showCategory=rv$ora_desc, main ="Barplot (qscore)")
        
        print(paste0("barplot1","_", ora_cond))
        #rv$dotplot = dotplot(rv$ora, showCategory=rv$ora_desc)
        rv[[paste0("barplot1","_", ora_cond)]] = dotplot(rv$ora, showCategory=rv$ora_desc)
        
        print("pairwise_termsim")
        p_t = pairwise_termsim(rv$ora)
        
        print(paste0("emapplot","_", ora_cond))
        #rv$emapplot = emapplot(p_t, showCategory = rv$ora_desc)
        rv[[paste0("emapplot","_", ora_cond)]] = emapplot(p_t, showCategory = rv$ora_desc)
        
        print(paste0("cnetplot","_", ora_cond))
        #rv$cnetplot = cnetplot(rv$ora, categorySize="pvalue", showCategory = rv$ora_desc, foldChange = rv$fc )
        rv[[paste0("cnetplot","_", ora_cond)]] = cnetplot(rv$ora, categorySize="pvalue", showCategory = rv$ora_desc, foldChange = fc_vector )          
        
        print(paste0("heatplot","_", ora_cond))
        #rv$heatplot = heatplot(rv$ora, showCategory = rv$ora_desc, foldChange = rv$fc )
        rv[[paste0("heatplot","_", ora_cond)]] = heatplot(rv$ora, showCategory = rv$ora_desc, foldChange = fc_vector )
        
        print(paste0("upsetplot","_", ora_cond))
        #rv$upsetplot = upsetplot(rv$ora, showCategory = rv$ora_desc)
        rv[[paste0("upsetplot","_", ora_cond)]] = upsetplot(rv$ora, showCategory = rv$ora_desc)
        
        print(paste0("treeplot","_", ora_cond))
        #rv$treeplot = treeplot(p_t)
        rv[[paste0("treeplot","_", ora_cond)]] = treeplot(p_t)
        
        
        print("DONE plots")
        #shinyjs::hide("text_ora_plots")
        #shinyjs::show("ora_plots_tabsetPanel")
        shinyjs::hide("text_ora_plots_dynamic")
        shinyjs::show("ora_plots_dynamic_tabsetPanel")
        
        
      }else{
        html("text_ora_plots", "<strong>No enriched terms, try with a different threshold!</strong>")
      }
      
    }
      
    
  })

observeEvent(
  eventExpr = input$submit_custom_upset_plot, {  
    if (input$custom_upset_terms != "" & length(input$custom_upset_terms) > 0 ){
      gsea_df_list <- bind_rows(lapply(rv$gsea_list, as.data.frame), .id = "column_label") 
      #print(paste(head(gsea_df_list, 2), collapse = " ,") )
      rv$gsea_custom_upsetplot <- gsea_df_list %>%
        select(Description, core_enrichment) %>%
        filter(Description %in% input$custom_upset_terms) %>%
        separate_rows(core_enrichment, sep = "/") %>%
        mutate(Description = ifelse(str_length(Description) > 30, paste0(str_sub(Description,1,30), "..."), Description)) %>%
        group_by(core_enrichment) %>%
        summarise(terms = list(Description)) %>%
        ggplot(aes(x=terms)) +
        geom_bar() +
        scale_x_upset()
    }else{
      rv$gsea_custom_upsetplot = NULL
    }
  }
)




observeEvent(
  eventExpr = input$submit_number_of_clusters,{
    
    annotation_dol <- as.data.frame(cbind(rownames(rv$all_together_tcpm), 
                                          cluster = cutree(rv$all_cluster_heatmap$tree_row, k = input$number_of_clusters))[,"cluster"])
    colnames(annotation_dol) <- "cluster"
    
    cols <- colorRampPalette(brewer.pal(12, "Paired"))
    mycolors <- cols(length(unique(annotation_dol$cluster)))
    names(mycolors) <- unique(annotation_dol$cluster)
    mycolors <- list(time = mycolors)
  
    rv$all_cluster_heatmap <- pheatmap(rv$all_together_tcpm,
                                       fontsize_row = 3,
                                       angle_col=90,
                                       scale="row",
                                       cluster_cols=F,
                                       color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                                       annotation_row = annotation_dol,
                                       annotation_colors = mycolors,
                                       main = paste0("Expression of all DEGs(",length(degs), " with FDR < " ,fdr_filter ," and |logFC| > ",fc_filter_low ," ) for each sample"),
                                       width=8,
                                       height=10,
                                       show_rownames=F)

    
    clusters_heatmap <- rv$msi_df %>% 
      select(group_id , Name) %>%
      separate_rows(Name, sep = ",") %>%
      left_join(annotation_dol %>% 
                  mutate(Name = rownames(annotation_dol))) %>% 
      distinct(group_id, cluster) %>%
      mutate(cluster = ifelse(is.na(cluster), "NO_DEG", cluster))
    
    colnames(clusters_heatmap) <- c("group_id" , "clusters_heatmap") 
    
    rv$msi_df <- rv$msi_df  %>% 
      select(-any_of(c("clusters_heatmap"))) %>%
      left_join(clusters_heatmap) 
    
    # and update the group selection and the values
    update_group_selection()
    if(!is.null(input$gene_select_by_group)){
      if(length(input$gene_select_by_group) > 0 & input$gene_select_by_group == "clusters_heatmap"){
        update_group_selection_values("clusters_heatmap")
      }
    }
    
    
  })