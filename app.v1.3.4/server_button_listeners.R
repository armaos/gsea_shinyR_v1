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
  
  shinyjs::html("pageHeader", rv$condition)
  
  rv$msi_df = as.data.frame(read.table(paste0("MSI_data/", rv$condition, "_all_values.csv"), header = T , sep = "\t"))
  colnames(rv$msi_df)[1] = "Name"
  rv$msi_df <- rv$msi_df %>% select(colnames(rv$msi_df)[!grepl("X", colnames(rv$msi_df))])
  
  # add DE values to the annotations table
  rv$update_annotations =  rv$update_annotations + 1
  rv$update_msi_df =  rv$update_msi_df + 1
})




#if i define gene set by filters
observeEvent(
  eventExpr = input$submit_thresholds, {
    # this is the reading of MSI data
    fc_filter_low = input$slider_lfc[1] 
    fc_filter_high = input$slider_lfc[2] 
    
    #pval_filter = input$slider_pval[1] 
    fdr_filter = input$slider_fdr[1] 
    
    msi_df <- rv$msi_df %>%  
      filter(abs(logFC) >= fc_filter_low & abs(logFC) <= fc_filter_high & FDR <= 10^(-fdr_filter)) %>%
      select(Name, logFC) %>%
      separate_rows(Name, sep = ",")
    
    rv$input_data <- msi_df$Name
    rv$fc = msi_df$logFC
    names(rv$fc)= rv$input_data
    print(paste("filtered genes load succesfully", length(rv$input_data)[1]))
  })

#if i define gene set by input
observeEvent(
  eventExpr = input$submit_custom_gene_set, {
    
    input_genes = switch(input$choose_custom_gso,
                         File = read_file(input$gso_file_input$datapath),
                         Text = input$gso_text_input
    )
    
    input_genes = str_split(input_genes, "\n")[[1]] 
    msi_df <- rv$msi_df %>%  
      separate_rows(Name, sep = ",") %>%
      filter(Name %in% input_genes) %>%
      select(Name, logFC) 
    
    rv$input_data <- input_genes
    
    rv$fc = msi_df$logFC
    names(rv$fc)= msi_df$Name
    print(paste("custom genes load succesfully", length(rv$input_data)))
  })




###

observeEvent(
  eventExpr = input$run_gsea, {
    if (rv$condition == ""){
      shinyjs::show("text_configure_gsea")
    }else{
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
      
      rv$ridgeplot_gsea <- NULL
      rv$dotplot_gsea <- NULL
      rv$gsea <- data.frame()
      rv$okplot_KEGGmap_gsea <- FALSE
      rv$map_gsea <- ""
      rv$okplot_gseaplot <-  FALSE
      
      #####
      gsea_kegg_v = ko_annots %>% 
        inner_join(rv$msi_df %>% select(Name, logFC) %>% separate_rows(Name,sep = ",")) %>%
        #inner_join(msi_df %>% select(Name, logFC) %>% separate_rows(Name,sep = ",")) %>%
        group_by(KO) %>%
        summarise(logFC = max(logFC, na.rm = T)) %>%
        ungroup()
      rv$kegg_gsea_vector = gsea_kegg_v$logFC
      names(rv$kegg_gsea_vector) <- gsea_kegg_v$KO
      
      print("GSEA")
      # first run the KEGG gsea just to have it stored so that the maps can be loaded
      if(dim(rv$gsea_kegg)[1] == 0){
        rv$gsea_kegg = GSEA(rv$gsea_input, 
                            eps = 1e-40,
                            TERM2GENE=ko_Pathway, 
                            TERM2NAME = ko_Pathway_terms, 
                            )
        
        updateSelectInput(session, "choose_map_gsea", choices = rv$gsea_kegg$Description[!is.na(rv$gsea_kegg$Description)])
        rv$map_gsea = ifelse(length(rv$gsea_kegg[, 1]) > 0 , rv$gsea_kegg[1, 1], rv$map_gsea)
        shinyjs::show("gsea_kegg_pathway_dropdown")
        print(paste("rv$map_gsea ", rv$map_gsea ))
      }
      
      print(paste("length rv$gsea_kegg" , length(rv$gsea_kegg)))
      
      if (input$enrichment_term_gsea %in% c("biological_process", "molecular_function", "cellular_component")){
        term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == input$enrichment_term_gsea & ArabT == "0"))$GO_ID)
        term2name = go_terms %>% filter(GO_namespace == input$enrichment_term_gsea & ArabT == "0")
        
        rv$gsea <- GSEA(rv$gsea_input, 
                        eps = 1e-40,
                        TERM2GENE = term2gene,
                        TERM2NAME = term2name)
        

      }else if (input$enrichment_term_gsea %in% c("biological_process_arth", "molecular_function_arth", "cellular_component_arth")){
        term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term_gsea, "_arth") & ArabT == "1"))$GO_ID)
        term2name = go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term_gsea, "_arth") & ArabT == "1")
        
        rv$gsea <- GSEA(rv$gsea_input, 
                        eps = 1e-40,
                        TERM2GENE = term2gene,
                        TERM2NAME = term2name)
        
      }else if(input$enrichment_term_gsea == "KEGG_Pathways"){
        term2gene = ko_Pathway
        term2name = ko_Pathway_terms
        rv$gsea = rv$gsea_kegg
      }
      
      
      
      
      
      rv$dotplot_gsea = dotplot(rv$gsea)
      print("doing ridgeplot_gsea")
      rv$ridgeplot_gsea = ridgeplot(rv$gsea )
      print("ridgeplot_gsea DONE")
      updateSelectInput(session, "choose_GO", choices = rv$gsea$Description[!is.na(rv$gsea$Description)])
      #updateCheckboxGroupInput(session, "custom_upset_terms", choices = rv$gsea$Description[!is.na(rv$gsea$Description)])
      updatePickerInput(session, "custom_upset_terms", choices = rv$gsea$Description[!is.na(rv$gsea$Description)])
      
      rv$show_gsea_plots_mout = T
      
      shinyjs::show("gsea_plots_tabsetPanel")
      
      shinyjs::show("custom_upset_box")
      
      shinyjs::hide("text_KEGGscape_GSEA")
      shinyjs::show("gsea_kegg_pathway_table")
      
      
      
      
      html("text_gsea_table", "GSEA: DONE")
      print("GSEA: DONE")
      
      shinyjs::show("gsea_table_box")
      rv$okplot_KEGGmap_gsea <- TRUE
      rv$okplot_gseaplot <- TRUE
    }
  })



observeEvent(
  eventExpr = input$run_ora, {
    shinyjs::hide("box_ora_send_to_plot")
    shinyjs::hide("ora_table_box")
    print(" input$run_ora")
    
    #rv$okplot_KEGGmap_ora <- FALSE
    rv$barplot1 <- NULL
    rv$barplot2 <- NULL
    rv$dotplot <- NULL
    rv$emapplot <- NULL
    rv$cnetplot <- NULL
    rv$heatplot <- NULL
    rv$treeplot <- NULL
    rv$upsetplot <- NULL
    
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
    # we should also add the input data to the universe
    universe <- rbind(universe, data.frame(Name = rv$input_data)) %>% distinct(Name)
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
    
    rv$ora <- enricher(rv$input_data, 
                       universe = universe$Name,
                       TERM2GENE = term2gene,
                       TERM2NAME = term2name, 
                       minGSSize = input$ORA_minGSSize,
                       maxGSSize = input$ORA_maxGSSize)
    
    print(paste("ORA with enricher :DONE with n. input data:", length(rv$input_data)))
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
    html("text_ora_plots", "Processing...!")
    
    print(paste("input$submit_plot_ora", dim(as.data.frame(rv$ora))[1]))
    if(length(rv$ora ) > 0){
      #if(!is.null(rv$ora) & rv$ora != ""){
      
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
        print("barplot1")
        rv$barplot1 = barplot(rv$ora, showCategory=rv$ora_desc, main = "Barplot (count)" )
        rv$barplot2 = mutate(rv$ora,
                             qscore = -log(p.adjust, base=10)) %>%
          barplot(x="qscore", showCategory=rv$ora_desc, main ="Barplot (qscore)")
        print("dotplot")
        rv$dotplot = dotplot(rv$ora, showCategory=rv$ora_desc)
        print("pairwise_termsim")
        p_t = pairwise_termsim(rv$ora)
        print("emapplot")
        rv$emapplot = emapplot(p_t, showCategory = rv$ora_desc)
        print("cnetplot")
        rv$cnetplot = cnetplot(rv$ora, categorySize="pvalue", showCategory = rv$ora_desc, foldChange = rv$fc )
        print("heatplot")
        rv$heatplot = heatplot(rv$ora, showCategory = rv$ora_desc, foldChange = rv$fc )
        print("upsetplot")
        rv$upsetplot = upsetplot(rv$ora, showCategory = rv$ora_desc)
        print("treeplot")
        #rv$treeplot = treeplot(p_t)
        print("DONE plots")
        shinyjs::hide("text_ora_plots")
        shinyjs::show("ora_plots_tabsetPanel")
        
        
      }else{
        html("text_ora_plots", "<strong>No enriched terms, try with a different threshold!</strong>")
      }
      
    }
  })
