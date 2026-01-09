# button to fill in the demo genes
observeEvent(input$demo, {
  updateTextAreaInput(session, "gso_text_input",
                      value = paste(head(unique(demo_genes$Name) , 40), collapse = "\n") )
  updateSelectInput(session, "choose_gso", selected = "Custom Set")
  
})

#select_deselect all galaxy inputs
observeEvent(input$selectall_galaxy, {
  if(input$selectall_galaxy > 0 ){
    #print(paste("selectall_galaxy", input$selectall_galaxy, rv$galaxy_outputs$file, sep = "-"))
    if(input$selectall_galaxy %% 2 == 1 ){
      updateSelectInput(session, "galaxy_outputs", choices = galaxy_outputs_choices_list, selected = rv$galaxy_outputs$file2)
      rv$aws_files <- rv$galaxy_outputs$file2
    }else{
      updateSelectInput(session, "galaxy_outputs", choices = galaxy_outputs_choices_list, selected = c())
      rv$aws_files <- c()
    }
    
  }
})

#select_deselect all contrasts
observeEvent(input$selectall_contrasts, {
  if(input$selectall_contrasts > 0 ){
    if(input$selectall_contrasts %% 2 == 1 ){
      updateSelectInput(session, "choose_conditions", choices = rv$conditions, selected = rv$conditions)
    }else{
      updateSelectInput(session, "choose_conditions", choices = rv$conditions, selected = c())
    }
  }
})

# when i want to check the metadata
observeEvent(input$submit_metadata, {
  shinyjs::hide("metadata_table_box")
  shinyjs::hide("metadata_error_p")
  shinyjs::disable("submit_metadata")
  shinyjs::disable("submit_galaxy_outputs")
  rv$correct_metadata_input = T
  metadata_warning = ""
  
  print("Start submit_metadata")
  
  rv$metadata_df <- data.frame()
  if(input$metadata_input_type == "text"){
    lines = stri_split_lines1(gsub('\r','' , input$metadata_text))

    rv$metadata_df <- read.table(text = gsub('[" ", ","]', "\t", lines), sep = "\t", fill = NA)
    if(length(colnames(rv$metadata_df)) != 3){
      rv$correct_metadata_input = F
      metadata_warning = "ERROR!: Check the number of columns"
    }else{
      colnames(rv$metadata_df) <- c("names", "subgroup", "treatment")
      rv$metadata_df <- rv$metadata_df %>% 
        mutate(values = paste(subgroup, treatment, sep = ".")) %>%
        select(names, values)
      
      if(length(unique(rv$metadata_df$names)) != length(unique(rv$metadata_df$values))){
        rv$correct_metadata_input = F
        metadata_warning = "ERROR!: Subgroup values should be distinct among Samples of same Treatment"
      }
    }
   
    
  }else if(input$metadata_input_type == "manual"){
    rv$metadata_df <- get_manual_metadata()
    if(length(unique(rv$metadata_df$names)) != length(unique(rv$metadata_df$values))){
      rv$correct_metadata_input = F
      metadata_warning = "ERROR!: Subgroup values should be distinct among Samples of same Treatment"
    }
    
  }else if(input$metadata_input_type == "galaxy"){
    rv$metadata_df <- data.frame()
    
  }else if(input$metadata_input_type == "file"){
    file <- input$metadata_file
    rv$metadata_df <- read.table(text = gsub('[" ", ","]', "\t", readLines(file$datapath)), sep = "\t", header = input$metadata_header)
    if(length(colnames(rv$metadata_df)) != 3){
      metadata_warning = "ERROR!: Check the number of columns"
      rv$correct_metadata_input = F
    }else{
      colnames(rv$metadata_df) <- c("names", "subgroup", "treatment")
      rv$metadata_df <- rv$metadata_df %>% 
        mutate(values = paste(subgroup, treatment, sep = ".")) %>%
        select(names, values)
      if(length(unique(rv$metadata_df$names)) != length(unique(rv$metadata_df$values))){
        rv$correct_metadata_input = F
        metadata_warning = "ERROR!: Subgroup values should be distinct among Samples of same Treatment"
      }
    }
   
    
  }
  
  print("Mid submit_metadata")
  
  if(rv$correct_metadata_input){
    rownames(rv$metadata_df ) = NULL
    
    metadata <- rv$metadata_df  %>% 
      mutate(names = clean_RunAcc(names)) 
    
    rv$metadata_df <- data.frame(names = input$galaxy_outputs) %>%
      mutate(SampleID = clean_RunAcc(names)) %>%
      left_join(metadata, by = c("SampleID" = "names") )
    
    shinyjs::show("metadata_table_box")
    shinyjs::enable("submit_galaxy_outputs")
    shinyjs::enable("submit_metadata")
  }else{
    showNotification(metadata_warning, duration = 5)
    html("metadata_error_message", metadata_warning)
    shinyjs::show("metadata_error_p")
  }
  
  print("END submit_metadata")
  
})

# this returns the filtered table from the tcpm input table
observeEvent(input$submit_custom_data_preprocessing, {
  shinyjs::disable("submit_custom_data_preprocessing")
  print("eventReactive(submit_custom_data_preprocessing")
  #shinyjs::show("overview_text")
  html("overview_text", "Reading Count Data...")
  
  print("Reading...")
  data0 <- switch(input$choose_custom_tcmp,
                  testdata = rename_dol(read.table(file = 'tcpm.csv', sep = ',', header = TRUE) %>% select(-"Undetermined_S0")),
  )
  
  msi_precompiled_runs = c( "B1", "B1", "B1D", "B1D", "B1L", 
                            "B1L", "B2D", "B2D", "B2Dayo", "B2Dayo", "B2L", "B2L",
                            "B3", "B3", "B3D", "B3D", "B3L", "B3L",
                            "B4", "B4", "B4D", "B4D", "B4L", "B4L", 
                            "G1", "G1", "G1D", "G1D", "G1L", "G1L", 
                            "G2", "G2", "G2D", "G2D", "G2L", "G2L", 
                            "G3", "G3", "G3D", "G3D", "G3L", "G3L", 
                            "G4", "G4", "G4D", "G4D", "G4L", "G4L")
  fastqc_runs = msi_precompiled_runs
  source(file = file.path(version, "server_fastqc_reports.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  
  colnames(data0)[colnames(data0)=="gene_name"]="GeneSymbol"
  
  #data<- head(data0,1000)
  data <- data0

  
  colnames(data)[colnames(data)=="Geneid"]="Geneid"
  
  source(file = file.path(version, "server_data_preprocessing.R"),
         local = TRUE,
         encoding = "UTF-8")
  shinyjs::enable("submit_custom_data_preprocessing")
})



# run the dif exp analysis from Meik : this is When i select the Submit DE analysis
observeEvent(input$run_DE, {
  
  shinyjs::disable("run_DE")
  shinyjs::disable("submit_custom_data_preprocessing")
  create_msi_df = T
  rv$msi_df <- data.frame()
  rv$condition <- input$choose_conditions
  
  updateSelectInput(session, "choose_conditions", selected = c())
  
  norm_and_de_tab_list %>%
    walk(~removeTab("volcano_norm_de_dynamic", .x)) %>%
    walk(~removeTab("pca_norm_de_dynamic", .x)) %>%
    walk(~removeTab("heatmap_norm_de_dynamic", .x))
  norm_and_de_tab_list <<- NULL
  volcano_plot_data_list <<- list()
  pca_data_list <<- list()
  heatmap_plot_data_list <<- list()
  
  for(cond_global in rv$condition){
    local({
      cond = cond_global
      
      rv[[paste0("volcano0_",cond)]] <- NULL
      rv[[paste0("pca_",cond)]] <- NULL
      rv[[paste0("pheatmap_",cond)]] <- NULL
      
      norm_and_de_tab_list <<- c(norm_and_de_tab_list, cond)
      
      
      appendTab(inputId = "volcano_norm_de_dynamic",
                tabPanel(cond,
                         withSpinner(plotOutput(paste0("volcano0_",cond), height = "600px"), type = 6 ))
      )
      
      appendTab(inputId = "heatmap_norm_de_dynamic",
                tabPanel(cond,
                         withSpinner(plotOutput(paste0("pheatmap_", cond), height = "600px"), type = 6 ))
      )
      
      appendTab(inputId = "pca_norm_de_dynamic",
                tabPanel(cond,
                         withSpinner(plotOutput(paste0("pca_", cond), height = "600px"), type = 6 ))
      )

      source(file = file.path(version, "server_NORM_plots_dynamic.R"),
             local = TRUE,
             encoding = "UTF-8")
    })
  }
  
  
  shinyjs::hide("run_gsea_box")
  shinyjs::hide("menu_norm_and_de_panel")
  source(file = file.path(version, "server_run_DE_Meik.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  # for(cond_1 in rv$condition){
  #   local({
  #     cond = cond_1
  #     source(file = file.path(version, "server_NORM_plots_dynamic.R"),
  #                 local = TRUE,
  #                 encoding = "UTF-8")
  #   })}
  
  shinyjs::enable("run_DE")
  shinyjs::enable("submit_custom_data_preprocessing")
  

})


# submit contrast
observeEvent(input$submit_contrast, {
  #shinyjs::disable("submit_contrast")
  #shinyjs::disable("submit_norm_and_de")
  
  shinyjs::show("submitted_contrast1")
  #shinyjs::show("submitted_contrast2")
  rv$condition  <- unname(sapply(input$Condition, rename_condition_dol))
  
  cond = input$Condition[1]
  cond_label = rename_condition_dol(cond)
  
  rv$msi_df <- data.frame()
  rv$msi_df = as.data.frame(read.table(paste0(MSI_data_folder, cond, "_all_values.csv"), header = T , sep = "\t"))
  
  rv$msi_df <- rv$msi_df %>% select(-c(genelabels, change))
  rv$msi_df <- rename_dol(rv$msi_df)
  colnames(rv$msi_df)[1] = "Name"
  rv$msi_df <- rv$msi_df %>% select(colnames(rv$msi_df)[!grepl("X", colnames(rv$msi_df))])

  for(col in c("logCPM", "logFC", "LR", "PValue", "PValue", "FDR", "genelabels", "change")){
    colnames(rv$msi_df) = rename_msi_df_columns(colnames(rv$msi_df), col, paste(cond_label, col, sep = "."))
  }

  
  if(length(input$Condition) > 1){
    for(cond in input$Condition[2:length(input$Condition)]){
      
      cond_label = rename_condition_dol(cond)
      
      msi_tmp = as.data.frame(read.table(paste0(MSI_data_folder, cond, "_all_values.csv"), header = T , sep = "\t"))
      print(paste("msi_tmp", cond, dim(msi_tmp)[1]))
      msi_tmp <- rename_dol(msi_tmp)
      colnames(msi_tmp)[1] = "Name"
      rv$msi_df <- rv$msi_df %>% select(colnames(rv$msi_df)[!grepl("X", colnames(rv$msi_df))])
      for(col in c("logCPM", "logFC", "LR", "PValue", "PValue", "FDR", "genelabels", "change")){
        colnames(msi_tmp) = rename_msi_df_columns(colnames(msi_tmp), col, paste(cond_label, col, sep = "."))
      }
      rv$msi_df <- rv$msi_df %>%
        full_join(msi_tmp)
    }  
  }
  

  #make a PCA:
  
  filtered <- rv$msi_df %>% select(-contains(c("logFC","logCPM", "LR", "PValue", "FDR","genelabels", "change")))
  rownames(filtered) = filtered$Name
  filtered$Name = NULL
  filtered <- as.matrix(filtered)
  rv$filtered = filtered
  
  #donor <- factor(gsub("\\.\\w+$","",colnames(filtered),perl=T))
  #group <- factor(gsub("^\\w+\\.","",colnames(filtered),perl=T))
  subgroup <- factor(gsub("\\.\\w+$","",colnames(filtered),perl=T))  # this is the subgroup
  subgroup_class <- factor(gsub("\\w\\.\\w+$","",colnames(filtered),perl=T))  # this is the subgroup without numbering
  treatment <- factor(gsub("^\\w+\\.","",colnames(filtered),perl=T))  # this is the treatment
  
  
   #### add the fastqc files for the pre compiled)
  msi_precompiled_runs = c( "B1", "B1", "B1D", "B1D", "B1L", 
                            "B1L", "B2D", "B2D", "B2Dayo", "B2Dayo", "B2L", "B2L",
                            "B3", "B3", "B3D", "B3D", "B3L", "B3L",
                            "B4", "B4", "B4D", "B4D", "B4L", "B4L", 
                            "G1", "G1", "G1D", "G1D", "G1L", "G1L", 
                            "G2", "G2", "G2D", "G2D", "G2L", "G2L", 
                            "G3", "G3", "G3D", "G3D", "G3L", "G3L", 
                            "G4", "G4", "G4D", "G4D", "G4L", "G4L")
  fastqc_runs = msi_precompiled_runs
  source(file = file.path(version, "server_fastqc_reports.R"),
         local = TRUE,
         encoding = "UTF-8")
  #### end fastqc
 
  print("i am just putting the treatment here for now this should be customized to the input contrast")
  group = treatment
  create_msi_df = F
  
  y <- calculate_y(filtered, group)
  source(file = file.path(version, "server_PCA_plot.R"),
         local = TRUE,
         encoding = "UTF-8")
  rv$pca_plot = pca_plot
  
  for(tab_to_hide in c("Volcano plot","Heat-Map plot (overview)",  "PCA plot (pairwise)", "Volcano plot (pairwise)",  "Heat-Map plot (pairwise)")){
    hideTab(inputId = "norm_de_plot_tabsetpanel", target = tab_to_hide)
  }
 
  
  rv$show_norm_and_de_mout = T  
  shinyjs::show("menu_norm_and_de_panel")
  # add DE values to the annotations table
  rv$update_annotations =  rv$update_annotations + 1
  rv$update_msi_df =  rv$update_msi_df + 1
  rv$update_msi_clustering =  rv$update_msi_clustering + 1
  
})




#if i define gene set by filters
observeEvent(
  eventExpr = input$submit_thresholds, {
    rv$submission_by = input$choose_gso
    protein_type <<- input$choose_protein_type
    print(paste("rv$submission_by", rv$submission_by))
    
    fc_filter_low <<- input$slider_lfc[1] 
    fc_filter_high <<- input$slider_lfc[2] 
    
    #pval_filter = input$slider_pval[1] 
    #fdr_filter <<- get_fdr_threshold(input$slider_fdr[1])
    
    print(paste("submit thresholds ", fdr_filter))
    
    #### server clustering does the heatmap but also calculate the deg (are the Kals that are above the threshold in any of the conditions)
    # degs is difened in that script too
    source(file = file.path(version, "server_clustering.R"),
           local = TRUE,
           encoding = "UTF-8")
    
    rv$input_data <- data.frame(primary_name = degs)
    
    # non_redundant_input_data_size  <- rv$msi_df %>%  
    #   filter(Name %in% degs) %>%
    #   select(Name, contains("FDR"), contains("logFC"))
    # 
    # rv$non_redundant_input_data_size = dim(non_redundant_input_data_size)[1]
    # rv$input_data <- non_redundant_input_data_size %>%
    #   separate_rows(Name, sep = ",")
      
    shinyjs::show("selected_data_overview_box")
    shinyjs::show("selected_data_overview_table_box")
    print(paste("filtered genes load succesfully by thresholds", length(rv$input_data$primary_name)))
  })


#if i define gene set by groups
observeEvent(
  eventExpr = input$submit_gene_select_by_group, {
    rv$submission_by = input$choose_gso
    print(paste("rv$submission_by", rv$submission_by))
    
    rv$input_data <- get_ora_vector(ora_cond = "") %>% distinct(primary_name)
    
    shinyjs::show("selected_data_overview_box")
    shinyjs::show("selected_data_overview_table_box")
    print(paste("filtered genes load succesfully by group", length(rv$input_data$primary_name)))
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
    
    rv$input_data <- get_ora_vector(ora_cond = "") %>% distinct(primary_name)
    
    # rv$input_data <- rv$msi_df %>%  
    #   separate_rows(Name, sep = ",") %>%
    #   filter(Name %in% input_genes) %>%
    #   select(Name, contains("logFC"), contains("FDR")) 

    shinyjs::show("selected_data_overview_box")
    shinyjs::show("selected_data_overview_table_box")
    print(paste("custom genes load succesfully by input names", length(rv$input_data$primary_name)))
  })



####
# Downloadable csv of selected dataset ----
output$download_selected_Data <- downloadHandler(
  filename = function() {
    "selection_data.csv"
  },
  content = function(file) {
    write.csv(selected_data_table(), file, row.names = FALSE)
  }
)

###

observeEvent(
  eventExpr = input$run_gsea, {
    if (rv$condition == ""){
      shinyjs::show("text_configure_gsea")
    }else{
      withProgress(message = 'GSEA', value = 0.0, {
        gsea_cond = input$choose_gsea_condition
        shinyjs::disable("run_gsea")
        shinyjs::hide("text_configure_gsea")
        shinyjs::hide("gsea_table_box")
        shinyjs::hide(paste0("choose_GO_box_",gsea_cond ))
        
        html("text_gsea_plot_overview", "Processing...")
        shinyjs::show("text_gsea_plot_overview")
        
        html("text_gsea_plot", "Processing...")
        shinyjs::show("text_gsea_plot")
        
        html("text_gsea_table", "Processing...")
        shinyjs::show("text_gsea_table")
        
        html("text_KEGGscape_GSEA", "Processing...")
        shinyjs::show("text_KEGGscape_GSEA")
        
        ######
        incProgress(0.1)
        # rv$ridgeplot_gsea <- NULL
        # rv$dotplot_gsea <- NULL
        
        #rv$gsea <- data.frame()
        rv$okplot_KEGGmap_gsea <- FALSE
        rv$map_gsea <- ""
        rv$okplot_gseaplot <-  FALSE
        
        rv$gsea_list[[gsea_cond]] = data.frame()
        updateSelectInput(session, paste0("choose_GO_", gsea_cond ), choices = c())
        
        # print("run gsea")
        # print(paste(dim(rv$msi_df), collapse = ", "))
        # print(paste(colnames(rv$msi_df), collapse = ", "))
        # print(paste(head(rv$msi_df,1), collapse = ", "))
        # 
        incProgress(0.1, detail = paste("Preparing GSEA input"))
        
        gsea_input = get_filtered_msi_df(rv$msi_df %>% select(primary_name, contains("logFC"), contains("FDR")), gsea_cond, 1, 0 ) %>%
          arrange(desc(logFC)) 
          
        gsea_input_vector = gsea_input$logFC
        names(gsea_input_vector) = gsea_input$primary_name
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
        
        incProgress(0.1, detail = paste("Preparing GSEA output tabs"))
        
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
        
        incProgress(0.2, detail = paste("Preparing GSEA annotation lists"))
        
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
        
        incProgress(0.1, detail = paste("Running GSEA ..."))
        rv$gsea_list[[gsea_cond]] <- GSEA(gsea_input_vector, 
                                         eps = 1e-40,
                                         TERM2GENE = term2gene, 
                                         TERM2NAME = term2name)
        
        
        
        rv$current_gsea <- rv$gsea_list[[gsea_cond]] 
        
        incProgress(0.2, detail = paste("Preparing GSEA plots..."))
        
        print("doing dotplot_gsea")
        #rv$dotplot_gsea = dotplot(rv$gsea_list[[gsea_cond]])
        rv[[paste0("dotplot_gsea","_", gsea_cond)]] = dotplot(rv$gsea_list[[gsea_cond]])
        
        print("doing ridgeplot_gsea")
        #rv$ridgeplot_gsea = ridgeplot(rv$gsea_list[[gsea_cond]] )
        rv[[paste0("ridgeplot_gsea","_", gsea_cond)]] = ridgeplot(rv$gsea_list[[gsea_cond]] )
        
        incProgress(0.1, detail = paste("Preparing GSEA plots..."))
        
        valid_terms = sort(rv$gsea_list[[gsea_cond]]$Description[!is.na(rv$gsea_list[[gsea_cond]]$Description)])
        updateSelectInput(session, paste0("choose_GO_", gsea_cond ), choices = valid_terms)
        print("GSEA vlaid terms")
        print(paste(head(valid_terms), collapse = ","))
        
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
      })
      shinyjs::enable("run_gsea")
      
    }
  })



observeEvent(
  eventExpr = input$run_ora, {
    withProgress(message = 'Running ORA', value = 0, {
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
      
      incProgress(0.2, detail = "Constructing universe...")
      
      universe = switch(input$choose_universe,
                        all_t = annotations() %>%
                          distinct(primary_name),
                        all_annot_t = annotations() %>%
                          filter_at(annotations_columns[annotations_columns != "Name" & annotations_columns != "mcl_cluster"], any_vars(. != "-")) %>%
                          distinct(primary_name),
                        all_de = rv$msi_df %>% select(primary_name) %>% separate_rows(primary_name, sep = ",") 
                        
      )
      
      ora_cond = input$choose_ora_condition
      ora_input_data <- get_ora_vector(ora_cond)
      
      print(paste("in ORA, lenght of data", dim(ora_input_data)))
      
      # we should also add the input data to the universe
      universe <- rbind(universe, data.frame(primary_name = ora_input_data$primary_name)) %>% distinct(primary_name)
      
      html("text_ora_table", paste("Processing...", 
                                   "\n",
                                   "Set size:",dim(ora_input_data)[1],
                                   "\n",
                                   "Universe size:", dim(universe)[1])
           )
      
      incProgress(0.2, detail = "Finding annotations...")
      
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
      
    
      
      incProgress(0.2, detail = "Running ORA...")
      
      rv$ora <- enricher(ora_input_data$primary_name, 
                         universe = universe$primary_name,
                         TERM2GENE = term2gene,
                         TERM2NAME = term2name, 
                         minGSSize = input$ORA_minGSSize,
                         maxGSSize = input$ORA_maxGSSize)
      
      
      print(paste("ORA with enricher dim", dim(rv$ora)[1]))
      
      if(dim(rv$ora)[1] > 0 ){
        shinyjs::show("box_ora_send_to_plot")
        shinyjs::show("ora_table_box")
        html("text_ora_table", paste(
          "Set size:",dim(ora_input_data)[1],
          "<br/>",
          "Universe size", dim(universe)[1],
          "<br/>",
          "ORA completed!")
          )
        
      }else{
        
        html("text_ora_table", paste("No enrichment!",
                                     "<br/>",
                                     "Try again by setting a different gene set or changing ORA configuration.",
                                     "<br/>",
                                     "Set size:", length(ora_input_data$primary_name),
                                     "<br/>",
                                     "Universe size", dim(universe)[1],
                                     "<br/>"
                                     ))
        
      }
      setProgress(0.9, detail = "Finishing process...")
      print(paste("ORA button listener DONE"))
    })
  })


observeEvent(
  eventExpr = input$submit_plot_ora, {
    withProgress(message = 'Plotting ORA', value = 0, {
      html("text_ora_plots_dynamic", "Processing...!")
      
      html("text_tabs_plots_ora_dynamic", 
           paste("ORA:", rv$submission_by, ifelse(rv$submission_by == "By group", paste0(input$gene_select_by_group,": ", input$gene_select_by_group_values), 
                                                  ifelse(rv$submission_by == "Upon filtering", paste("logFC:", fc_filter_low, "FDR:", fdr_filter),
                                                         ifelse(rv$submission_by == "Custom Set","", "")))))
      
      ora_cond = input$choose_ora_condition
  
      if(dim(rv$ora )[1] > 0){
        ora_input_data <- get_ora_vector(ora_cond)
        # ora_input_data = get_filtered_msi_df(rv$input_data, ora_cond, fdr_filter, fc_filter_low)
        fc_vector = ora_input_data$logFC
        names(fc_vector) = ora_input_data$primary_name
        
        setProgress(0.1, detail = "Initiating Dynamic plots...")
        
        source(file = file.path(version, "server_create_ora_tab_dynamic.R"),
               local = TRUE,
               encoding = "UTF-8")
        
        
        #print(max_to_show)
        p.adjust_low = 0
        p.adjust_high = 1
        if(!is.null(input$ORA_p.adjust)){
          p.adjust_low = as.numeric(str_split(input$ORA_p.adjust, ' ')[[1]][1]) 
          p.adjust_low = ifelse(is.numeric(p.adjust_low), p.adjust_low, 0)
          
          p.adjust_high = as.numeric(str_split(input$ORA_p.adjust, ' ')[[1]][3])
          p.adjust_high = ifelse(is.numeric(p.adjust_high), p.adjust_high, 1)
          
        }
        
        setProgress(0.2, detail = "Collecting material...")
        
        print(paste(p.adjust_low , p.adjust_high))
        
        if(dim(as.data.frame(rv$ora))[1] > 0 ){
          
          #get the descriptors
          rv$ora_desc_1 = head(
            subset(
              as.data.frame(rv$ora) %>%
                filter(p.adjust >= p.adjust_low & p.adjust <= p.adjust_high )
            )$Description,
            input$ORA_max_elements_1 )
          
          rv$ora_desc_2 = head(
            subset(
              as.data.frame(rv$ora) %>%
                filter(p.adjust >= p.adjust_low & p.adjust <= p.adjust_high )
            )$Description,
            input$ORA_max_elements_2 )
          
          setProgress(0.3, detail = "Barplot...")
          
          print(paste0("barplot1","_", ora_cond))
          rv[[paste0("barplot1","_", ora_cond)]] = barplot(rv$ora, showCategory=rv$ora_desc_1, main = "Barplot (count)" )
          
          print(paste0("barplot2","_", ora_cond))
          rv[[paste0("barplot2","_", ora_cond)]] = mutate(rv$ora,
                                                   qscore = -log(p.adjust, base=10)) %>%
            barplot(x="qscore", showCategory=rv$ora_desc_1, main ="Barplot (qscore)")
          
          setProgress(0.4, detail = "Dotplot...")
          
          print(paste0("dotplot","_", ora_cond))
          rv[[paste0("dotplot","_", ora_cond)]] = dotplot(rv$ora, showCategory=rv$ora_desc_1)
          
          print("pairwise_termsim")
          p_t = pairwise_termsim(rv$ora)
          
          setProgress(0.4, detail = "Emaplot...")
          
          print(paste0("emapplot","_", ora_cond))
          rv[[paste0("emapplot","_", ora_cond)]] = emapplot(p_t, showCategory = rv$ora_desc_2)
          
          setProgress(0.5, detail = "Cneplot...")
          print(paste0("cnetplot","_", ora_cond))
          rv[[paste0("cnetplot","_", ora_cond)]] = cnetplot(rv$ora, categorySize="pvalue", showCategory = rv$ora_desc_2, foldChange = fc_vector )          
          
          setProgress(0.6, detail = "Heatmap-plot...")
          print(paste0("heatplot","_", ora_cond))
          rv[[paste0("heatplot","_", ora_cond)]] = heatplot(rv$ora, showCategory = rv$ora_desc_1, foldChange = fc_vector )
          
          setProgress(0.7, detail = "UpSetplot...")
          print(paste0("upsetplot","_", ora_cond))
          rv[[paste0("upsetplot","_", ora_cond)]] = upsetplot(rv$ora, showCategory = rv$ora_desc_1)
          
          setProgress(0.8, detail = "Treeplot...")
          print(paste0("treeplot","_", ora_cond))
          
          rv[[paste0("treeplot","_", ora_cond)]] <- tryCatch({
            treeplot(p_t, showCategory = rv$ora_desc_1)
            },
            error=function(cond) {
              hideTab(inputId = paste0("ora_plot_tabsetpanel_",ora_cond), target = "Tree Plot")
              
              message("treeplot error")
              # Choose a return value in case of error
              return(NULL)
            }
          )
          
          
          setProgress(0.9, detail = "Finalizing plots..")
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
      select(group_id , primary_name) %>%
      #separate_rows(Name, sep = ",") %>%
      left_join(annotation_dol %>% 
                  mutate(primary_name = rownames(annotation_dol))) %>% 
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