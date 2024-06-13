

  source(file = file.path(version, "server_initialize.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  
  html("overview_text", "Initialization complete.")
  shinyjs::show("data_import_box")
  
  
  source(file = file.path(version, "server_button_listeners.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  
  source(file = file.path(version, "server_retrieve_runs_from_aws.R"),
         local = TRUE,
         encoding = "UTF-8")

  #dropdown for selecting text or file for custom genes
  observeEvent(input$choose_custom_gso, {
    if(input$choose_custom_gso == "File"){
      if(is.null(input$gso_file_input$datapath) || input$gso_file_input$datapath == ""  ){
        shinyjs::disable("submit_custom_gene_set")
      }else{
        shinyjs::enable("submit_custom_gene_set")
      }
    }else if(input$choose_custom_gso == "text"){
      if(is.null(input$gso_text_input) || input$gso_text_input == ""){
        shinyjs::disable("submit_custom_gene_set")
      }else{
        shinyjs::enable("submit_custom_gene_set")
      }
    }
  })

  #enabling the submition with text for custom gene sets
  observeEvent(input$gso_text_input, {
    if(input$choose_custom_gso == "Text"){
      if(!is.null(input$gso_text_input) && input$gso_text_input != ""){
        shinyjs::enable("submit_custom_gene_set")
      }else{
        shinyjs::disable("submit_custom_gene_set")
      }
    }
  })
  
  #loading the count data for custom contrast
  # observeEvent(input$data_file$datapath, {
  #   if(input$choose == "file"){
  #     if(!is.null(input$data_file$datapath) && input$data_file$datapath != "" ){
  #       print("file loaded")
  #       shinyjs::enable("submit_norm_and_de")
  #     }else{
  #       shinyjs::disable("submit_norm_and_de")
  #       
  #     }
  #   }
  # })
  
  
  
  # add DE values to the annotations table
  annotations <- eventReactive(
    rv$update_annotations, {
      if(rv$update_annotations == 0){
        data.frame()
      }else{
        
        annotations_default %>% 
          select(all_of(annotations_columns)) %>% 
          full_join(rv$msi_df) %>% 
          mutate(In_DE_data = ifelse(!is.na(logFC), "Yes", "No"))
      }
      
  })
  
  
  
  observeEvent(
    rv$update_msi_df, {
      if(rv$update_msi_df > 0){
        
        rv$msi_df <- rv$msi_df %>% arrange(desc(logFC))
      # annotations <- annotations %>% select(all_of(annotations_columns))
      # annotations <- annotations %>% full_join(rv$msi_df) %>% mutate(In_DE_data = ifelse(!is.na(logFC), "Yes", "No"))
      # print(paste("dim annotations in submit contrast", dim(annotations), collapse = " "))
      # 
      # prepare gsea input
      rv$gsea_input = rv$msi_df$logFC
      rv$gsea_kegg = data.frame()
      names(rv$gsea_input) <- rv$msi_df$Name
      
      # summarise the DE table to reduce redundancy
      print(paste("msi df dim before group.", dim(rv$msi_df)))
      rv$msi_df <- rv$msi_df %>% group_by(across(c(-Name))) %>% summarise(Name = paste(Name, sep = ",", collapse = ",")) %>% ungroup() 
      rv$msi_df <- rv$msi_df %>% select(Name, everything())
      print(paste("msi df dim After group.", dim(rv$msi_df)))
      # add annotations to the DE table
      rv$msi_df = rv$msi_df %>% mutate(primary_name = Name) %>% separate(primary_name, c("primary_name"), ",") %>% select(Name, everything() ) #%>% select(-condition, )
      rv$msi_df <- rv$msi_df %>% left_join(annotations_default, by = c("primary_name" = "Name"))
      
      #rv$msi_df = as.data.frame(read.table(paste0("MSI_Kal_t_18_june22_annotated/", rv$condition, "_all_values.csv"), header = T , sep = "\t", quote = ""))
      #rv$msi_df = rv$msi_df %>% mutate(primary_name = Name) %>% separate(primary_name, c("primary_name"), ",") %>% select(Name, everything() ) #%>% select(-condition, )
      #
      print("OK msi df")
     
      
      default_columns = c("In_DE_data", colnames(rv$msi_df)) #[1:which(colnames(rv$msi_df) == "Ontology_term") - 1]
      #annotations_columns = colnames(rv$msi_df)[which(colnames(rv$msi_df) == "Ontology_term" ) : ( length(colnames(rv$msi_df) ) -1 )] 
      updatePickerInput(session, "overview_table_Default_values", choices = default_columns, selected = c("Name", "In_DE_data", "logFC", "logCPM", "FDR"))
      #rv$msi_df <- rv$msi_df %>% left_join(annotations)
      updatePickerInput(session, "overview_table_terms", choices = annotations_columns[annotations_columns != "Name"])
      
      
      shinyjs::show("run_gsea_box")
      shinyjs::show("data_overview_box")
      shinyjs::show("data_overview_table_box")
      
      print("msi_df and contrast load succesfully!")
      p(id = "overview_text", "Data loaded succesfully. You can head on to the downstream analysis on the left panel.")
      rv$show_gso_mout = T
      rv$show_gsea_mout = T
    }
  })
  
  
  ##### render dynamic menus
  source(file = file.path(version, "server_render_menus.R"),
         local = TRUE,
         encoding = "UTF-8")

  ######## end render dynamic menus
  
  # observeEvent(input$submit_norm_and_de, {
  #   rv$condition  <- input$Condition
  #   rv$norm_and_de_signal =  rv$norm_and_de_signal + 1 
  #   print("observeEvent(input$submit_norm_and_de")
  #   #shinyjs::disable("submit_contrast")
  #   #shinyjs::disable("submit_norm_and_de")
  # })
  
  
  #picker for KEGGscape in GSO
  observeEvent(input$choose_map_gso, {
    rv$map <- subset(ko_Pathway_terms %>% filter(V2 == input$choose_map_gso))$V1 #str_replace(input$choose_map_ora, "ko", "")
    #print(paste0("observeEvent choose_map_ora rv$map-", rv$map,"-", is.null(rv$map),"-", is.na(rv$map), "-",length(rv$map)))
  })
  
  #picker for KEGGscape in ORA
  # observeEvent(input$choose_map_ora, {
  #   rv$map <- subset(ko_Pathway_terms %>% filter(V2 == input$choose_map_ora))$V1 #str_replace(input$choose_map_ora, "ko", "")
  #   #print(paste0("observeEvent choose_map_ora rv$map-", rv$map,"-", is.null(rv$map),"-", is.na(rv$map), "-",length(rv$map)))
  # })
  # 
  # #picker for KEGGscape in GSEA
  # observeEvent(input$choose_map_gsea, {
  #   rv$map_gsea <- subset(ko_Pathway_terms %>% filter(V2 == input$choose_map_gsea))$V1 
  #   #print(paste0("observeEvent choose_map_gsea rv$map_gsea-", rv$map_gsea,"-", is.null(rv$map_gsea),"-", is.na(rv$map_gsea), "-",length(rv$map_gsea)))
  # })
  
  
  observeEvent(input$choose_GO, {
    # rv$okplot_gseaplot <- TRUE
    rv$go_term <- input$choose_GO
    #rv$gseaplot = gseaplot(gsea, by = "all", title = rv$go_term, geneSetID = 1)
    #print(paste("observeEvent choose_GO rv$go_term",rv$okplot_gseaplot , rv$go_term))
  })
  
  
  #when i have defined the input data either by filters or by input
  observeEvent(
    eventExpr = rv$input_data, {
      
      if (rv$condition == ""){
        shinyjs::show("text_configure_ora")
      }else{
        shinyjs::hide("text_configure_ora")
        #shinyjs::disable("run")
        
        html("text_ora_plots", "You should submit your ploting configuration first and Plot your table!!")
        shinyjs::show("text_ora_plots")
        
        # html("text_KEGGscape_ORA", "Processing...")
        # shinyjs::show("text_KEGGscape_ORA")
        html("text_KEGGscape_gso", "Processing...")
        shinyjs::show("text_KEGGscape_gso")
        
        # rv$data = switch(input$choose,
        #                  file = read_file(input$data_file$datapath),
        #                  text = input$input_text
        # ) 
        print(paste("dim annotations() in rv$input_data", dim(annotations()), collapse = " "))
        
        shinyjs::show("summary_gene_set_box")
        # rv$input_genes_with_annotations = annotations() %>%
        #   filter(Name %in% rv$input_data) %>%
        #   filter_at(annotations_columns[annotations_columns != "Name" & annotations_columns != "mcl_cluster"], any_vars(. != "-")) %>%
        #   distinct(Name)
        
        # print(paste0("N genes: ", length(rv$input_data), 
        #              "<br/>",
        #              "N annotated genes: ", dim(rv$input_genes_with_annotations)[1]))
        
        output$summary_gene_set <- renderUI({ 
          HTML(paste0("N genes: ", length(rv$input_data)
                      #"<br/>",
                      #"N annotated genes: ", dim(rv$input_genes_with_annotations)[1]
                      ))
          })
        
        rv$IP = input$getIP
        
        rv$run_counter = rv$run_counter +1
        

        
        rv$map <- ""
        rv$go_term <- ""
        
        
        
        #input_data = str_split(rv$data, "\n")[[1]] 
        #print(paste("demo inut data", length(input_data)))
        
        ### this is the diana process
        #genes_to_GO = go2gene_all %>% select(Name)
        # rv$go <- enrich_go_custom(
        #   input_data,
        #   universe = genes_to_GO[, 1],
        #   genes_to_GO,
        #   qvalue = 0.1,
        #   pvalue = 0.05,
        #   GO_type = "BP"
        # )
        
       
        print(paste("MSI input data", length(rv$input_data)))
        rv$dim_msi_df <- length(rv$input_data)
       
       
        
        ## UPDATE the GSO kegg map picker
        gso_kegg_maps = ko_Pathway %>% 
          filter(Name %in% rv$input_data) %>% 
          distinct(val) %>% 
          left_join(ko_Pathway_terms, by = c("val" = "V1")) %>%
          arrange(V2) %>%
          filter(!is.na(V2))
        updateSelectInput(session, "choose_map_gso", choices = gso_kegg_maps$V2)
        rv$map <- ifelse(length(gso_kegg_maps$V2) > 0 , gso_kegg_maps[1,2], rv$map)
        
        #update the select overlaying data:
        updateSelectInput(session, "kegg_gso_plotting_values", choices = rv$condition, selected = "NULL")
        shinyjs::show("gso_kegg_pathway_dropdown")
        shinyjs::show("gso_keggscape_table_box")
        
        
        # find the ko node id terms
        # print("FIND ko_annots")
        # gso_kegg <- ko_annots %>%
        #   filter(Name %in% rv$input_data) %>%
        #   left_join(rv$msi_df %>% select(Name, logFC) %>% separate_rows(Name, sep = ",") ) %>%
        #   group_by(val) %>%
        #   summarise(logFC = max(logFC, na.rm = T)) %>%
        #   ungroup()
        # 
        # rv$kegg_gso_vector <- gso_kegg$logFC
        # names(rv$kegg_gso_vector) <- gso_kegg$val
        # print(paste("kegg_gso_vector",dim(rv$kegg_gso_vector)[1]))
        
        #rv$okplot <- TRUE
        #rv$okplot_KEGGmap_ora <- TRUE 
        rv$okplot_KEGGmap_gso <- TRUE 
        
        print("and now show!")
        shinyjs::show("box_configure_term_ora")
        shinyjs::show("enrichment_term")
        shinyjs::show("submit_enrichment")
        #shinyjs::hide("text_KEGGscape_ORA")
        shinyjs::hide("text_KEGGscape_gso")
        #shinyjs::show("ora_keggscape_table_box")
        
        rv$show_KeggScape_mout = T
        rv$show_ora_mout = T
      }
    })
  
  
  ##### plotly_functions
  source(file = file.path(version, "server_plotly_functions.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  ######## end plotly_functions
 
  
  #TABLES
  source(file = file.path(version, "server_tables.R"),
         local = TRUE,
         encoding = "UTF-8")


  # ORA plots  
  source(file = file.path(version, "server_ORA_plots.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  
  ### NORM plots
  source(file = file.path(version, "server_NORM_plots.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  
  #GSEA plots
  source(file = file.path(version, "server_GSEA_plots.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  
  #GSSO KEGG_map_
  source(file = file.path(version, "server_KEGG_map_gso.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  #GSEA KEGG_map_
  # source(file = file.path(version, "server_KEGG_map_gsea.R"),
  #        local = TRUE,
  #        encoding = "UTF-8")
  