  
  source(file = file.path(version, "utils.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  source(file = file.path(version, "server_initialize.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  # source(file = file.path(version, "server_load_files.R"),
  #        local = TRUE,
  #        encoding = "UTF-8")
  
  
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

  
  
  set_fdr_threshold <- function(input_fdr, input_by5 ){
    fdr_filter <<- ifelse(input_by5 == T,  10^(-input_fdr) * 5,  10^(-input_fdr) )
    #fdr_filter <<- 0.05
    output$fdr_value <- renderText({ fdr_filter })
  }
  
  
  observeEvent(input$slider_fdr, {
    set_fdr_threshold(input$slider_fdr, input$fdr_by5)
  })
  
  observeEvent(input$fdr_by5, {
    set_fdr_threshold(input$slider_fdr, input$fdr_by5)
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
  
  # observeEvent(input$choose_protein_type, {
  #   if(is.null(input$choose_protein_type) | length(input$choose_protein_type) == 0){
  #     showNotification(paste("Nuclear Genes automatically selected"), duration = 3, type = "warning")
  #   }else if("Chloroplast" %in% input$choose_protein_type | "Mitochondrial" %in% input$choose_protein_type){
  #     showNotification(paste("Chloroplast and Mitochondrial Genes haven't been defined yet. Nuclear Genes automatically selected"), duration = 3, type = "warning")
  #   }
  # })
  
  
  
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
          mutate(In_DE_data = !(if_all(contains("logFC"), is.na))) %>%
          mutate(In_DE_data = ifelse(In_DE_data, "Yes", "No")) %>%
          arrange(desc(In_DE_data)) 
      }
      
  })
  
  
  selected_data_table <- reactive({
    colnames_to_show = c(input$selected_overview_table_Default_values,  
                         input$selected_annotatations_ips,
                         input$selected_annotatations_em,
                         input$selected_annotatations_kfk,
                         input$selected_annotatations_clusters,
                         input$selected_annotatations_b2g,
                         input$selected_annotatations_others
    )
    print(paste("dim annotations() - render data_overview_table"))
    print(paste(dim(annotations()), collapse = " "))
    #print(paste(colnames(annotations()), collapse = " "))
    #print(paste(colnames_to_show, collapse = ", "))
    
    annotations() %>% 
      filter(primary_name %in% rv$input_data$primary_name) %>%
      rename(Name_aliases = Name) %>%
      select(colnames_to_show) 
  })
  

  observeEvent(
    rv$update_msi_df, {
      if(rv$update_msi_df > 0){
        print("start update_msi_df observer")
        #rv$msi_df <- rv$msi_df %>% arrange(desc(logFC))
      # annotations <- annotations %>% select(all_of(annotations_columns))
      # annotations <- annotations %>% full_join(rv$msi_df) %>% mutate(In_DE_data = ifelse(!is.na(logFC), "Yes", "No"))
      # print(paste("dim annotations in submit contrast", dim(annotations), collapse = " "))
      # 
      # prepare gsea input
      
      # rv$gsea_input = rv$msi_df %>% select(c(Name, contains("logFC") ))
      #names(rv$gsea_input) <- rv$msi_df$Name
      
      rv$gsea_kegg = data.frame()
      
      # summarise the DE table to reduce redundancy
      print(paste("msi df dim before Redundancy reduction", dim(rv$msi_df)))
      rv$msi_df <- rv$msi_df %>% 
        group_by(across(c(-Name))) %>% 
        summarise(Name = paste(Name, sep = ",", collapse = ",")) %>%
        ungroup() %>%
        mutate(group_id = 1:n()) 
      
      rv$msi_df <- rv$msi_df %>% select(Name, everything())
      print(paste("msi df dim After Redundancy reduction", dim(rv$msi_df)))
      # add annotations to the DE table
      rv$msi_df = rv$msi_df %>% 
        mutate(primary_name = Name) %>% 
        separate(primary_name, c("primary_name"), ",") %>% 
        select(primary_name, Name, everything() ) #%>% select(-condition, )
      
      default_columns = c("In_DE_data", colnames(rv$msi_df)) #[1:which(colnames(rv$msi_df) == "Ontology_term") - 1]
      default_columns = default_columns[default_columns != "Name"]
      default_columns = c(default_columns, "Name_aliases")
      rv$msi_df <- rv$msi_df %>% left_join(annotations_default, by = c("primary_name" = "Name")) 
      
      #rv$msi_df = as.data.frame(read.table(paste0("MSI_Kal_t_18_june22_annotated/", rv$condition, "_all_values.csv"), header = T , sep = "\t", quote = ""))
      #rv$msi_df = rv$msi_df %>% mutate(primary_name = Name) %>% separate(primary_name, c("primary_name"), ",") %>% select(Name, everything() ) #%>% select(-condition, )
      #
      print("OK msi df")
     
      
      ##### UPDATE COLUMNS IN THE DATA overview
      #annotations_columns = colnames(rv$msi_df)[which(colnames(rv$msi_df) == "Ontology_term" ) : ( length(colnames(rv$msi_df) ) -1 )] 
      contrast_cols = colnames(head(rv$msi_df, 1) %>% select(contains(c("logFC", "logCPM" ,"FDR"))))
      updatePickerInput(session, "overview_table_Default_values", choices = default_columns, selected = c("Name_aliases", "primary_name", "In_DE_data", contrast_cols))
      #rv$msi_df <- rv$msi_df %>% left_join(annotations)
      updatePickerInput(session, "annotatations_ips", choices = annotations_columns[startsWith(annotations_columns, prefix = "ips_")])
      updatePickerInput(session, "annotatations_em", choices = annotations_columns[startsWith(annotations_columns, prefix = "em_")])
      updatePickerInput(session, "annotatations_kfk", choices = annotations_columns[startsWith(annotations_columns, prefix = "koala_")])
      updatePickerInput(session, "annotatations_clusters", choices = annotations_columns[annotations_columns=="mcl_cluster" | startsWith(annotations_columns, prefix = "cluster_")])
      updatePickerInput(session, "annotatations_b2g", choices = annotations_columns[startsWith(annotations_columns, prefix = "b2go_")])
      updatePickerInput(session, "annotatations_dbCAN2", choices = annotations_columns[startsWith(annotations_columns, prefix = "dbCAN2_")])
      
      updatePickerInput(session, "annotatations_others", choices = annotations_columns[!startsWith(annotations_columns, prefix = "ips_") & 
                                                                                         !startsWith(annotations_columns, prefix = "em_") &
                                                                                         !startsWith(annotations_columns, prefix = "koala_") &
                                                                                         !startsWith(annotations_columns, prefix = "cluster_") &
                                                                                         !startsWith(annotations_columns, prefix = "b2go_") &
                                                                                         !startsWith(annotations_columns, prefix = "dbCAN2_") &
                                                                                         annotations_columns != "mcl_cluster" &
                                                                                         annotations_columns != "Name" & 
                                                                                         annotations_columns != "primary_name" 
                                                                                       
                                                                                         ])
      ##### UPDATE COLUMNS IN THE SELECTED DATA overview
      #annotations_columns = colnames(rv$msi_df)[which(colnames(rv$msi_df) == "Ontology_term" ) : ( length(colnames(rv$msi_df) ) -1 )] 
      contrast_cols = colnames(head(rv$msi_df, 1) %>% select(contains(c("logFC", "logCPM" ,"FDR"))))
      updatePickerInput(session, "selected_overview_table_Default_values", choices = default_columns, selected = c("Name_aliases", "primary_name", "In_DE_data", contrast_cols))
      #rv$msi_df <- rv$msi_df %>% left_join(annotations)
      updatePickerInput(session, "selected_annotatations_ips", choices = annotations_columns[startsWith(annotations_columns, prefix = "ips_")])
      updatePickerInput(session, "selected_annotatations_em", choices = annotations_columns[startsWith(annotations_columns, prefix = "em_")])
      updatePickerInput(session, "selected_annotatations_kfk", choices = annotations_columns[startsWith(annotations_columns, prefix = "koala_")])
      updatePickerInput(session, "selected_annotatations_clusters", choices = annotations_columns[annotations_columns=="mcl_cluster" | startsWith(annotations_columns, prefix = "cluster_")])
      updatePickerInput(session, "selected_annotatations_b2g", choices = annotations_columns[startsWith(annotations_columns, prefix = "b2go_")])
      updatePickerInput(session, "selected_annotatations_dbCAN2", choices = annotations_columns[startsWith(annotations_columns, prefix = "dbCAN2_")])
      
      updatePickerInput(session, "selected_annotatations_others", choices = annotations_columns[!startsWith(annotations_columns, prefix = "ips_") & 
                                                                                                  !startsWith(annotations_columns, prefix = "em_") &
                                                                                                  !startsWith(annotations_columns, prefix = "koala_") &
                                                                                                  !startsWith(annotations_columns, prefix = "cluster_") &
                                                                                                  !startsWith(annotations_columns, prefix = "b2go_") &
                                                                                                  !startsWith(annotations_columns, prefix = "dbCAN2_") &
                                                                                                  annotations_columns != "mcl_cluster" &
                                                                                                  annotations_columns != "Name" & 
                                                                                                  annotations_columns != "primary_name" 
                                                                                         
                                                                                           
                                                                                           
      ])
      
      
      
      
      shinyjs::show("run_gsea_box")
      shinyjs::show("data_overview_box")
      shinyjs::show("data_overview_table_box")
      
      print("msi_df and contrast load succesfully!")
      p(id = "overview_text", "Data loaded succesfully. You can head on to the downstream analysis on the left panel.")
      rv$show_gso_mout = T
      rv$show_gsea_mout = T
      
      print("start update_msi_df observer: OK")
    }
  })

 
  
  #
  
  
  
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
  
 
  
  # observeEvent(input$choose_GO, {
  #  
  #   if( input$choose_GO != "" && length(input$choose_GO) > 0 && !is.na(input$choose_GO)){  
  #    
  #     shinyjs::hide("text_gsea_plot")
  #    
  #     
  #     print("input$choose_GO   HERE I SHOULD THINK OF SOMETHING")
  #     
  #     rv[[paste0("gseaplot","_", input$choose_gsea_condition)]] = gseaplot2(rv$gsea, geneSetID = which(rv$gsea$Description == input$choose_GO) , title = input$choose_GO)
  #     
  #     
  #   }
  #   
  #  
  # })
  
  
  ########
  # When i select gene set by group. This will load the column sin the picker to choose from
  update_group_selection <- function(){
    choices_by_group_selection = c(annotations_columns[annotations_columns=="mcl_cluster" | startsWith(annotations_columns, prefix = "cluster_")],
                                   annotations_columns[endsWith(annotations_columns, suffix = "_id")],
                                   annotations_columns[endsWith(annotations_columns, suffix = "_OGs")],
                                   annotations_columns[endsWith(annotations_columns, suffix = "_KEGG_Pathway")]
    )
    if("clusters_heatmap" %in% colnames(rv$msi_df)){
      choices_by_group_selection <- c("clusters_heatmap", choices_by_group_selection)
    } 
    
    updatePickerInput(session, "gene_select_by_group", choices = choices_by_group_selection, selected = "")
  }
  #
  observeEvent(
    eventExpr = input$choose_gso, {
      if(input$choose_gso == 'By group'){
        update_group_selection()
      }
    })
  
  # When i select gene set by group and pick the groups/annotations i want, this will load up the picker with the suitable values
  update_group_selection_values <- function(select_group){
    print(paste("in update_group_selection_values ", select_group))
    values = rv$msi_df %>% 
      select(all_of(select_group)) %>% 
      filter(!is.na(!!as.symbol(select_group)) &  !!as.symbol(select_group) != "-") %>%
      separate_rows(!!as.symbol(select_group), sep = "[,]")
    print(paste(dim(values)))
    updatePickerInput(session, "gene_select_by_group_values", choices = unique(values[,]), selected = NULL)
    print(paste("in update_group_selection_values ", "OK"))
  }
  
  observeEvent(
    eventExpr = input$gene_select_by_group, {
      if(input$choose_gso == 'By group'){
        if(input$gene_select_by_group != ""){
          update_group_selection_values(input$gene_select_by_group)  
        }
      }
    })
  
  #
  observeEvent(
    eventExpr = input$gene_select_by_group_values, {
      if(length(input$gene_select_by_group_values) > 0 ){
        shinyjs::enable("submit_gene_select_by_group")
      }else{
        shinyjs::disable("submit_gene_select_by_group")
      }
      
    })
  ########
  observeEvent(
    eventExpr = rv$condition, {
      if(length(rv$condition) > 0  & rv$condition != "" ){
        #shinyjs::html("pageHeader", paste(collapse = ",", rv$condition))
        update_notifications(rv$condition)
        updateSelectInput(session, "choose_ora_condition", choices = rv$condition)
        updateSelectInput(session, "choose_gsea_condition", choices = rv$condition)
        source(file = file.path(version, "server_create_dynamic_tabs.R"),
               local = TRUE,
               encoding = "UTF-8")
      }
      
    })
  
  #enable and disable submit contrast button
  observeEvent(
    eventExpr = input$Condition, {
      if(!is.null(input$Condition)){
        if(length(input$Condition) > 0  & input$Condition != "" ){
          shinyjs::enable("submit_contrast")
        }else{
          shinyjs::disable("submit_contrast")
        }
      }else{
        shinyjs::disable("submit_contrast")
      }
    })
  
  output$summary_gene_set <- renderUI({ 
    HTML(paste0(#"N (non-redundant) genes: ", rv$non_redundant_input_data_size,
                #"<br/>",
                "N (non-redundant) genes: ", dim(rv$input_data)[1]
    ))
  })
  
  #when i have defined the input data either by filters or by input
  observeEvent(
    eventExpr = rv$input_data, {
      print("in the observeEvent of rv$input_data")
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
        
       
        #print(paste("MSI input data", length(rv$input_data)))
        rv$dim_msi_df <- dim(rv$input_data)[1]
       
       
        
        ## UPDATE the GSO kegg map picker
        gso_kegg_maps = ko_Pathway %>% 
          filter(Name %in% rv$input_data$primary_name) %>% 
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


  
  ### NORM plots
  source(file = file.path(version, "server_NORM_plots.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  
  
  #GSSO KEGG_map_
  source(file = file.path(version, "server_KEGG_map_gso.R"),
         local = TRUE,
         encoding = "UTF-8")
  