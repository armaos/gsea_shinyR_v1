library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(DT)
#library(plotly)
library(shinyjs)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(ggupset )
#library(ggpubr)

if (!interactive()) sink(stderr(), type = "output")
wd = "/srv/shiny-server/gsea_shinyR/"



see_pathview <- function(..., save_image = FALSE)
{
  arguments <- list(...)
  pathview(gene.data = arguments$gene.data,
           gene.idtype = arguments$gene.idtype,
           pathway.id = arguments$pathway.id,
           species = arguments$species,
           bins = arguments$bins,
           limit = arguments$limit,
           low = arguments$low,
           mid = arguments$mid,
           high = arguments$high,
           na.col = arguments$na.col,
           kegg.native = arguments$kegg.native,
           kegg.dir = arguments$kegg.dir,
           same.layer = arguments$same.layer,
           sign.pos = arguments$sign.pos,
           new.signature =  arguments$new.signature ,
           out.suffix	= arguments$out.suffix
           
               
  )
  
  msg <- capture.output(pathview::pathview(...), type = "message")
  
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  file.remove(filename)
  img <- png::readPNG(filename)
  grid::grid.raster(img)
  
}

empty_image <- function(..., save_image = FALSE)
{
  print("empty")
  grid::grid.raster(png::readPNG("empty.png"))
}





server <- function(input, output, session) {
  
  # the annotation table 
  #annotations = as.data.frame(read.delim2("Ka.3.1_annotation_table_test.tsv", header = T, sep = "\t", fill = T))
  setwd(wd)
  
  annotations = as.data.frame(read.delim2("Ka.3.1_annotation_table.tsv", header = T, sep = "\t", fill = T))
  

  # these are the kegg table thingy
  
  koala_thresholds = as.data.frame(read.table("koala_genes_vector.txt"))
  koala_thresholds <- koala_thresholds %>%
    group_by(V1 ) %>% 
    summarise(V2 = max(V2))
  
  
  # MSI_D_L =  as.data.frame(read.table("MSI_D_L.txt", header = T))
  
  
  go_terms =  as.data.frame(read.delim2("GOIDs_GOterms.txt", header = T, sep = "\t", 
                                        stringsAsFactors = F, fill = T))
  go2gene_all <- as.data.frame(read.table("kal_v3_all_GO_test.txt", header = T))
  
  # here AnnotationDbi::Ontology("GO:0006486") for go ontology grouping
  #go2gene_all <- as.data.frame(read.table("kal_v3_all_GO.txt", header = T))
  #demo_genes <- as.data.frame(read.table("genes_to_test.txt", header = T))
  
  #### read KO annotations (i have manipulated the data to load faster)
  # ko_annots <- as.data.frame(read.table("Ka.2_plus_3.1_KO_terms.txt", header = T)) 
  # ko_annots %>%
  #   gather(key, val, c(koala_KEGG_ko, em_KEGG_ko)) %>%
  #   separate_rows(val, sep = "[,]") %>%
  #   select(-key) %>%
  #   filter(val != "-") %>%
  #   write.table("Ka.2_plus_3.1_KO_terms_manipulate.txt", quote = F, row.names = F)
  ko_annots <- as.data.frame(read.table("Ka.2_plus_3.1_KO_terms_manipulate.txt", header = T)) 
  
  
  ko_symbols <- as.data.frame(read.delim2("KO_titles_new.txt", header = F, sep = ";", quote = "")) 
  
  
  #### read kegg mapp annotations (i have manipulated the data to load faster
  # ko_Pathway = "Ka.2_plus_3.1_KO_Pathways.txt"
  # ko_Pathway <- as.data.frame(read.table(ko_Pathway, header = T))
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
    mutate(val = str_replace(val ,"ko", "")) 
  
  ######
  
  ko_Pathway_terms <- as.data.frame(read.table("Ka_KEGGpathways.txt", header = F, sep = "\t", stringsAsFactors = F, 
                                               colClasses = c("character", "character")))
  
  ko_nodes <- as.data.frame(read.table("KEGGmaps_and_nodes.txt", header = F, sep = "\t") ) %>% 
    rename(map = V1, nodes = V2) %>%
    mutate(map = str_replace(map, "ko", ""))
  
  
  
  rv <- reactiveValues(gsea_plots = NULL,
                       IP = NULL, 
                       
                       condition = "",
                       
                       data = NULL, 
                       text = NULL, 
                       okplot = FALSE, 
                       okplot_KEGGmap_ora = FALSE,
                       okplot_KEGGmap_gsea = FALSE,
                       okplot_gseaplot = FALSE,
                       map = "",
                       map_gsea = "",
                       go_term = "",
                       run_counter = 0, 
                       
                       gse = NULL,
                       gse_desc = NULL,
                       
                       barplot1 = NULL,
                       barplot2 = NULL,
                       dotplot = NULL,
                       emapplot = NULL,
                       cnetplot = NULL,
                       heatplot = NULL,
                       treeplot = NULL,
                       upsetplot = NULL,
                       
                       dotplot_gsea = NULL,
                       ridgeplot_gsea = NULL,
                       
                       gsea = NULL,
                       gseaplot = NULL,
                       
                       enriched_kegg_pathways = NULL,
                       #enriched_kegg_pathways_desc = NULL,
                       #enriched_kegg_pathways_id = NULL,
                       
                       input_data = NULL,
                       kegg_ora_vector = NULL,
                       kegg_gsea_vector = NULL,
                       go_term = NULL,
                       
                       ora = data.frame(), 
                       msi_df = data.frame(), 
                       gsea_input = NULL,
                       gsea = data.frame(),
                       gsea_kegg = data.frame(),
                       go = NULL,
                       enricher_plots = NULL,
                       
                       fc = c(),
                       
                       dim_msi_df = NULL ,

                       kegg_ora_col_low = "red",
                       kegg_ora_col_high = "green",
                       kegg_gsea_col_low =  "red",
                       kegg_gsea_col_high = "green"
                       
                       )
  
  # observeEvent(input$demo, {
  #   updateTextAreaInput(session, "input_text", 
  #                       value = paste(head(unique(demo_genes$Name) , 40), collapse = "\n") )
  #   updateSelectInput(session, "choose", selected = "text")
  #   
  # })
  # 
  # observeEvent(input$choose, {
  #   if(input$choose == "file"){
  #     if(is.null(input$data_file$datapath) || input$data_file$datapath == ""  ){
  #       shinyjs::disable("run")
  #     }else{
  #       shinyjs::enable("run")
  #     }
  #   }else if(input$choose == "text"){
  #     if(is.null(input$input_text) || input$input_text == ""){
  #       shinyjs::disable("run")
  #     }else{
  #       shinyjs::enable("run")
  #     }
  #   }
  # })
  # 
  observeEvent(input$data_file$datapath, {
    if(input$choose == "file"){
      if(!is.null(input$data_file$datapath) && input$data_file$datapath != "" ){
        print("file loaded")
        shinyjs::enable("submit_norm_and_de")
      }else{
        shinyjs::disable("submit_norm_and_de")
        
      }
    }
  })
  # 
  # observeEvent(input$input_text, {
  #   if(input$choose == "text"){
  #     if(!is.null(input$input_text) && input$input_text != ""){
  #       shinyjs::enable("run")
  #     }else{
  #       shinyjs::disable("run")
  #     }
  #   }
  # })
  
  observeEvent(input$submit_contrast, {
    shinyjs::disable("submit_contrast")
    shinyjs::disable("submit_norm_and_de")
    shinyjs::show("submitted_contrast1")
    shinyjs::show("submitted_contrast2")
    rv$condition  <- input$Condition
    
    rv$msi_df = as.data.frame(read.table(paste0("MSI_data/", rv$condition, "_all_values.csv"), header = T , sep = "\t"))
    colnames(rv$msi_df)[1] = "Name"
    rv$msi_df <- rv$msi_df %>% arrange(desc(logFC))
    rv$gsea_input = rv$msi_df$logFC
    rv$gsea_kegg = data.frame()
    names(rv$gsea_input) <- rv$msi_df$Name
    
    updatePickerInput(session, "overview_table_Default_values", choices = colnames(rv$msi_df), selected = c("Name", "logFC", "logCPM", "FDR"))
    rv$msi_df <- rv$msi_df %>% left_join(annotations)
    updatePickerInput(session, "overview_table_terms", choices = colnames(annotations)[29:length(colnames(annotations))])
    

    shinyjs::show("run_gsea_box")
    shinyjs::show("data_overview_box")


  })
  
  observeEvent(input$submit_norm_and_de, {
    rv$condition  <- input$Condition
    #shinyjs::disable("submit_contrast")
    #shinyjs::disable("submit_norm_and_de")
  })
  
  
  
  observeEvent(input$choose_map_ora, {
    rv$map <- subset(ko_Pathway_terms %>% filter(V2 == input$choose_map_ora))$V1 #str_replace(input$choose_map_ora, "ko", "")
    #print(paste0("observeEvent choose_map_ora rv$map-", rv$map,"-", is.null(rv$map),"-", is.na(rv$map), "-",length(rv$map)))
  })
  
  observeEvent(input$choose_map_gsea, {
    rv$map_gsea <- subset(ko_Pathway_terms %>% filter(V2 == input$choose_map_gsea))$V1 
    #print(paste0("observeEvent choose_map_gsea rv$map_gsea-", rv$map_gsea,"-", is.null(rv$map_gsea),"-", is.na(rv$map_gsea), "-",length(rv$map_gsea)))
  })
  
  observeEvent(input$choose_GO, {
    # rv$okplot_gseaplot <- TRUE
    rv$go_term <- input$choose_GO
    #rv$gseaplot = gseaplot(gsea, by = "all", title = rv$go_term, geneSetID = 1)
    #print(paste("observeEvent choose_GO rv$go_term",rv$okplot_gseaplot , rv$go_term))
  })
  
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
  
        # find the ko node id terms
        print("FIND ko_annots")
        print(paste("rv$msi_df", dim(rv$msi_df)))
        ko_annots <- ko_annots %>%
          inner_join(rv$msi_df %>% select(Name, logFC))
  
        print("FIND ko_annots DONE")
  
        gsea_kegg_v = ko_annots %>% distinct(val, logFC )
        rv$kegg_gsea_vector = gsea_kegg_v$logFC
        names(rv$kegg_gsea_vector) <- gsea_kegg_v$val
  
        print("GSEA")
        # first run the KEGG gsea just to have it stored so that the maps can be loaded
        if(dim(rv$gsea_kegg)[1] == 0){
          rv$gsea_kegg = GSEA(rv$gsea_input, TERM2GENE=ko_Pathway, TERM2NAME = ko_Pathway_terms)
          print(length(rv$gsea_kegg))
          updateSelectInput(session, "choose_map_gsea", choices = rv$gsea_kegg$Description[!is.na(rv$gsea_kegg$Description)])
          rv$map_gsea = ifelse(length(rv$gsea_kegg[, 1]) > 0 , rv$gsea_kegg[1, 1], rv$map_gsea)
          shinyjs::show("gsea_kegg_pathway_dropdown")
          print(paste("rv$map_gsea ", rv$map_gsea ))
        }

        print(paste("length rv$gsea_kegg" , length(rv$gsea_kegg)))

        if (input$enrichment_term_gsea %in% c("biological_process", "molecular_function", "cellular_component")){
          term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == input$enrichment_term_gsea))$GO_ID)
          term2name = go_terms %>% filter(GO_namespace == input$enrichment_term_gsea)
          
          rv$gsea <- GSEA(rv$gsea_input, 
                             TERM2GENE = term2gene,
                             TERM2NAME = term2name)
          
        }else if(input$enrichment_term_gsea == "KEGG_Pathways"){
          term2gene = ko_Pathway
          term2name = ko_Pathway_terms
          rv$gsea = rv$gsea_kegg
        }else if(input$enrichment_term_gsea == "Reactome"){
          term2gene = ko_Pathway
          term2name = ko_Pathway_terms
        }
        
        
        


        rv$dotplot_gsea = dotplot(rv$gsea)
        print("doing ridgeplot_gsea")
        rv$ridgeplot_gsea = ridgeplot(rv$gsea )
        print("ridgeplot_gsea DONE")
        updateSelectInput(session, "choose_GO", choices = rv$gsea$Description[!is.na(rv$gsea$Description)])
        #updateCheckboxGroupInput(session, "custom_upset_terms", choices = rv$gsea$Description[!is.na(rv$gsea$Description)])
        updatePickerInput(session, "custom_upset_terms", choices = rv$gsea$Description[!is.na(rv$gsea$Description)])
        shinyjs::show("gsea_plots_tabsetPanel")
        
        shinyjs::show("custom_upset_box")
        
        shinyjs::hide("text_KEGGscape_GSEA")
        shinyjs::show("gsea_kegg_pathway_table")
        
        
        
  
       
        print("GSEA: DONE")
  
        rv$okplot_KEGGmap_gsea <- TRUE
        rv$okplot_gseaplot <- TRUE
      }
    })
  
  observeEvent(
    eventExpr = input$run, {
      
      if (rv$condition == ""){
        shinyjs::show("text_configure_ora")
      }else{
        shinyjs::hide("text_configure_ora")
        shinyjs::disable("run")
        
        html("text_ora_plots", "You should submit your ploting configuration first and Plot your table!!")
        shinyjs::show("text_ora_plots")
        
        html("text_KEGGscape_ORA", "Processing...")
        shinyjs::show("text_KEGGscape_ORA")
        
        # rv$data = switch(input$choose,
        #                  file = read_file(input$data_file$datapath),
        #                  text = input$input_text
        # ) 
        
        rv$IP = input$getIP
        rv$run_counter = rv$run_counter +1
        #rv$okplot <- FALSE
        rv$okplot_KEGGmap_ora <- FALSE
  
        rv$barplot1 <- NULL
        rv$barplot2 <- NULL
        rv$dotplot <- NULL
        rv$emapplot <- NULL
        rv$cnetplot <- NULL
        rv$heatplot <- NULL
        rv$treeplot <- NULL
        rv$upsetplot <- NULL
  
        rv$input_data <- NULL
        rv$ora <- data.frame()
        rv$ora_desc <- NULL
        rv$map <- ""
        rv$go_term <- ""
        rv$gsea <- data.frame()
        rv$fc = c()
        
        

        #shinyjs::hide("kegg_pval")
        
        
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
        
        # this is the reading of MSI data
        fc_filter_low = input$slider_lfc[1] 
        fc_filter_high = input$slider_lfc[2] 
        
        #pval_filter = input$slider_pval[1] 
        fdr_filter = input$slider_fdr[1] 
        
        #print(paste(fc_filter_low, fc_filter_high, pval_filter, paste0("MSI_data/", input$Condition, "_all_values.csv")))
        
        msi_df <- rv$msi_df %>%  filter(abs(logFC) >= fc_filter_low & abs(logFC) <= fc_filter_high & FDR <= 10^(-fdr_filter))
        
        # msi_df = subset(as.data.frame(read.table(paste0("MSI_data/", isolate(input$Condition), "_all_values.csv"), header = T , sep = "\t")))$Name
        
        rv$input_data <- msi_df$Name
        rv$fc = msi_df$logFC
        names(rv$fc)= rv$input_data
        
        # print(paste("rv$msi_df ", head(subset(as.data.frame(read.table(paste0("MSI_data/", isolate(input$Condition), "_all_values.csv"), header = T , sep = "\t")))$Name ,1)))
        # print(paste("dim_msi_df", rv$dim_msi_df))
        print(paste("MSI input data", length(rv$input_data)))
        rv$dim_msi_df <- length(rv$input_data)
        # GSEA with enricher
        
        
        #rv$ora <- enricher(rv$input_data, TERM2GENE=go2gene_all, TERM2NAME = go_terms)
        # if(!is.null(rv$ora)){
        #   rv$ora_desc = head(rv$ora[rv$ora$pvalue <= 10^(-input$go_pval), 2] ,20)
        #   if(length(rv$ora_desc) > 0 ){
        #     
        #     rv$barplot1 = barplot(rv$ora, showCategory=rv$ora_desc, main = "Barplot (count)" )
        #     rv$barplot2 = mutate(rv$ora,
        #                          qscore = -log(p.adjust, base=10)) %>%
        #       barplot(x="qscore", showCategory=rv$ora_desc, main ="Barplot (qscore)")
        #     rv$dotplot = dotplot(rv$ora, showCategory=rv$ora_desc)
        #     p_t = pairwise_termsim(rv$ora)
        #     
        #     rv$emapplot = emapplot(p_t, showCategory = rv$ora_desc)
        #     rv$cnetplot = cnetplot(rv$ora, categorySize="pvalue", showCategory = rv$ora_desc, foldChange = fc )
        #     rv$heatplot = heatplot(rv$ora, showCategory = rv$ora_desc, foldChange = fc )
        #     rv$upsetplot = upsetplot(rv$ora, showCategory = rv$ora_desc)
        #     rv$treeplot = treeplot(p_t)
        #     
        #   }    
        #   
        # }
        
        print("ORA with enricher: DONE")
        
        # find the enriched kegg paths
        print("KEGG with enricher")
        
        rv$enriched_kegg_pathways = enricher(rv$input_data, TERM2GENE=ko_Pathway, TERM2NAME = ko_Pathway_terms)
        #rv$enriched_kegg_pathways_desc = rv$enriched_kegg_pathways[rv$enriched_kegg_pathways$pvalue <= 10^(-input$kegg_pval), 2]
        #rv$enriched_kegg_pathways_id = rv$enriched_kegg_pathways[rv$enriched_kegg_pathways$pvalue <= 10^(-input$kegg_pval), 1]
        
        updateSelectInput(session, "choose_map_ora", choices = rv$enriched_kegg_pathways[, 2])
        rv$map <- ifelse(length(rv$enriched_kegg_pathways[, 1]) > 0 , rv$enriched_kegg_pathways[1, 1], rv$map)
        shinyjs::show("ora_kegg_pathway_dropdown")

        print("KEGG with enricher: DONE")
        
        
        ### ########
        
     
        
        # find the ko node id terms
        print("FIND ko_annots")
        ora_kegg <- ko_annots %>%
          inner_join(msi_df %>% select(Name, logFC)) %>% 
          filter(Name %in% rv$input_data) %>%
          distinct(val, logFC)
        rv$kegg_ora_vector <- ora_kegg$logFC
        names(rv$kegg_ora_vector) <- ora_kegg$val
        #print(paste("kegg_ora_vector",head(rv$kegg_ora_vector)))
        
        #rv$okplot <- TRUE
        rv$okplot_KEGGmap_ora <- TRUE 
        print("and now show!")
        shinyjs::show("box_configure_term_ora")
        shinyjs::show("enrichment_term")
        shinyjs::show("submit_enrichment")
        shinyjs::hide("text_KEGGscape_ORA")
        shinyjs::show("ora_kegg_pathway_table")
        
        
      }
    })
  
  output$txt <- renderText({
    paste("#Records:",  rv$dim_msi_df)
  })
  
  
  
  observeEvent(
    eventExpr = input$submit_enrichment, {
      print(" input$submit_enrichment")
      
      html("text_ora_go_table", "Processing...")
      shinyjs::show("text_ora_go_table")
      print("ORA with enricher")
      if (input$enrichment_term %in% c("biological_process", "molecular_function", "cellular_component")){
        term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == input$enrichment_term))$GO_ID)
        term2name = go_terms %>% filter(GO_namespace == input$enrichment_term)
        
        rv$ora <- enricher(rv$input_data, 
                           TERM2GENE = term2gene,
                           TERM2NAME = term2name)
        
      }else if(input$enrichment_term == "KEGG_Pathways"){
        rv$ora = rv$enriched_kegg_pathways
      }else if(input$enrichment_term == "Reactome"){
        rv$ora = rv$enriched_kegg_pathways
      }
     
      
      print("ORA with enricher :DONE")
      shinyjs::show("box_ora_send_to_plot")
    })
  
  observeEvent(
    eventExpr = input$submit_plot_ora, {
      html("text_ora_plots", "Processing...!")
      
      print(" input$submit_plot_ora")
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

          rv$barplot1 = barplot(rv$ora, showCategory=rv$ora_desc, main = "Barplot (count)" )
          rv$barplot2 = mutate(rv$ora,
                               qscore = -log(p.adjust, base=10)) %>%
            barplot(x="qscore", showCategory=rv$ora_desc, main ="Barplot (qscore)")
          rv$dotplot = dotplot(rv$ora, showCategory=rv$ora_desc)
          p_t = pairwise_termsim(rv$ora)
          
          rv$emapplot = emapplot(p_t, showCategory = rv$ora_desc)
          rv$cnetplot = cnetplot(rv$ora, categorySize="pvalue", showCategory = rv$ora_desc, foldChange = rv$fc )
          rv$heatplot = heatplot(rv$ora, showCategory = rv$ora_desc, foldChange = rv$fc )
          rv$upsetplot = upsetplot(rv$ora, showCategory = rv$ora_desc)
          rv$treeplot = treeplot(p_t)
          shinyjs::hide("text_ora_plots")
          shinyjs::show("ora_plots_tabsetPanel")
          
         
        }else{
          html("text_ora_plots", "<strong>No enriched terms, try with a different threshold!</strong>")
        }

      }
    })
  
  
  # output$enriched_go <- renderPlotly({
  #   if(rv$okplot){
  #     shinyjs::enable("run")
  #     if (!is.null(rv$data ) && rv$data != ""){
  #       draw_enrich_go(rv$go, max_go = 20)
  #       shinyjs::hide("text1")
  #     }else{
  #       html("text1", "<strong>bold</strong> WRONG input data")
  #     }
  #   }
  # })
  
  
  # output$enriched_map <- renderPlot({
  #   if(rv$okplot){
  #     shinyjs::enable("run")
  #     if (!is.null(rv$data ) && rv$data != ""){
  #       draw_enrich_go_map(rv$go)
  #       shinyjs::hide("text1")
  #     }else{
  #       html("text1", "<strong>bold</strong> WRONG input data")
  #     }
  #   }
  # })
  
  
  # output$enriched_map <- renderPlot({
  #   if(rv$okplot){
  #     shinyjs::enable("run")
  #     if (!is.null(rv$data ) && rv$data != ""){
  #       
  #       hist(rnorm(100, 4, 1),breaks = 50)
  #       shinyjs::hide("text1")
  #       shinyjs::hide("text_gsea_diana")
  #     }else{
  #       html("text1", "<strong>bold</strong> WRONG input data")
  #       html("text_gsea_diana", "<strong>bold</strong> WRONG input data")
  #     }
  #   }
  # })
  
  #TABLES

  output$data_overview <- DT::renderDT({
    colnames_to_show = c(input$overview_table_Default_values,  input$overview_table_terms)
    rv$msi_df %>% select(colnames_to_show)
  })
  

  output$ora_go_table <- DT::renderDT({
    if(length(rv$ora ) > 0){
    #if(!is.null(rv$ora) & rv$ora != ""){
      shinyjs::enable("run")
      shinyjs::hide("text_ora_go_table")
      shinyjs::show("submit_plot_ora")
      shinyjs::show("ORA_max_elements")
      

      as.data.frame(rv$ora) %>% select(-geneID) 
      #as.data.frame(rv$ora[rv$ora$pvalue <= 10^(-isolate(input$go_pval))])  %>% select(-geneID)
      }
    
  },
  callback = JS(
    # "table.on( 'search.dt', function () { Shiny.setInputValue( 'search', table.search() ); } );",
    "$( 'input').on( 'input', function(){",
    "var value = this.value;",
    "var clicked_td = $(this).closest('td');",
    "var td_index = clicked_td.index(); ",
    "var tr_index = clicked_td.parent().index();",
    #"console.log( 'search11', td_index , tr_index, table.column(td_index).header().textContent);",
    "Shiny.setInputValue( 'ORA_' + table.column( td_index ).header().textContent , value);",
    "});"
  ),
  filter="top",
  options = list(
    dom = '<"top" pif>',
    pageLength = 10)
  )
  
  
  # output$ora_kegg_table <- DT::renderDT({
  #   if(!is.null(rv$enriched_kegg_pathways)){
  #     shinyjs::enable("run")
  #     shinyjs::hide("text_ora_kegg_table")
  #     shinyjs::show("kegg_pval")
  #     as.data.frame(rv$enriched_kegg_pathways[rv$enriched_kegg_pathways$pvalue <= 10^(-input$kegg_pval)])  %>% select(-geneID)}
  # },
  # filter="top",
  # options = list(
  #   dom = '<"top" pif>',
  #   pageLength = 10)
  # )
  
  
  output$gsea_table <- DT::renderDT({
    if(length(rv$gsea ) > 0 ){
      shinyjs::enable("run")
      shinyjs::hide("text_gsea_table")
      as.data.frame(rv$gsea) %>% select(-core_enrichment)}
  },
  filter="top",
  options = list(
    dom = '<"top" pif>',
    pageLength = 10)
  )
  
  
  # ORA plots  
  output$barplot1 <- renderPlot({
    if (!is.null(rv$barplot1) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_barplot1_text")
      rv$barplot1
    }else{
      empty_image()
    }
  })
  
  output$barplot2 <- renderPlot({
    if (!is.null(rv$barplot2) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_barplot2_text")      
      rv$barplot2
    }else{
      empty_image()
    }
  })

    output$dotplot <- renderPlot({
    if(!is.null(rv$dotplot) & !is.null(rv$ora_desc)){
      shinyjs::enable("run")
      shinyjs::hide("ora_dotplot_text")
      rv$dotplot
    }else{
      shinyjs::enable("run")
      html("text_ora_plots", "<strong>Make your configurarion for plotting from the previous tab panel first!</strong>")
      empty_image()
    }
  })
  
  output$emapplot <- renderPlot({
    if (!is.null(rv$emapplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_emapplot_text")

      rv$emapplot
    }else{
      empty_image()
    }
  })
  
  output$cnetplot <- renderPlot({
    if (!is.null(rv$cnetplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_cnetplot_text")

      rv$cnetplot
    }else{
      empty_image()
    }
  })
  output$heatplot <- renderPlot({
    if (!is.null(rv$heatplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_heatplot_text")

      rv$heatplot
    }else{
      empty_image()
    }
  })
  
  output$treeplot <- renderPlot({
    if (!is.null(rv$treeplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_treeplot_text")

      rv$treeplot
    }else{
      empty_image()
    }
  })
  
  
  output$upsetplot <- renderPlot({
    if (!is.null(rv$upsetplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_upsetplot_text")

      rv$upsetplot
    }else{
      empty_image()
    }
  })
  
  #GSEA overview
  output$dotplot_gsea <- renderPlot({
    if (!is.null(rv$dotplot_gsea) & !is.null(rv$gsea)){
      shinyjs::hide("text_gsea_plot_overview")
      shinyjs::enable("run")
      shinyjs::hide("gsea_dotplot_text")

      rv$dotplot_gsea
    }else{
      empty_image()
    }
  })
  
  output$ridgeplot_gsea <- renderPlot({
    if (!is.null(rv$ridgeplot_gsea) & !is.null(rv$gsea)){
      shinyjs::hide("text_gsea_plot_overview")
      shinyjs::enable("run")
      shinyjs::hide("gsea_ridgeplot_text")
      
      rv$ridgeplot_gsea
      
    }else{
      empty_image()
    }
  })
  # output$treeplot <- renderPlot({
  #   if (!is.null(rv$treeplot)){
  #     shinyjs::hide("text_ora_plots")
  #     rv$treeplot
  #   }
  # })
  
  #GSEA PLOT
  output$gseaplot <- renderPlot({
    if(rv$okplot_gseaplot){
      shinyjs::enable("run")
      #rv$gsea_plots[rv$go_term]
      shinyjs::hide("text_gsea_plot")
      shinyjs::show("choose_GO_box")
      # if ((!is.null(rv$data ) && rv$data != "" ) && (!is.null(rv$go_term ) && rv$go_term != ""  && rv$go_term != "NA" && !is.na(rv$go_term ) ) ){
      #   ggarrange(plotlist=gsea_plots, col = 1)
      #   print("in gseaplot")
      #   shinyjs::hide("text_gsea_plot")
      # }
      if( rv$go_term != "" && length(rv$go_term) > 0 && !is.na(rv$go_term )){  
        print(paste("output$gseaplot" ,rv$go_term ))
        print(rv$gsea[1,1])
        print(which(rv$gsea$Description == rv$go_term))
        gseaplot2(rv$gsea, geneSetID = which(rv$gsea$Description == rv$go_term) , title = rv$go_term)
      }
    }
  })
  
  
  
  my_custom_upset_plot <- eventReactive(input$submit_custom_upset_plot,{
    as.data.frame(isolate(rv$gsea)) %>% 
      filter(Description %in% isolate(input$custom_upset_terms)) %>%
      select(Description, core_enrichment) %>% 
      separate_rows(core_enrichment, sep = "/") %>% 
      mutate(Description = ifelse(str_length(Description) > 30, paste0(str_sub(Description,1,30), "..."), Description)) %>%
      group_by(core_enrichment) %>% 
      summarise(terms = list(Description)) %>% 
      ggplot(aes(x=terms)) +
      geom_bar() +
      scale_x_upset()
  })
  output$gsea_custom_upsetplot <- renderPlot({
    my_custom_upset_plot()
  })
  
  
  ## ENRICHER kegg table and plot
  output$ora_keggscape_table <- DT::renderDT({
    #shinyjs::hide("text_KEGGscape_ORA")
    #shinyjs::show("ora_kegg_pathway_table")
    if( rv$map != "" && length(rv$map) > 0 && !is.na(rv$map )){  
      ko_nodes %>% 
        filter(map == rv$map) %>%
        #filter(map == "00630") %>%
        separate_rows(nodes, sep = ",") %>%
        left_join(ko_symbols, by = c("nodes" = "V1")) %>%
        inner_join(ko_annots, by=c("nodes" = "val")) %>%
        filter(nodes %in% names(rv$kegg_ora_vector)) %>%
        filter(Name %in% rv$input_data) %>%
        group_by(nodes, V2, V3) %>% 
        summarise(Transcripts = paste(unique(Name), collapse= "; ")) %>%
        rename(Symbol = V2, KEGG_ko = nodes) 
    }
  },
  filter="top",
  options = list(
    pageLength = 10)
  )
  

   observeEvent(
    eventExpr = input$submit_kegg_ora_toggle, {
      rv$kegg_ora_col_low = ifelse(rv$kegg_ora_col_low == "red", "green", "red")
      rv$kegg_ora_col_high = ifelse(rv$kegg_ora_col_high == "red", "green", "red")
    })
  
  observeEvent(
    eventExpr = input$submit_kegg_gsea_toggle, {
      rv$kegg_gsea_col_low = ifelse(rv$kegg_gsea_col_low == "red", "green", "red")
      rv$kegg_gsea_col_high = ifelse(rv$kegg_gsea_col_high == "red", "green", "red")
    })

  output$KEGG_map <- renderImage({
    outfile <- tempfile(fileext='.png')
    png(outfile, width=400, height=300)
    frame()
    dev.off()
    blank <- list(src = outfile,
                  contentType = 'image/png',
                  width = 400,
                  height = 300,
                  alt = " ")
    
    if(rv$okplot_KEGGmap_ora){
      shinyjs::enable("run")
      #if( !is.null(rv$map) && rv$map != ""  && rv$map != "NA" && !is.na(rv$map ) ){
      if( rv$map != "" && length(rv$map) > 0 && !is.na(rv$map )){  
        # hist(rnorm(100, 4, 1),breaks = 50, main = rv$run_counter)
        
        tempFolder <- tempdir() # tempFolder = "temp";
        
        kegg_ora_vector = rv$kegg_ora_vector[names(rv$kegg_ora_vector) %in% str_split(subset(ko_nodes %>% filter(map == rv$map))$nodes, pattern = "[,]")[[1]]]
        
        max_range = ceiling(max(abs(min(kegg_ora_vector)), abs( max(kegg_ora_vector)) ))
        
        outfile <- paste( tempFolder,"/ko", rv$map, ".", rv$IP, ".ora.png",sep="")
        wd = getwd()
        setwd(tempFolder)
        
        pvout = pathview(gene.data = kegg_ora_vector,
                         gene.idtype = "KEGG",
                         pathway.id = rv$map,
                         species = "ko",
                         bins = list(gene = 8),
                         #limit = list(gene = c(0,1)),
                         limit = list(gene = c(-max_range, max_range)),
                         low = list(gene = rv$kegg_ora_col_low  ),
                         #low = list(gene = "red"  ),
                         mid = list(gene = "yellow"),
                         #high = list(gene =  "green"),
                         high = list(gene =  rv$kegg_ora_col_high),
                         na.col = "transparent",
                         kegg.dir = file.path(wd, "kegg_data"),
                         #out.suffix = suffix,
                         kegg.native = T,
                         same.layer = T,
                         sign.pos = "bottomleft",
                         new.signature = F,
                         #key.pos = "topleft",
                         out.suffix = paste0(rv$IP,".ora"))
        setwd(wd)
        shinyjs::show("keggscape_ora_conf")        

        print(outfile)
        return( list(src = outfile,
                     contentType = 'image/png',
                     idth = "150%",
                     height = "150%",
                     alt = "This is alternate text") )
        
        
      }else{
        print("empty images KEGG")
        html("text_KEGGscape_ORA", "<strong>No enriched KEGG terms! Try again with different filters</strong>")
        return (blank)
        #empty_image()
      }
    }else{
      return (blank)
    }
  }, deleteFile = T)
  
  
  ### and the gsea kegg table and kegg scape
  output$gsea_keggscape_table <- DT::renderDT({
    if( rv$map_gsea != "" && length(rv$map_gsea) > 0 && !is.na(rv$map_gsea )){  
      ko_nodes %>% 
        filter(map == rv$map_gsea) %>%
        separate_rows(nodes, sep = ",") %>%
        left_join(ko_symbols, by = c("nodes" = "V1")) %>%
        inner_join(ko_annots, by=c("nodes" = "val")) %>%
        filter(nodes %in% names(rv$kegg_gsea_vector)) %>%
        group_by(nodes, V2, V3) %>% 
        summarise(Transcripts = paste(unique(Name), collapse= "; ")) %>%
        rename(Symbol = V2, KEGG_ko = nodes) 
    }
  },
  filter="top",
  options = list(
    pageLength = 10)
  )
  
  #the pathview GSEA
  output$KEGG_map_gsea <- renderImage({
    outfile <- tempfile(fileext='.png')
    png(outfile, width=400, height=300)
    frame()
    dev.off()
    blank <- list(src = outfile,
                  contentType = 'image/png',
                  width = 400,
                  height = 300,
                  alt = " ")
    
    if(rv$okplot_KEGGmap_gsea){
      shinyjs::enable("run")
      print(paste0("output$KEGG_map_gsea  rv$map_gsea:",rv$map_gsea,"-",is.null(rv$map_gsea)))
      if( rv$map_gsea != "" && length(rv$map_gsea) > 0 && !is.na(rv$map_gsea )){  
        # hist(rnorm(100, 4, 1),breaks = 50, main = rv$run_counter)
        
        tempFolder <- tempdir() # tempFolder = "temp";
        #pathID = rv$map_gsea
        #myIP = rv$IP
        #kegg_gsea_vector = rv$kegg_gsea_vector
        
        kegg_gsea_vector = rv$kegg_gsea_vector[names(rv$kegg_gsea_vector) %in% str_split(subset(ko_nodes %>% filter(map == rv$map_gsea))$nodes, pattern = "[,]")[[1]]]
        
        max_range = ceiling(max(abs(min(kegg_gsea_vector)), abs( max(kegg_gsea_vector)) ))
        
        outfile <- paste( tempFolder,"/ko", rv$map_gsea, ".", rv$IP, ".gsea.png",sep="")
        setwd(tempFolder)
        
        pvout = pathview(gene.data = kegg_gsea_vector,
                         gene.idtype = "KEGG",
                         pathway.id = rv$map_gsea,
                         species = "ko",
                         bins = list(gene = 8),
                         #limit = list(gene = c(0,1)),
                         limit = list(gene = c(-max_range, max_range)),
                         low = list(gene = rv$kegg_gsea_col_low  ),
                         mid = list(gene = "yellow"),
                         high = list(gene =  rv$kegg_gsea_col_high),
                         na.col = "transparent",
                         kegg.dir = file.path(wd, "kegg_data"),
                         #out.suffix = suffix,
                         kegg.native = T,
                         same.layer = T,
                         sign.pos = "bottomleft",
                         new.signature = F,
                         #key.pos = "topleft",
                         out.suffix = paste0(rv$IP,".gsea"))
        setwd(wd)
        shinyjs::show("keggscape_gsea_conf")
        
        print(outfile)
        return( list(src = outfile,
                     contentType = 'image/png',
                     idth = "150%",
                     height = "150%",
                     alt = "This is alternate text") )
        
        
      }else{
        print("empty images KEGG")
        html("text_KEGGscape_GSEA", "<strong>No enriched KEGG terms! Try again with different filters</strong>")
        return (blank)
        #empty_image()
      }
    }else{
      return (blank)
    }
  }, deleteFile = T)
  
}

