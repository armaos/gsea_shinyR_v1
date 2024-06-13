## ENRICHER kegg table and plot
output$gso_keggscape_table <- DT::renderDT({
  #shinyjs::hide("text_KEGGscape_ORA")
  #shinyjs::show("ora_keggscape_table_box")
  if( rv$map != "" && length(rv$map) > 0 && !is.na(rv$map )){  
    ko_nodes %>% 
      filter(map == rv$map) %>%
      #filter(map == "00590") %>%
      separate_rows(nodes, sep = ",") %>%
      left_join(ko_symbols, by = c("nodes" = "V1")) %>%
      inner_join(ko_annots, by=c("nodes" = "KO")) %>%
      filter(Name %in% representatives$Name) %>%
      mutate(Name = paste0(Name, " (KFK:", koala_threshold, ", em:", em, ")" )) %>%
      #filter(nodes %in% names(rv$kegg_gso_vector)) %>%
      #filter(Name %in% rv$input_data) %>%
      group_by(nodes, Symbol, Node_Description) %>% 
      summarise(Transcripts = paste(unique(Name), collapse= "; ")) %>%
      ungroup() %>%
      rename(KEGG_ko = nodes) 
  }
},
filter="top",
options = list(
  scrollX = TRUE, 
  scrollY = "400px",
  dom = '<"top" lpift>')
)


observeEvent(
  eventExpr = input$submit_kegg_gso_toggle, {
    rv$kegg_gso_col_low = ifelse(rv$kegg_gso_col_low == "red", "green", "red")
    rv$kegg_gso_col_high = ifelse(rv$kegg_gso_col_high == "red", "green", "red")
  })


observeEvent(
  eventExpr = input$submit_kegg_gsea_toggle, {
    rv$kegg_gsea_col_low = ifelse(rv$kegg_gsea_col_low == "red", "green", "red")
    rv$kegg_gsea_col_high = ifelse(rv$kegg_gsea_col_high == "red", "green", "red")
  })

output$KEGG_map_gso <- renderImage({
  outfile <- tempfile(fileext='.png')
  png(outfile, width=400, height=300)
  frame()
  dev.off()
  blank <- list(src = outfile,
                contentType = 'image/png',
                width = 400,
                height = 300,
                alt = " ")
  
  if(rv$okplot_KEGGmap_gso){
    
    
    if( rv$map != "" && length(rv$map) > 0 && !is.na(rv$map )){  
      tempFolder <- tempdir() # tempFolder = "temp";
      #tempFolder <- "~/OneDrive - Fondazione Istituto Italiano Tecnologia/B_CRO_projects/re_algae_c026/gsea_shinyR"
      # this is the previous implementation  
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
      # 
      # 
      # kegg_gso_vector = rv$kegg_gso_vector[names(rv$kegg_gso_vector) %in% str_split(subset(ko_nodes %>% filter(map == rv$map))$nodes, pattern = "[,]")[[1]]]
      # 
      # max_range = ceiling(max(abs(min(kegg_gso_vector)), abs( max(kegg_gso_vector)) ))
      # ##
      #
      kegg_gso_vector <- ko_nodes %>% 
        #filter(map == "05221") %>% 
        filter(map == rv$map) %>% 
        select(nodes) %>% 
        separate_rows(nodes) %>%
        distinct(nodes) %>%
        left_join(ko_annots %>% select(KO, Name), by = c("nodes" = "KO")) %>%
        #filter(!is.na(koala_threshold)) %>%
        inner_join(rv$msi_df %>% select(Name, contains("logFC")) %>% separate_rows(Name, sep = ",") ) %>%
        select(-Name) %>%
        group_by(nodes) %>%
        summarise_all(~max(., na.rm = T)) %>%
        left_join(koala_thresholds , by = c("nodes" = "KO")) 
      
      
      
      # kegg_gso_vector <- ko_nodes %>% 
      #   filter(map == rv$map) %>% 
      #   #filter(map == "05221") %>% 
      #   select(nodes) %>% 
      #   separate_rows(nodes) %>%
      #   distinct(nodes) %>%
      #   left_join(ko_annots, by = c("nodes" = "KO")) %>%
      #   left_join(rv$msi_df %>% select(Name, contains("logFC")) %>% separate_rows(Name, sep = ",") ) %>%
      #   #left_join(msi_df %>% select(Name, logFC) %>% separate_rows(Name, sep = ",") ) %>%
      #   group_by(nodes) %>%
      #   summarise(logFC = max(logFC, na.rm = T)) %>%
      #   ungroup() %>%
      #   mutate(logFC = ifelse(is.finite(logFC), logFC, NA)) %>%
      #   left_join(koala_thresholds , by = c("nodes" = "KO")) 
      
      kegg_gso_vector$koala_threshold = kegg_gso_vector$koala_threshold / 100
      print(paste("length kegg_gso_vector", dim(kegg_gso_vector)[1]))
      print(paste("input$kegg_gso_plotting_radio", input$kegg_gso_plotting_radio))
      
      #select the values to overlay in the kegg pathway
      if (input$kegg_gso_plotting_radio == "annotations"){
        #kegg_gso_vector$logFC = kegg_gso_vector$logFC / max(abs(kegg_gso_vector$logFC), na.rm = T)
        #gene_vector =  cbind(matrix(kegg_gso_vector$logFC), kegg_gso_vector$koala_threshold, kegg_gso_vector$em)
        gene_vector =  cbind(matrix(kegg_gso_vector$koala_threshold), kegg_gso_vector$em)
        rownames(gene_vector) =  kegg_gso_vector$nodes
        
        print(paste("server kegg map gso", max(kegg_gso_vector$koala_threshold, na.rm = T)))
        max_range = 1
        min_range = 0
      }else if (input$kegg_gso_plotting_radio == "conditions"){
        gene_vector =  cbind(as.matrix(kegg_gso_vector  %>% select(contains("logFC")) ))
        rownames(gene_vector) = kegg_gso_vector$nodes
        max_range = ceiling(max(abs(min(gene_vector, na.rm = T)), abs( max(gene_vector, na.rm = T)) ))
        min_range = -max_range
        
        print(head(gene_vector))
      }
      
      # if (input$kegg_gso_plotting_values == "logfc"){
      #   gene_vector =  kegg_gso_vector$logFC
      #   names(gene_vector) = kegg_gso_vector$nodes
      #   max_range = ceiling(max(abs(min(kegg_gso_vector$logFC, na.rm = T)), abs( max(kegg_gso_vector$logFC, na.rm = T)) ))
      #   min_range = -max_range
      # }else if (input$kegg_gso_plotting_values == "em"){
      #   gene_vector =  kegg_gso_vector$em
      #   names(gene_vector) = kegg_gso_vector$nodes
      #   max_range = 1
      #   min_range = 0
      # }else if (input$kegg_gso_plotting_values == "kfk"){
      #   gene_vector =  kegg_gso_vector$koala_threshold
      #   names(gene_vector) = kegg_gso_vector$nodes
      #   max_range = 1
      #   min_range = 0
      # }else if (input$kegg_gso_plotting_values == "all"){
      #   kegg_gso_vector$logFC = kegg_gso_vector$logFC / max(abs(kegg_gso_vector$logFC), na.rm = T)
      #   gene_vector =  cbind(matrix(kegg_gso_vector$logFC), kegg_gso_vector$koala_threshold, kegg_gso_vector$em)
      #   rownames(gene_vector) =  kegg_gso_vector$nodes
      #   max_range = 1
      #   min_range = -1
      # }
      
      print(paste("gene_vector", dim(gene_vector), collapse = ","))
      
      
      #print(gene_vector)
      
      
      outfile <- ifelse(dim(gene_vector)[2] > 1,
                        paste( tempFolder,"/ko", rv$map, ".", rv$IP, ".gso.multi.png",sep=""),
                        paste( tempFolder,"/ko", rv$map, ".", rv$IP, ".gso.png",sep=""))
      #outfile <- paste( tempFolder,"/ko", rv$map, ".", rv$IP, ".gso.png",sep="")
      wd = getwd()
      setwd(tempFolder)
      
      pvout = pathview(gene.data = gene_vector,
                       gene.idtype = "KEGG",
                       pathway.id = rv$map,
                       species = "ko",
                       bins = list(gene = 8),
                       #limit = list(gene = c(0,1)),
                       limit = list(gene = c(min_range, max_range)),
                       low = list(gene = rv$kegg_gso_col_low  ),
                       #low = list(gene = "red"  ),
                       mid = list(gene = "yellow"),
                       #high = list(gene =  "green"),
                       high = list(gene =  rv$kegg_gso_col_high),
                       na.col = "transparent",
                       kegg.dir = file.path(wd, "kegg_data"),
                       #out.suffix = suffix,
                       kegg.native = T,
                       same.layer = T,
                       sign.pos = "bottomleft",
                       new.signature = F,
                       #key.pos = "topleft",
                       out.suffix = paste0(rv$IP,".gso"))
      setwd(wd)
      shinyjs::show("keggscape_gso_conf")        
      
      print(outfile)
      return( list(src = outfile,
                   contentType = 'image/png',
                   width = "150%",
                   height = "150%",
                   alt = "This is alternate text") )
      
      
    }else{
      print("empty images KEGG")
      html("text_KEGGscape_gso", "<strong>No KEGG terms! Try again with different filters</strong>")
      return (blank)
      #empty_image()
    }
  }else{
    return (blank)
  }
}, deleteFile = F)
