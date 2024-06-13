### and the gsea kegg table and kegg scape
# output$gsea_keggscape_table <- DT::renderDT({
#   if( rv$map_gsea != "" && length(rv$map_gsea) > 0 && !is.na(rv$map_gsea )){  
#     ko_nodes %>% 
#       filter(map == rv$map_gsea) %>%
#       separate_rows(nodes, sep = ",") %>%
#       left_join(ko_symbols, by = c("nodes" = "V1")) %>%
#       inner_join(ko_annots, by=c("nodes" = "KO")) %>%
#       filter(nodes %in% names(rv$kegg_gsea_vector)) %>%
#       mutate(Name = paste0(Name, " (KFK:", koala_threshold, ", em:", em, ")" )) %>%
#       group_by(nodes, Symbol, Node_Description) %>% 
#       summarise(Transcripts = paste(unique(Name), collapse= "; ")) %>%
#       rename(KEGG_ko = nodes) 
#   }
# },
# filter="top",
# options = list(
#   scrollX = TRUE, 
#   scrollY = "400px",
#   dom = '<"top" lpift>')
# )

#the pathview GSEA
# output$KEGG_map_gsea <- renderImage({
#   outfile <- tempfile(fileext='.png')
#   png(outfile, width=400, height=300)
#   frame()
#   dev.off()
#   blank <- list(src = outfile,
#                 contentType = 'image/png',
#                 width = 400,
#                 height = 300,
#                 alt = " ")
#   
#   if(rv$okplot_KEGGmap_gsea){
#     shinyjs::enable("run")
#     print(paste0("output$KEGG_map_gsea  rv$map_gsea:",rv$map_gsea,"-",is.null(rv$map_gsea)))
#     if( rv$map_gsea != "" && length(rv$map_gsea) > 0 && !is.na(rv$map_gsea )){  
#       # hist(rnorm(100, 4, 1),breaks = 50, main = rv$run_counter)
#       
#       tempFolder <- tempdir() # tempFolder = "temp";
#       #pathID = rv$map_gsea
#       #myIP = rv$IP
#       #kegg_gsea_vector = rv$kegg_gsea_vector
#       
#       kegg_gsea_vector = rv$kegg_gsea_vector[names(rv$kegg_gsea_vector) %in% str_split(subset(ko_nodes %>% filter(map == rv$map_gsea))$nodes, pattern = "[,]")[[1]]]
#       
#       max_range = ceiling(max(abs(min(kegg_gsea_vector)), abs( max(kegg_gsea_vector)) ))
#       
#       outfile <- paste( tempFolder,"/ko", rv$map_gsea, ".", rv$IP, ".gsea.png",sep="")
#       setwd(tempFolder)
#       
#       pvout = pathview(gene.data = kegg_gsea_vector,
#                        gene.idtype = "KEGG",
#                        pathway.id = rv$map_gsea,
#                        species = "ko",
#                        bins = list(gene = 8),
#                        #limit = list(gene = c(0,1)),
#                        limit = list(gene = c(-max_range, max_range)),
#                        low = list(gene = rv$kegg_gsea_col_low  ),
#                        mid = list(gene = "yellow"),
#                        high = list(gene =  rv$kegg_gsea_col_high),
#                        na.col = "transparent",
#                        kegg.dir = file.path(wd, "kegg_data"),
#                        #out.suffix = suffix,
#                        kegg.native = T,
#                        same.layer = T,
#                        sign.pos = "bottomleft",
#                        new.signature = F,
#                        #key.pos = "topleft",
#                        out.suffix = paste0(rv$IP,".gsea"))
#       setwd(wd)
#       shinyjs::show("keggscape_gsea_conf")
#       
#       print(outfile)
#       return( list(src = outfile,
#                    contentType = 'image/png',
#                    width = "150%",
#                    height = "150%",
#                    alt = "This is alternate text") )
#       
#       
#     }else{
#       print("empty images KEGG")
#       html("text_KEGGscape_GSEA", "<strong>No enriched KEGG terms! Try again with different filters</strong>")
#       return (blank)
#       #empty_image()
#     }
#   }else{
#     return (blank)
#   }
# }, deleteFile = T)